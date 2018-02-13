#include "kernel2.hpp"
#include <iostream>
#include <x86intrin.h>


void poisson( const double xsources[], const double ysources[], const double wsources[], const int nsources, double phi[], const int M, const double xmesh[], const double ymesh[], const double theta)
{
    //data storage for tree
    const double *x, *y, *w;
    double *xsorted, *ysorted, *wsorted;
    double *expansions;

    x = xsources;
    y = ysources;
    w = wsources;
    const int n = nsources;
    //align memory for sorted data
    posix_memalign((void **)&xsorted, 32, sizeof(double) * n);
    posix_memalign((void **)&ysorted, 32, sizeof(double) * n);
    posix_memalign((void **)&wsorted, 32, sizeof(double) * n);

    const int k = 32*1;	// leaf capacity
    const int maxnodes = (n + k - 1) / k * 60; //maximum number of nodes
    Node* nodes;
    //align memory for tree (nodes and expansions)
    posix_memalign((void **)&nodes, 32, sizeof(Node) * maxnodes);
    posix_memalign((void **)&expansions, 32, sizeof(double) * 2 * ORDER * maxnodes);

    //build tree and expansions
    build(x, y, w, n, k, xsorted, ysorted, wsorted, nodes, expansions);


    const double thetasquared = theta*theta;

    //Timer tm; tm.start();

    //evaluate phi (in parallel) at mesh points using the tree
    #pragma omp parallel for schedule(static,1)
    for(int i = 0; i < M*M; ++i)
    {
        streamfunction(nodes, expansions, xsorted, ysorted, wsorted, thetasquared, phi + i, xmesh[i], ymesh[i]);
    }


    //double t = tm.elapsed();
    //printf("Evaluation took %.3f ms (%.3f us per target)\n", t*1e-6, t*1e-3 / NDST);

    //free memory
    free(xsorted);
    free(ysorted);
    free(wsorted);
    free(nodes);
    free(expansions);
}




void streamfunction(const Node* nodes, const double* expansions, const double *xdata, const double *ydata, const double *mdata,
                    const double thetasquared, double * const result, const double xt, const double yt)
{
    enum { BUFSIZE = 16 }; //buffer size (chosen arbitrarily, only multiple of 4?)

    int stack[LMAX * 3];

    int bufcount = 0;
    double rzs[BUFSIZE], izs[BUFSIZE], masses[BUFSIZE]; //space for temporarily saving parameters
    const double *rxps[BUFSIZE], *ixps[BUFSIZE];

    int stackentry = 0, maxentry = 0;

    stack[0] = 0;
    *result = 0;
    while(stackentry > -1)
    {
        //take first node from stack
        const int nodeid = stack[stackentry--];
        const Node * const node = nodes + nodeid;

        //assert if node is valid (either is leaf or has a child with higher nodeID)
        assert(nodeid < node->child_id || node->child_id == 0);

        //compute squared distance between target point and center of mass of current node
        const double r2 = (xt - node->xcom)*(xt - node->xcom) + (yt - node->ycom)*(yt - node->ycom);


        //if distance is high enough, proceed with e2p kernel
        if (node->r * node->r < thetasquared * r2)
        {
            //fill in data in buffer
            rzs[bufcount] = xt - node->xcom; //x-distance
            izs[bufcount] = yt - node->ycom; //y-distance
            masses[bufcount] = node->mass; //mass
            rxps[bufcount] = expansions + ORDER * (2 * nodeid + 0); //expansions
            ixps[bufcount] = expansions + ORDER * (2 * nodeid + 1);
            ++bufcount; //increase bufcount = number of objects in buffer

            if (bufcount == BUFSIZE) //if buffer is filled completely
            {
                bufcount = 0; //reset bufcount
                *result += e2p(rzs, izs, masses, rxps, ixps, BUFSIZE); //apply e2p kernel on buffered data

            }
        }
        //if distance between target and center of mass of node is too small
        else
        {
            //if node is a leaf
            if (node->child_id == 0)
            {
                //extract start and end of node
                const int s = node->part_start;
                const int e = node->part_end;

                //apply p2p kernel on particle in leaf
                *result += p2p(xdata + s, ydata + s, mdata + s, e - s, xt, yt);
            }
            //if node is not a leaf
            else
            {
                //put all 4 child nodes on stack for further search
                for(int c = 0; c < 4; ++c)
                    stack[++stackentry] = node->child_id + c;

                //adjust maxentry
                maxentry = std::max(maxentry, stackentry);
            }
        }
    }
    //if no further nodes: if there is data left in buffer, apply e2p kernel
    if (bufcount)
        *result += e2p(rzs, izs, masses, rxps, ixps, bufcount);

    *result *= (-1.); //negative sign in poisson equation for vorticity!

}



void velocityfield_openBC(const double phi[], const int M, const double dx, double ux[], double uy[])
{

    // finite difference for velocity field
    #pragma omp parallel for
    for(int i = 0; i < M; ++i)
    {
        for(int j = 0; j < M; ++j)
        {
            if(i == 0)  //lower border: contribution from outside the mesh is 0
            {
                ux[i*M+j] = (phi[(i+1)*M+j]-phi[(i)*M+j])/dx;
            }
            else if(i == M-1)  //upper border: contribution from outside the mesh is 0
            {
                ux[i*M+j] = (phi[(i)*M+j]-phi[(i-1)*M+j])/dx;
            }
            else  //inside mesh: contribution from both sides
            {
                ux[i*M+j] = (phi[(i+1)*M+j]-phi[(i-1)*M+j])/(2*dx);
            }
            if(j == 0)  //left border: contribution from outside the mesh is 0
            {
                uy[i*M+j] = -(phi[i*M+j+1]-phi[i*M])/(dx);
            }
            else if(j == M-1)  //right border: contribution from outside the mesh is 0
            {
                uy[i*M+j] = -(phi[i*M+j]-phi[i*M+j-1])/(dx);
            }

            else  //inside mehs: contribution from both sides
            {
                uy[i*M+j] = -(phi[i*M+j+1]-phi[i*M+j-1])/(2*dx);
            }
        }
    }

}


void vorticityfield_by_curl_openBC(const double ux[], const double uy[], const int M, const double dx, double vorticity[])
{

    #pragma omp parallel for
    for(int i = 0; i < M; ++i)
    {
        for(int j = 0; j < M; ++j)
        {
            double duyDx = 0.;
            double duxDy = 0.;
            if(i == 0)  //lower border: contribution from outside the mesh is 0
            {
                duxDy = (ux[(i+1)*M+j]-ux[(i)*M+j])/(dx);
            }
            else if(i == M-1)  //upper border: contribution from outside the mesh is 0
            {
                duxDy = (ux[(i)*M+j]-ux[(i-1)*M+j])/(dx);
            }
            else  //inside mesh: both contributions
            {
                duxDy = (ux[(i+1)*M+j]-ux[(i-1)*M+j])/(2*dx);
            }
            if(j == 0)  //left border: contribution from outside the mesh is 0
            {
                duyDx = (uy[i*M+j+1]-uy[i*M+j])/(dx);
            }
            else if(j == M-1)  //right border: contribution from outside the mesh is 0
            {
                duyDx = (uy[i*M+j]-uy[i*M+j-1])/(dx);
            }
            else  //inside mesh: both contributions
            {
                duyDx = (uy[i*M+j+1]-uy[i*M+j-1])/(2*dx);
            }
            vorticity[i*M+j] = duyDx-duxDy;
        }
    }
}


/*
For the convection step, we treat each grid node as a particle and update its location with a
first order explicit Euler step, where for example
The particles carry the vorticity computed at the grid node and constitute the sources for the first
step of the algorithm.
For accuracy, use the Lagrangian time step criterion computed from the velocity field according to LFC=0.1
*/

void convection(double xparticle[], double yparticle[], double wparticle[], const double ux[], const double uy[], const double xmesh[], const double ymesh[], double omega[], const int M, const double dt)
{
//cutting of veolcities at border to make sure poisson solver will work....
    for (int i=0; i<M; ++i)
    {
        for (int j=0; j<M; ++j)
        {
            if( i<0.1*M || i>0.9*M || j<0.1*M || j>0.9*M )
            {
                omega[i*M+j] = 0.0;
            }
        }
    }

    const double factor = 1.0/M/M/M_PI*0.5; //for poisson solver scaling
    double * const xptr = &xparticle[0];
    double * const yptr = &yparticle[0];
    const double * const xmsh = &xmesh[0];
    const double * const ymsh = &ymesh[0];
    const double * const uxptr = &ux[0];
    const double * const uyptr = &uy[0];
    const double * const omegamsh = &omega[0];
    double * const wptr = &wparticle[0];

    const __m256d d_t = _mm256_set1_pd(dt);
    const __m256d fp = _mm256_set1_pd(factor);


    for (int i=0; i<M*M/4; ++i)
    {
        const __m256d w = _mm256_loadu_pd(omegamsh+i*4);
        const __m256d tmpfact = _mm256_mul_pd(fp, w);

        const __m256d x = _mm256_loadu_pd(xmsh+i*4);
        const __m256d y = _mm256_loadu_pd(ymsh+i*4);

        const __m256d u_x = _mm256_loadu_pd(uxptr+i*4);
        const __m256d u_y = _mm256_loadu_pd(uyptr+i*4);

        const __m256d tmp1 = _mm256_mul_pd(d_t, u_x);
        const __m256d tmp2 = _mm256_mul_pd(d_t, u_y);

        const __m256d xp = _mm256_add_pd(tmp1, x);
        const __m256d yp = _mm256_add_pd(tmp2, y);


        _mm256_store_pd(xptr+i*4, xp);
        _mm256_store_pd(yptr+i*4, yp);

        _mm256_store_pd(wptr+i*4, tmpfact);


    }
//"translation"
    /*for (int i=0; i<M; ++i) {
        for (int j=0; j<M; ++j)  {
             xparticle[i*M+j] = xmesh[i*M+j] + dt* ux[i*M+j];
             yparticle[i*M+j] = ymesh[i*M+j] + dt* uy[i*M+j];
         wparticle[i*M+j] = omega[i*M+j];
             }}
    std::swap(wparticle, omega);*/

}

void particlexchangemethod_openBC(const double xmesh[], const double ymesh[], double wgrid[],
                                  double wtemp[], const int M, const double dx, const double visc, const double dt, Celllist cells)
{
//kernel:

    const double eps= 2.0*dx;
    const double PI16DX = 4.0/M_PI/M_PI*visc*dt; //(assuming that dx=dy), epsilon
    const int ncell = cells.getNumberOfCells();

    for ( int c=0; c<ncell*ncell; ++c)
    {
        double tmp=0.0;
        for(int n = 0; n < 9; ++n)
        {
            int index = cells.cellNeighbour_openBC(c,n);
            if(index > -1)  //index == -1 means that there is no such neighbour/that neighbour is empty ghost cell. These can be ignored.
            {
                for(auto i: cells.celllist_[c])
                {
                    for(auto p : cells.celllist_[index])
                    {
                        double distquad = ( ( xmesh[p]-xmesh[i] ) *( xmesh[p]-xmesh[i] ) + (ymesh[p]- ymesh[i]) *(ymesh[p]- ymesh[i]));

                        if(distquad >= 1e-12 && distquad< 25.0*eps*eps)
                        {
                            tmp += (wgrid[p]-wgrid[i]) /(distquad*distquad*distquad*distquad +1.0);
                        }
                    }
                    wtemp[i] = wgrid[i] + tmp*PI16DX;
                }
            }
        }

    }
    std::swap(wtemp, wgrid);

}


 // ################## OLD/UNUSED VERSIONS ####################

/*

//old version, do not use
void vorticityfield(const double xsources[], const double ysources[], const double wsources[], const int nsources, const int M, const double dx, const double xmesh[], const double ymesh[], double vorticity[]){
	
	#pragma omp parallel for
	for(int k = 0; k < M*M; ++k){
		double sum = 0.;
		for(int p = 0; p < nsources; p++){
			double lambdax = (xsources[p] - xmesh[k])/dx;
			double lambday = (ysources[p] - ymesh[k])/dx;
			sum += interpol(lambdax)*interpol(lambday)*wsources[p];
		}
		vorticity[k] = sum;
	}	
	
}
//old version, do not use
double interpol(double lambda){
	double abslambda = std::abs(lambda);
	if(abslambda >= 2.){
		return 0.;
	}
	else if(abslambda >= 1.){
		return 0.5*(2.-abslambda)*(2.-abslambda)*(1.-abslambda);
	}
	else{
		return 1. - 2.5*lambda*lambda - 1.5*abslambda*abslambda*abslambda;
	}
}



void particlexchangemethod(const double xmesh[], const double ymesh[], double wgrid[],
             double wtemp[], const int M, const double dx, const double visc, const double dt, Celllist cells)
{
//kernel:

const double eps= 2.0*dx;
const double PI16DX = 4.0/M_PI/M_PI*visc*dt; //(assuming that dx=dy), epsilon
//double rhotmp[M*M];

const int ncell = cells.getNumberOfCells();


	for ( int c=0; c<ncell*ncell; ++c){
		double tmp=0.0;
		for(int n = 0; n < 9; ++n){
			int index = cells.cellNeighbour(c,n);
			for(auto i: cells.celllist_[c]){
				for(auto p : cells.celllist_[index]){			
					double distquad = ( ( xmesh[p]-xmesh[i] ) *( xmesh[p]-xmesh[i] ) + (ymesh[p]- ymesh[i]) *(ymesh[p]- ymesh[i]));

					if(distquad >= 1e-12 && distquad< 25.0*eps*eps) {
						tmp += (wgrid[p]-wgrid[i]) /(distquad*distquad*distquad*distquad +1.0);

					}
				}
				wtemp[i] = wgrid[i] + tmp*PI16DX;		
			}
		}

	}	
	std::swap(wtemp, wgrid);

}

void velocityfield(const double phi[], const int M, const double dx, double ux[], double uy[]){
	
	// finite difference for velocity field
	#pragma omp parallel for
	for(int i = 0; i < M; ++i){
		for(int j = 0; j < M; ++j){
			if(i == 0){
				ux[i*M+j] = (phi[(i+1)*M+j]-phi[(M-1)*M+j])/(2*dx);
			}
			else if(i == M-1){
				ux[i*M+j] = (phi[0*M+j]-phi[(i-1)*M+j])/(2*dx);
			}
			else{
				ux[i*M+j] = (phi[(i+1)*M+j]-phi[(i-1)*M+j])/(2*dx);
			}
			if(j == 0){
				uy[i*M+j] = -(phi[i*M+j+1]-phi[i*M+M-1])/(2*dx);
			}
			else if(j == M-1){
				uy[i*M+j] = -(phi[i*M+0]-phi[i*M+j-1])/(2*dx);
			}
			else{			
				uy[i*M+j] = -(phi[i*M+j+1]-phi[i*M+j-1])/(2*dx);
			} 
		}
	}

	//what about only one for loot from 0 to M*M?	 --> slower!	
	
	#pragma omp parallel for
	for(int i = 0; i < M*M; ++i){
		if(i < M){
			ux[i] = -(phi[i+M]-phi[(M-1)*M+i])/dx;
		}
		else if(i > M*(M-1)){
			ux[i] = -(phi[i-M*(M-1)]-phi[i-M])/dx;
		}
		else{
			ux[i] = -(phi[i+M]-phi[i-M])/dx;
		}
		if(i % M == 0){
			uy[i] = (phi[i+1]-phi[i+M-1])/dx;
		}
		else if(i % M == M-1){
			uy[i] = (phi[i-M+1]-phi[i-1])/dx;
		}
		else{			
			uy[i] = (phi[i+1]-phi[i-1])/dx;
		} 
	}
	
		
}


void vorticityfield_by_curl(const double ux[], const double uy[], const int M, const double dx, double vorticity[]){
	
	#pragma omp parallel for
	for(int i = 0; i < M; ++i){
		for(int j = 0; j < M; ++j){
			double duyDx = 0.;
			double duxDy = 0.;
			if(i == 0){
				duxDy = (ux[(i+1)*M+j]-ux[(M-1)*M+j])/(2*dx);
			}
			else if(i == M-1){
				duxDy = (ux[0*M+j]-ux[i-1*M+j])/(2*dx);
			}
			else{
				duxDy = (ux[(i+1)*M+j]-ux[(i-1)*M+j])/(2*dx);
			}
			if(j == 0){
				duyDx = (uy[i*M+j+1]-uy[i*M+M-1])/(2*dx);
			}
			else if(j == M-1){
				duyDx = (uy[i*M+0]-uy[i*M+j-1])/(2*dx);
			}
			else{			
				duyDx = (uy[i*M+j+1]-uy[i*M+j-1])/(2*dx);
			} 
			vorticity[i*M+j] = duyDx-duxDy;
		}
	}
}

*/

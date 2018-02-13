#include <iostream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <cstdio>
#include <cassert>
#include <cstring>
#include<fstream>
#include <stdio.h>

/* 
generates INPUT DATA FILE for one vortex 
*/
void getmesh(const double xmin, const double xmax, const double ymin, const double ymax, double* xmesh, double* ymesh, const int M)
{
     double extx = xmax-xmin;
     double exty = ymax-ymin;

     //generate mesh from (xmin, ymin) to (xmax, ymax)
     for (int i=0; i<M; ++i)
     {
         for (int j=0; j<M; ++j)
         {
         xmesh[i*M+j] = xmin + j*extx/(M-1);
         ymesh[i*M+j] = ymin + i*exty/(M-1);
         }
     }
}


int main(int argc, char *argv[])
{
int n = 300;
const double T =0.220433;
const double coresize = 0.1;
const double strength = 5.0 ; //0.0;
const double xpos = 0.5;
const double ypos = 0.5;
const double xmin = 0.0;
const double xmax = 1.0;
const double nsources = n;
const double VISC=0.1;
const double coresize_t = std::sqrt(4.0*VISC*T+coresize*coresize);

double sum=0.0;

double *xmesh, *ymesh, *wmesh, *all, *all_analytical, *wmesh_analytical;
	posix_memalign((void **)&xmesh, 32, sizeof(double) * n*n);
	posix_memalign((void **)&ymesh, 32, sizeof(double) * n*n);
	posix_memalign((void **)&wmesh, 32, sizeof(double) * n*n);
	posix_memalign((void **)&wmesh_analytical, 32, sizeof(double) * n*n);
	posix_memalign((void **)&all, 32, sizeof(double) * (3*n*n+3));
	posix_memalign((void **)&all_analytical, 32, sizeof(double) * (3*n*n+3));

getmesh(0.,1.,0.,1., xmesh, ymesh, n);

std::string filenameh = "vortex_C01xS5xN300xT02.csv";
std::ofstream outfileh;
outfileh.open(filenameh);
//outfileh<< " X " << "\t"  << " Y "<< "\t" << " W "<< std::endl;

std::string filenamead = "AV01xC01xS5xN300xT0220344.csv.1";
std::ofstream outfilead;
outfilead.open(filenamead);
//outfilead<< " X " << "\t"  << " Y "<< "\t" << " W "<< std::endl;

std::FILE* outfilep;
outfilep = fopen("vortex_C01xS5xN300xT02.bin", "wb");

std::FILE* outfilea;
outfilea = fopen("analyticalV01xC01xS5xN300xT0220344.bin", "wb");


for (int p=0; p<n*n ; ++p)
{
double tmp = -( (xmesh[p]-xpos)*(xmesh[p]-xpos)/(coresize*coresize) + ( ymesh[p]-ypos)*(ymesh[p]-ypos)/(coresize*coresize) );
wmesh[p] = strength/M_PI/coresize/coresize*exp(tmp);
sum += wmesh[p];
wmesh[p] = wmesh[p]/n/n/M_PI*0.5; //scale for poisson solver

double tmp_20 =  -( (xmesh[p]-xpos)*(xmesh[p]-xpos)/(coresize_t*coresize_t) + ( ymesh[p]-ypos)*(ymesh[p]-ypos)/(coresize_t*coresize_t) );
wmesh_analytical[p] = strength/coresize_t/coresize_t/M_PI*exp(tmp_20);


if(xmesh[p]<0.1 || ymesh[p]<0.1 || xmesh[p]>0.9 || ymesh[p]>0.9)
{
wmesh[p] =0.0;
}
//wmesh[p] = strength/coresize/coresize/M_PI*exp(tmp)/n/n;
sum += wmesh[p];
outfileh<< xmesh[p] << " , " << ymesh[p] << " , " << wmesh[p] << std::endl;
outfilead<< xmesh[p] << " , " << ymesh[p] << " , " << wmesh_analytical[p] << std::endl;
}
outfileh.close();
outfilead.close();
 
all[0] = xmin;
all[1] = xmax;
all[2] = nsources;
for (int i=3; i< n*n+3 ; ++i) { all[i] = xmesh[i-3]; }
for (int i=n*n+3; i< n*n+n*n+3; ++i){ all[i] = ymesh[i-n*n-3]; }
for (int i=n*n+n*n+3; i< (3*n*n+3) ; ++i) { all[i] = wmesh[i-n*n-n*n-3];}

fwrite(&all[0], sizeof(double),(3*n*n+3), outfilep);
std::fclose(outfilep);

all_analytical[0] = xmin;
all_analytical[1] = xmax;
all_analytical[2] = nsources;
for (int i=3; i< n*n+3 ; ++i) { all_analytical[i] = xmesh[i-3]; }
for (int i=n*n+3; i< n*n+n*n+3; ++i){ all_analytical[i] = ymesh[i-n*n-3]; }
for (int i=n*n+n*n+3; i< (3*n*n+3) ; ++i) { all_analytical[i] = wmesh_analytical[i-n*n-n*n-3];}

fwrite(&all_analytical[0], sizeof(double),(3*n*n+3), outfilea);
std::fclose(outfilea);

std::cout<<" sum over all : " << sum/n/n<<std::endl;
//test:
/*
FILE * f = fopen("pcgaussdist33.bin", "r");

	
	double *xsources, *ysources, *wsources, *min, *max, *nn;
	posix_memalign((void **)&xsources, 32, sizeof(double) * nsources);
	posix_memalign((void **)&wsources, 32, sizeof(double) * nsources);
	posix_memalign((void **)&ysources, 32, sizeof(double) * nsources);
	posix_memalign((void **)&min, 32, sizeof(double));
	posix_memalign((void **)&max, 32, sizeof(double));
	posix_memalign((void **)&nn, 32, sizeof(double));


	fread(min, sizeof(double), 1, f);
	fread(max, sizeof(double), 1, f);
	fread(nn, sizeof(double), 1, f);
	fread(xsources, sizeof(double), nsources, f);
	fread(ysources, sizeof(double), nsources, f);
	fread(wsources, sizeof(double), nsources, f);

std::cout<< " xmin elem originial x : " << xmin << " writ : " << all[0] <<" read " << *min <<std::endl;
std::cout<< " n  elem originial x : " << n << " write&read : " << *nn <<std::endl;
std::cout<< " 2. elem originial x : " << xmesh[100] << " write&read : " << xsources[100] <<std::endl;
std::cout<< " 5. elem originial w : " << wmesh[100] << " write&read : " << wsources[100] <<std::endl;
std::cout<<" sum over all : " << sum/n/n <<std::endl;
*/
  return 0;}

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
generates INPUT DATA FILE for two vortices
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
const double coresize = 0.05;
const double strength1 = 5.0 ; 
const double xpos1 = 0.5-coresize*2.5;
const double ypos1 = 0.5;
const double strength2 = 5.0 ; 
const double xpos2 =  0.5+coresize*2.5;
const double ypos2 = 0.5;
const double T = 2.0;
const double xmin = 0.0;
const double xmax = 1.0;
const double nsources = n;

const double VISC=0.01;
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

std::string filenameh = "double_C005xS5xN300xT2.csv";
std::ofstream outfileh;
outfileh.open(filenameh);
//outfileh<< " X " << "\t"  << " Y "<< "\t" << " W "<< std::endl;

std::string filenamead = "doubleanalytical_V001xC005xS5xN300xT2.csv";
std::ofstream outfilead;
outfilead.open(filenamead);
//outfilead<< " X " << "\t"  << " Y "<< "\t" << " W "<< std::endl;

std::FILE* outfilep;
outfilep = fopen("double_C005xS5xN300xT2.bin", "wb");

std::FILE* outfilea;
outfilea = fopen("doubleanalytical_V001xC005xS5xN300xT2.bin", "wb");

for (int p=0; p<n*n ; ++p)
{
double tmp1 = -( (xmesh[p]-xpos1)*(xmesh[p]-xpos1)/(coresize*coresize) + ( ymesh[p]-ypos1)*(ymesh[p]-ypos1)/(coresize*coresize) );
double tmp2 = -( (xmesh[p]-xpos2)*(xmesh[p]-xpos2)/(coresize*coresize) + ( ymesh[p]-ypos2)*(ymesh[p]-ypos2)/(coresize*coresize) );

wmesh[p] = 1.0/coresize/coresize*M_PI*( strength1*exp(tmp1) + strength2*exp(tmp2) );
wmesh[p] = wmesh[p]/n/n/M_PI*0.5; //scale for poisson solver

double tmp20 =  -( (xmesh[p]-xpos2)*(xmesh[p]-xpos2)/(coresize_t*coresize_t) + ( ymesh[p]-ypos2)*(ymesh[p]-ypos2)/(coresize_t*coresize_t) );
double tmp10 =  -( (xmesh[p]-xpos1)*(xmesh[p]-xpos1)/(coresize_t*coresize_t) + ( ymesh[p]-ypos1)*(ymesh[p]-ypos1)/(coresize_t*coresize_t) );

wmesh_analytical[p] = 1.0/coresize/coresize*M_PI*( strength1*exp(tmp10) + strength2*exp(tmp20) );


if(xmesh[p]<0.1 || ymesh[p]<0.1 || xmesh[p]>0.9 || ymesh[p]>0.9)
{
wmesh[p] =0.0;
}

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

  return 0;}

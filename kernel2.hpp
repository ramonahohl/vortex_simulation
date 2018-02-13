/*
ETH ZÜRICH
HPCSE II : FS2016 : VORTEX PROJECT
Authors: Franziska Krummenacher (11-917-580) & Ramona Hohl (13-816-004)
*/

#ifndef KERNEL2
#define KERNEL2

#define LMAX 15
#include "kernels.h"	
#include "tree.h"
#include "Celllist.hpp"
#include <unistd.h>
#include <omp.h>

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <limits>


//Poisson Solver from Solution5
void poisson( const double xsources[], const double ysources[], const double wsources[], const int nsources, double phi[], const int M, const double xmesh[], const double ymesh[], const double theta);

//evaluate stream function
void streamfunction(const Node* nodes, const double* expansions, const double *xdata, const double *ydata, const double *mdata,
		const double thetasquared, double * const result, const double xt, const double yt);

//evaluate velocity field
void velocityfield_openBC(const double phi[], const int M, const double dx, double ux[], double uy[]);

//evaluate vorticity field
void vorticityfield_by_curl_openBC(const double ux[], const double uy[], const int M, const double dx, double vorticity[]);

//Convection step
void convection(double xparticle[], double yparticle[], double wparticle[], const double ux[], const double uy[], const double xmesh[], const double ymesh[], double omega[], const int M, const double dt);

// Particle Exchange Method
void particlexchangemethod_openBC(const double xmesh[], const double ymesh[], double wgrid[],
             double wtemp[], const int M, const double dx, const double visc, const double dt, Celllist cells);


/*
double interpol(double lambda);
void velocityfield(const double phi[], const int M, const double dx, double ux[], double uy[]);
void particlexchangemethod(const double xmesh[], const double ymesh[], double wgrid[],
             double wtemp[], const int M, const double dx, const double visc, const double dt, Celllist cells);
void vorticityfield(const double xsources[], const double ysources[], const double wsources[], const int nsources, const int M, const double dx, const double xmesh[], const double ymesh[], double vorticity[]);
*/          

#endif

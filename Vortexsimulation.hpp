/*
ETH ZÜRICH
HPCSE II : FS2016 : VORTEX PROJECT
Authors: Franziska Krummenacher (11-917-580) & Ramona Hohl (13-816-004)
*/

#ifndef VORTEXSIM
#define VORTEXSIM

#include "kernel2.hpp"
#include "kernels.h"
#include "timer.hpp"
#include "tree.h"
#include "Celllist.hpp"

#include <cstring>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <unistd.h>
#include <cstdio>
#include <string>



class Vortexsimulation{
	public:
	Vortexsimulation(const char filename[]);
	~Vortexsimulation();
	
	void initialize(const char filename[]);
	void getMesh();
	void advance();
	double getTime();
	double getdt();
	void advanceandoutput(std::string );
	void analyze();
	
	private:
	//const values
	const double LCFL = 0.1;
	const double VISC = 0.1;
	const double theta = 0.5;
	//data
	double xmin_, xmax_;
	int nsources_;
	double *xsources_;
	double *ysources_;
	double *wsources_;
	int M_;
	double *xmesh_;
	double *ymesh_;
	double dx_;
	double *phi_;
	double *ux_;
	double *uy_;
	double *vorticity_, *vorticity_temp_;
	double dt_; 
	double T_;
	Celllist celllist_;
	double cutoff_;
	int outnumber_;
};

#endif

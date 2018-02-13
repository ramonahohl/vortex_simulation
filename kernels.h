/*
 *  kernels.h
 *  ex5 solution
 *
 *  Created by Dmitry Alexeev on May 30, 2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */


#pragma once

void p2e(
	const  double xsources[],
	const  double ysources[],
	const  double qsources[],
	const  int nsources,
	const  double xcom,
	const  double ycom,
        double rexpansions[],
        double iexpansions[]);


double e2p(
       const double rzs[],
       const double izs[],
       const double masses[],
       const double * const  rxps[],
       const double * const  ixps[],
       const int ndst);


double p2p(
    const double *  _xsrc,
    const double *  _ysrc,
    const double *  _vsrc,
    const int nsources,
    const double xt,
    const double yt);



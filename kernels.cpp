#include <stdio.h>
#include <math.h>
#include <stdlib.h>




#define EPS (10 * __DBL_EPSILON__)



void p2e(
	const  double xsources[],
	const  double ysources[],
	const  double qsources[],
	const  int nsources,
	const  double xcom,
	const  double ycom,
  double *__restrict rexpansions,
  double *__restrict iexpansions  )
{
	double 
		rxp_0 = 0, ixp_0 = 0,
		rxp_1 = 0, ixp_1 = 0,
		rxp_2 = 0, ixp_2 = 0,
		rxp_3 = 0, ixp_3 = 0,
		rxp_4 = 0, ixp_4 = 0,
		rxp_5 = 0, ixp_5 = 0,
		rxp_6 = 0, ixp_6 = 0,
		rxp_7 = 0, ixp_7 = 0,
		rxp_8 = 0, ixp_8 = 0,
		rxp_9 = 0, ixp_9 = 0,
		rxp_10 = 0, ixp_10 = 0,
		rxp_11 = 0, ixp_11 = 0;

        for(int i = 0; i < nsources; ++i)
        {
		const double rprod_0 = xsources[i] - xcom;
		const double iprod_0 = ysources[i] - ycom;

		const double src = qsources[i];
                double rtmp = rprod_0 * src;
                double itmp = iprod_0 * src;

                rxp_0 -= rtmp;
                ixp_0 -= itmp;

                double rprod = rprod_0, iprod = iprod_0;

                
                rtmp = rprod * rprod_0 - iprod * iprod_0;
                itmp = rprod * iprod_0 + iprod * rprod_0;

                const double term_1 = src * (double)(1 / 2.);

                rprod = rtmp;
                iprod = itmp;

                rtmp = rprod * term_1;
                itmp = iprod * term_1;

                rxp_1 -= rtmp;
                ixp_1 -= itmp;
                
                rtmp = rprod * rprod_0 - iprod * iprod_0;
                itmp = rprod * iprod_0 + iprod * rprod_0;

                const double term_2 = src * (double)(1 / 3.);

                rprod = rtmp;
                iprod = itmp;

                rtmp = rprod * term_2;
                itmp = iprod * term_2;

                rxp_2 -= rtmp;
                ixp_2 -= itmp;
                
                rtmp = rprod * rprod_0 - iprod * iprod_0;
                itmp = rprod * iprod_0 + iprod * rprod_0;

                const double term_3 = src * (double)(1 / 4.);

                rprod = rtmp;
                iprod = itmp;

                rtmp = rprod * term_3;
                itmp = iprod * term_3;

                rxp_3 -= rtmp;
                ixp_3 -= itmp;
                
                rtmp = rprod * rprod_0 - iprod * iprod_0;
                itmp = rprod * iprod_0 + iprod * rprod_0;

                const double term_4 = src * (double)(1 / 5.);

                rprod = rtmp;
                iprod = itmp;

                rtmp = rprod * term_4;
                itmp = iprod * term_4;

                rxp_4 -= rtmp;
                ixp_4 -= itmp;
                
                rtmp = rprod * rprod_0 - iprod * iprod_0;
                itmp = rprod * iprod_0 + iprod * rprod_0;

                const double term_5 = src * (double)(1 / 6.);

                rprod = rtmp;
                iprod = itmp;

                rtmp = rprod * term_5;
                itmp = iprod * term_5;

                rxp_5 -= rtmp;
                ixp_5 -= itmp;
                
                rtmp = rprod * rprod_0 - iprod * iprod_0;
                itmp = rprod * iprod_0 + iprod * rprod_0;

                const double term_6 = src * (double)(1 / 7.);

                rprod = rtmp;
                iprod = itmp;

                rtmp = rprod * term_6;
                itmp = iprod * term_6;

                rxp_6 -= rtmp;
                ixp_6 -= itmp;
                
                rtmp = rprod * rprod_0 - iprod * iprod_0;
                itmp = rprod * iprod_0 + iprod * rprod_0;

                const double term_7 = src * (double)(1 / 8.);

                rprod = rtmp;
                iprod = itmp;

                rtmp = rprod * term_7;
                itmp = iprod * term_7;

                rxp_7 -= rtmp;
                ixp_7 -= itmp;
                
                rtmp = rprod * rprod_0 - iprod * iprod_0;
                itmp = rprod * iprod_0 + iprod * rprod_0;

                const double term_8 = src * (double)(1 / 9.);

                rprod = rtmp;
                iprod = itmp;

                rtmp = rprod * term_8;
                itmp = iprod * term_8;

                rxp_8 -= rtmp;
                ixp_8 -= itmp;
                
                rtmp = rprod * rprod_0 - iprod * iprod_0;
                itmp = rprod * iprod_0 + iprod * rprod_0;

                const double term_9 = src * (double)(1 / 10.);

                rprod = rtmp;
                iprod = itmp;

                rtmp = rprod * term_9;
                itmp = iprod * term_9;

                rxp_9 -= rtmp;
                ixp_9 -= itmp;
                
                rtmp = rprod * rprod_0 - iprod * iprod_0;
                itmp = rprod * iprod_0 + iprod * rprod_0;

                const double term_10 = src * (double)(1 / 11.);

                rprod = rtmp;
                iprod = itmp;

                rtmp = rprod * term_10;
                itmp = iprod * term_10;

                rxp_10 -= rtmp;
                ixp_10 -= itmp;
                
                rtmp = rprod * rprod_0 - iprod * iprod_0;
                itmp = rprod * iprod_0 + iprod * rprod_0;

                const double term_11 = src * (double)(1 / 12.);

                rprod = rtmp;
                iprod = itmp;

                rtmp = rprod * term_11;
                itmp = iprod * term_11;

                rxp_11 -= rtmp;
                ixp_11 -= itmp;
                
        }
	
         double rsum_0 = rxp_0;
         double isum_0 = ixp_0;
         double rsum_1 = rxp_1;
         double isum_1 = ixp_1;
         double rsum_2 = rxp_2;
         double isum_2 = ixp_2;
         double rsum_3 = rxp_3;
         double isum_3 = ixp_3;
         double rsum_4 = rxp_4;
         double isum_4 = ixp_4;
         double rsum_5 = rxp_5;
         double isum_5 = ixp_5;
         double rsum_6 = rxp_6;
         double isum_6 = ixp_6;
         double rsum_7 = rxp_7;
         double isum_7 = ixp_7;
         double rsum_8 = rxp_8;
         double isum_8 = ixp_8;
         double rsum_9 = rxp_9;
         double isum_9 = ixp_9;
         double rsum_10 = rxp_10;
         double isum_10 = ixp_10;
         double rsum_11 = rxp_11;
         double isum_11 = ixp_11;

        
        rexpansions[0] = rsum_0;
        iexpansions[0] = isum_0;
        
        rexpansions[1] = rsum_1;
        iexpansions[1] = isum_1;
        
        rexpansions[2] = rsum_2;
        iexpansions[2] = isum_2;
        
        rexpansions[3] = rsum_3;
        iexpansions[3] = isum_3;
        
        rexpansions[4] = rsum_4;
        iexpansions[4] = isum_4;
        
        rexpansions[5] = rsum_5;
        iexpansions[5] = isum_5;
        
        rexpansions[6] = rsum_6;
        iexpansions[6] = isum_6;
        
        rexpansions[7] = rsum_7;
        iexpansions[7] = isum_7;
        
        rexpansions[8] = rsum_8;
        iexpansions[8] = isum_8;
        
        rexpansions[9] = rsum_9;
        iexpansions[9] = isum_9;
        
        rexpansions[10] = rsum_10;
        iexpansions[10] = isum_10;
        
        rexpansions[11] = rsum_11;
        iexpansions[11] = isum_11;
        
}

double e2p(
       const double rzs[],
       const double izs[],
       const double masses[],
       const double * const  rxps[],
       const double * const  ixps[],
       const int ndst)
  {
	double result = 0;

	for(int i = 0;  i < ndst ; ++i)
	{
	   const double rz = rzs[i];
	   const double iz = izs[i];
	   const double r2 = rz * rz + iz * iz;
	   const double rinvz_1 = rz / r2;
	   const double iinvz_1 = -iz / r2;

	   
     	   const double rinvz_2 = rinvz_1 * rinvz_1 - iinvz_1 * iinvz_1;
     	   const double iinvz_2 = rinvz_1 * iinvz_1 + iinvz_1 * rinvz_1;
     	   const double rinvz_3 = rinvz_2 * rinvz_1 - iinvz_2 * iinvz_1;
     	   const double iinvz_3 = rinvz_2 * iinvz_1 + iinvz_2 * rinvz_1;
     	   const double rinvz_4 = rinvz_3 * rinvz_1 - iinvz_3 * iinvz_1;
     	   const double iinvz_4 = rinvz_3 * iinvz_1 + iinvz_3 * rinvz_1;
     	   const double rinvz_5 = rinvz_4 * rinvz_1 - iinvz_4 * iinvz_1;
     	   const double iinvz_5 = rinvz_4 * iinvz_1 + iinvz_4 * rinvz_1;
     	   const double rinvz_6 = rinvz_5 * rinvz_1 - iinvz_5 * iinvz_1;
     	   const double iinvz_6 = rinvz_5 * iinvz_1 + iinvz_5 * rinvz_1;
     	   const double rinvz_7 = rinvz_6 * rinvz_1 - iinvz_6 * iinvz_1;
     	   const double iinvz_7 = rinvz_6 * iinvz_1 + iinvz_6 * rinvz_1;
     	   const double rinvz_8 = rinvz_7 * rinvz_1 - iinvz_7 * iinvz_1;
     	   const double iinvz_8 = rinvz_7 * iinvz_1 + iinvz_7 * rinvz_1;
     	   const double rinvz_9 = rinvz_8 * rinvz_1 - iinvz_8 * iinvz_1;
     	   const double iinvz_9 = rinvz_8 * iinvz_1 + iinvz_8 * rinvz_1;
     	   const double rinvz_10 = rinvz_9 * rinvz_1 - iinvz_9 * iinvz_1;
     	   const double iinvz_10 = rinvz_9 * iinvz_1 + iinvz_9 * rinvz_1;
     	   const double rinvz_11 = rinvz_10 * rinvz_1 - iinvz_10 * iinvz_1;
     	   const double iinvz_11 = rinvz_10 * iinvz_1 + iinvz_10 * rinvz_1;
     	   const double rinvz_12 = rinvz_11 * rinvz_1 - iinvz_11 * iinvz_1;
     	   const double iinvz_12 = rinvz_11 * iinvz_1 + iinvz_11 * rinvz_1;

	   
	   double rsum_0 = rxps[i][0] * rinvz_1 - ixps[i][0] * iinvz_1;
	   double rsum_1 = rxps[i][1] * rinvz_2 - ixps[i][1] * iinvz_2;
	   double rsum_2 = rxps[i][2] * rinvz_3 - ixps[i][2] * iinvz_3;
	   double rsum_3 = rxps[i][3] * rinvz_4 - ixps[i][3] * iinvz_4;
	   double rsum_4 = rxps[i][4] * rinvz_5 - ixps[i][4] * iinvz_5;
	   double rsum_5 = rxps[i][5] * rinvz_6 - ixps[i][5] * iinvz_6;
	   double rsum_6 = rxps[i][6] * rinvz_7 - ixps[i][6] * iinvz_7;
	   double rsum_7 = rxps[i][7] * rinvz_8 - ixps[i][7] * iinvz_8;
	   double rsum_8 = rxps[i][8] * rinvz_9 - ixps[i][8] * iinvz_9;
	   double rsum_9 = rxps[i][9] * rinvz_10 - ixps[i][9] * iinvz_10;
	   double rsum_10 = rxps[i][10] * rinvz_11 - ixps[i][10] * iinvz_11;
	   double rsum_11 = rxps[i][11] * rinvz_12 - ixps[i][11] * iinvz_12;
 	   
rsum_0 += rsum_1; 
rsum_2 += rsum_3; 
rsum_4 += rsum_5; 
rsum_6 += rsum_7; 
rsum_8 += rsum_9; 
rsum_10 += rsum_11;  

rsum_0 += rsum_2; 
rsum_4 += rsum_6; 
rsum_8 += rsum_10;  

rsum_0 += rsum_4; 
  

rsum_0 += rsum_8;  

	   result += masses[i] * log(r2) / 2 + rsum_0;
	 }

	 return result;
  }

double p2p(
    const double *  _xsrc,
    const double *  _ysrc,
    const double *  _vsrc,
    const int nsources,
    const double xt,
    const double yt)
   {
   const double eps = EPS;

   double 
    s_0 = 0
   , s_1 = 0
   , s_2 = 0
   , s_3 = 0 ;
    const int nnice = 16 * (nsources / 16);

      for( int i = 0; i < nnice; i += 16)
      {
	const  double *  const xsrc = _xsrc + i;
      	const  double *  const ysrc = _ysrc + i;
      	const  double *  const vsrc = _vsrc + i;

	
	#pragma unroll
	for(int programIndex = 0; programIndex < 4; ++programIndex)
	{
	const double xr_0 = xt - xsrc[programIndex + 0];
	const double yr_0 = yt - ysrc[programIndex + 0];
	const double r2_0 = xr_0 * xr_0 + yr_0 * yr_0;
	const double factor_0 = fabs(r2_0) > EPS;
	s_0 += factor_0 * log(r2_0 + EPS) * vsrc[programIndex + 0];
	}
	
	#pragma unroll
	for(int programIndex = 0; programIndex < 4; ++programIndex)
	{
	const double xr_1 = xt - xsrc[programIndex + 4];
	const double yr_1 = yt - ysrc[programIndex + 4];
	const double r2_1 = xr_1 * xr_1 + yr_1 * yr_1;
	const double factor_1 = fabs(r2_1) > EPS;
	s_1 += factor_1 * log(r2_1 + EPS) * vsrc[programIndex + 4];
	}
	
	#pragma unroll
	for(int programIndex = 0; programIndex < 4; ++programIndex)
	{
	const double xr_2 = xt - xsrc[programIndex + 8];
	const double yr_2 = yt - ysrc[programIndex + 8];
	const double r2_2 = xr_2 * xr_2 + yr_2 * yr_2;
	const double factor_2 = fabs(r2_2) > EPS;
	s_2 += factor_2 * log(r2_2 + EPS) * vsrc[programIndex + 8];
	}
	
	#pragma unroll
	for(int programIndex = 0; programIndex < 4; ++programIndex)
	{
	const double xr_3 = xt - xsrc[programIndex + 12];
	const double yr_3 = yt - ysrc[programIndex + 12];
	const double r2_3 = xr_3 * xr_3 + yr_3 * yr_3;
	const double factor_3 = fabs(r2_3) > EPS;
	s_3 += factor_3 * log(r2_3 + EPS) * vsrc[programIndex + 12];
	}
	
    	}

	
s_0 += s_1; 
s_2 += s_3;  

s_0 += s_2;  

	for(int i = nnice; i < nsources; ++i)
    	{
	    const double xr = xt - _xsrc[i];
      	    const double yr = yt - _ysrc[i];
      	    const double r2 = xr * xr + yr * yr;
      	    const double f  = fabs(r2) > eps;

      	    s_0 += f * log(r2 + eps) * _vsrc[i];
    	}

	return s_0 * 0.5;
   }


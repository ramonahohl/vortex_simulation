#include <stdio.h>
#include <math.h>
#include <stdlib.h>
include(unroll.m4)

#define EPS (10 * __DBL_EPSILON__)
define(NACC, 4)
define(GSIZE, 4)

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
	double LUNROLL(n, 0, eval(ORDER - 1),`ifelse(n,0,,`,')
		TMP(rxp, n) = 0, TMP(ixp, n) = 0');

        for(int i = 0; i < nsources; ++i)
        {
		const double rprod_0 = xsources[i] - xcom;
		const double iprod_0 = ysources[i] - ycom;

		const double src = qsources[i];
                double rtmp = rprod_0 * src;
                double itmp = iprod_0 * src;

                TMP(rxp, 0) -= rtmp;
                TMP(ixp, 0) -= itmp;

                double rprod = rprod_0, iprod = iprod_0;

                LUNROLL(n, 1, eval(ORDER - 1),`
                rtmp = rprod * TMP(rprod, 0) - iprod * TMP(iprod, 0);
                itmp = rprod * TMP(iprod, 0) + iprod * TMP(rprod, 0);

                const double TMP(term, n) = src * (double)(1 / eval(n+1).);

                rprod = rtmp;
                iprod = itmp;

                rtmp = rprod * TMP(term, n);
                itmp = iprod * TMP(term, n);

                TMP(rxp, n) -= rtmp;
                TMP(ixp, n) -= itmp;
                ')
        }
	LUNROLL(i, 0, eval(ORDER - 1), `
         double TMP(rsum, i) = TMP(rxp, i);
         double TMP(isum, i) = TMP(ixp, i);')

        LUNROLL(i, 0, eval(ORDER - 1), `
        rexpansions[i] = TMP(rsum, i);
        iexpansions[i] = TMP(isum, i);
        ')
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

	   LUNROLL(j, 2, ORDER, `
     	   const double TMP(rinvz, j) = TMP(rinvz, eval(j - 1)) * rinvz_1 - TMP(iinvz, eval(j - 1)) * iinvz_1;
     	   const double TMP(iinvz, j) = TMP(rinvz, eval(j - 1)) * iinvz_1 + TMP(iinvz, eval(j - 1)) * rinvz_1;')

	   LUNROLL(j, 1, ORDER, `
	   double TMP(rsum, eval(j - 1)) = rxps[i][eval(j - 1)] * TMP(rinvz, j) - ixps[i][eval(j - 1)] * TMP(iinvz, j);')
 	   REDUCE(`+=', LUNROLL(i, 0, eval(ORDER - 1),`ifelse(i,0,,`,')TMP(rsum, i)'))

	   result += masses[i] * log(r2) / 2 + TMP(rsum, 0);
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

   double LUNROLL(`i', 0, eval(NACC - 1), `
   ifelse(i,0,,`,') TMP(s,i) = 0') ;
    const int nnice = eval(4 * NACC) * (nsources / eval(4 * NACC));

      for( int i = 0; i < nnice; i += eval(4 * NACC))
      {
	const  double *  const xsrc = _xsrc + i;
      	const  double *  const ysrc = _ysrc + i;
      	const  double *  const vsrc = _vsrc + i;

	LUNROLL(j, 0, eval(NACC - 1), `
	#pragma unroll
	for(int programIndex = 0; programIndex < 4; ++programIndex)
	{
	const double TMP(xr, j) = xt - xsrc[programIndex + eval(j * GSIZE)];
	const double TMP(yr, j) = yt - ysrc[programIndex + eval(j * GSIZE)];
	const double TMP(r2, j) = TMP(xr, j) * TMP(xr, j) + TMP(yr, j) * TMP(yr, j);
	const double TMP(factor, j) = fabs(TMP(r2, j)) > EPS;
	TMP(s, j) += TMP(factor, j) * log(TMP(r2, j) + EPS) * vsrc[programIndex + eval(j * GSIZE)];
	}
	')
    	}

	REDUCE(`+=', LUNROLL(`i', 0, eval(NACC - 1), `ifelse(i,0,,`,') TMP(s,i)'))

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


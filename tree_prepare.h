#include <omp.h>

#include <cstdio>

#if  !defined(__INTEL_COMPILER) && !defined(__clang__)
#include <parallel/algorithm>
#else
#include <algorithm>
#endif

#include <cmath>

#define LMAX 15

void minmax_vec(const double xsrc[], const double ysrc[], const int nsources, double xmin_xmax_ymin_ymax[])
{
	double lxmi = 1e13, lymi = 1e13, lxma = -1e13, lyma = -1e13;

	for(int i = 0; i < nsources; ++i)
	{
		const double xval = xsrc[i];
		const double yval = ysrc[i];

		lxmi = fmin(lxmi, xval);
		lxma = fmax(lxma, xval);

		lymi = fmin(lymi, yval);
		lyma = fmax(lyma, yval);
	}

	xmin_xmax_ymin_ymax[0] = lxmi;
	xmin_xmax_ymin_ymax[1] = lxma;
	xmin_xmax_ymin_ymax[2] = lymi;
	xmin_xmax_ymin_ymax[3] = lyma;
}

void extent(const int N, const double* const x, const double* const y,
			double& xmin, double& ymin, double& ext)
{
	double ext0, ext1;
	static const int chunksize = 1024 * 4;
	
	{
		const int nthreads = omp_get_max_threads();

		double xpartials[2][nthreads], ypartials[2][nthreads];

		#pragma omp parallel
		{
			const int tid = omp_get_thread_num();

			double lxmi = 1e13, lymi = 1e13, lxma = -1e13, lyma = -1e13;

			const int start = tid * chunksize;
			const int step = nthreads * chunksize;

			for(int i = start; i < N; i += step)
			{
				double result[4];

				minmax_vec(x + i, y + i, std::min((int)chunksize, N - i), result);

				lxmi = std::min(lxmi, result[0]);
				lxma = std::max(lxma, result[1]);
				lymi = std::min(lymi, result[2]);
				lyma = std::max(lyma, result[3]);
			}

			xpartials[0][tid] = lxmi;
			xpartials[1][tid] = lxma;
			ypartials[0][tid] = lymi;
			ypartials[1][tid] = lyma;
		}

		xmin = *std::min_element(xpartials[0], xpartials[0] + nthreads);
		ymin = *std::min_element(ypartials[0], ypartials[0] + nthreads);

		ext0 = (*std::max_element(xpartials[1], xpartials[1] + nthreads) - xmin);
		ext1 = (*std::max_element(ypartials[1], ypartials[1] + nthreads) - ymin);
	}

	// For numerical reasons, shift the domain boundaries out by epsilon
	const double eps = 10000 * std::numeric_limits<double>::epsilon();
	ext = std::max(ext0, ext1) * (1 + 2 * eps);
	xmin -= eps * ext;
	ymin -= eps * ext;
}

void morton(const int N, const double* const x, const double* const y,
			const double xmin, const double ymin, const double ext, int* index)
{
	#pragma omp parallel for
	for(int i = 0; i < N; ++i)
	{
		int xid = floor((x[i] - xmin) / ext * (1 << LMAX));
		int yid = floor((y[i] - ymin) / ext * (1 << LMAX));

		xid = (xid | (xid << 8)) & 0x00FF00FF;
		xid = (xid | (xid << 4)) & 0x0F0F0F0F;
		xid = (xid | (xid << 2)) & 0x33333333;
		xid = (xid | (xid << 1)) & 0x55555555;

		yid = (yid | (yid << 8)) & 0x00FF00FF;
		yid = (yid | (yid << 4)) & 0x0F0F0F0F;
		yid = (yid | (yid << 2)) & 0x33333333;
		yid = (yid | (yid << 1)) & 0x55555555;

		index[i] = xid | (yid << 1);
	}
}

void sort(const int N, int* index, int* keys)
{
	std::pair<int, int> * kv = NULL;
	posix_memalign((void **)&kv, 32, sizeof(*kv) * N);

	#pragma omp parallel for
	for(int i = 0; i < N; ++i)
	{
		kv[i].first = index[i];
		kv[i].second = keys[i];
	}

#if  !defined(__INTEL_COMPILER) && !defined(__clang__)
	__gnu_parallel::sort(kv, kv + N);
#else
	std::sort(kv, kv + N);
#endif


	#pragma omp parallel for
	for(int i = 0; i < N; ++i)
	{
		index[i] = kv[i].first;
		keys[i] = kv[i].second;
	}

	free(kv);
}

void reorder(const int N, const int* const keys, const double* const x, const double* const y, const double* const m, double* xsorted, double* ysorted, double *msorted)
{
	#pragma omp parallel for
	for(int i = 0; i < N; ++i)
	{
		const int entry = keys[i];
		
		xsorted[i] = x[entry];
		ysorted[i] = y[entry];
		msorted[i] = m[entry];
	}
}

inline void node_setup(const double xsources[], const double ysources[], const double msources[], const int nsources,
		double& mass, double& xcom, double& ycom, double& radius)
{
	mass = 0;
	double weight = 0, xsum = 0, ysum = 0;

	for(int i = 0; i < nsources; ++i)
	{
		const double x = xsources[i];
		const double y = ysources[i];
		const double m = msources[i];
		const double w = fabs(m);

		mass += m;
		weight += w;
		xsum += x * w;
		ysum += y * w;
	}

	xcom = weight ? xsum / weight : 0;
	ycom = weight ? ysum / weight : 0;

	double r2 = 0;
	for(int i = 0; i < nsources; ++i)
	{
		const double xr = xsources[i] - xcom;
		const double yr = ysources[i] - ycom;

		r2 = fmax(r2, xr * xr + yr * yr);
	}

	radius = sqrt(r2);
}

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

int lower_bound_vec(int s, int e, const int val, const int keys[])
{
	int c = e - s;

	if (keys[s] >= val)
		return s;

	if (keys[e - 1] < val)
		return e;

	while (c)
	{
		const int s0 = s, e0 = e;

		const double h = (e - s) * 1.f / 8;

		for(int programIndex = 0; programIndex < 8; ++programIndex)
		{
			//int candidate_s = s0, candidate_e = e0;
			const int i = MIN(e0 - 1, (int)(s0 + programIndex * h + 0.499999f));

			const bool isless = keys[i] < val;
			const int candidate_s = isless ? i : s0;
			const int candidate_e = isless ? e0 : i;

			s = MAX(s, candidate_s);
			e = MIN(e, candidate_e);
		}
	
		c = MIN(c / 8, e - s);
	}

	return s + 1;
}

int upper_bound_vec(int s, int e, const int val, const int keys[])
{
	int c = e - s;

	if (keys[s] > val)
		return s;

	if (keys[e - 1] <= val)
		return e;

	while (c)
	{
		const int s0 = s, e0 = e;

		const double h = (e - s) * 1.f / 8;

		for(int programIndex = 0; programIndex < 8; ++programIndex)
		{
			//int candidate_s = s0, candidate_e = e0;
			const int i = MIN(e0 - 1, (int)(s0 + programIndex * h + 0.499999f));
	
			const bool isless = keys[i] <= val;
			const int candidate_s = isless ? i : s0;
			const int candidate_e = isless ? e0 : i;

			s = MAX(s, candidate_s);
			e = MIN(e, candidate_e);
		}
	
		c = MIN(c / 8, e - s);
	}

	return s + 1;
}


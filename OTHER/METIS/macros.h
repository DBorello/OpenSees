/*
 * macros.h
 *
 * This file contains macros used in multilevel
 *
 * Started 9/25/94
 * George
 *
 * $Id: macros.h,v 1.2 2008-01-28 19:24:41 fmk Exp $
 *
 */

/*************************************************************************
* The following macro returns a random number in the specified range
**************************************************************************/
/*
#define RandomInRange(u) ((rand()>>2)%(u))
*/

#ifdef _WIN32

#define drand48(x) ((double)(rand()/RAND_MAX))
#define RandomInRange(u) ((int)(drand48()*((double)(u))))

#else

#define RandomInRange(u) ((int)(drand48()*((double)(u))))

#endif

#define HTVALUE(k, n) ((k)%(n))

#define FEwgtVtx(v) (((v)->cewgt + (v)->ewgtsum)*((v)->cewgt + (v)->ewgtsum))

#define SWAP(a, b, tmp) \
                        do {(tmp) = (a); (a) = (b); (b) = (tmp);} while(0)

#define INC_DEC(a, b, val) \
                        do {(a) += (val); (b) -= (val);} while(0)

#define TIMELVL(x) x;


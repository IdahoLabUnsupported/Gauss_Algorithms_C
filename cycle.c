/*
 * Copyright 2017 Battelle Energy Alliance
 */

/*
 * cycle.c - the part of the Gauss fit algorithm that calls lmder()
 */

#include <stdio.h>		/* for printf() */
#include <math.h>		/* for fabs(), sqrt() */
#include <float.h>		/* for DBL_EPSILON & DBL_MIN */
#include <stdlib.h>		/* for calloc() and NULL */
#include "GaussLib.h"
#include "GLprivate.h"		/* for lmder(), calc_resid(), calc_jacobian() */

/*
 * value for lmder_info to indicate error
 */

#define	CY_TERMINATE	-1

/*
 * definitions for dpmpar_()
 */

#ifdef DBL_EPSILON
#define	CY_DPMPAR_ONE	DBL_EPSILON
#else
#define	CY_DPMPAR_ONE	2.2204460492503131E-16
#endif

#ifdef DBL_MIN
#define CY_DPMPAR_TWO	DBL_MIN
#else
#define CY_DPMPAR_TWO	2.2250738585072014E-308
#endif

/*
 * global variables for use in fcn_()
 */

static GLChanRange	CY_FCN_chanrange;
static GLSpectrum	CY_FCN_spectrum;
static GPFitVary	*CY_FCN_fitvary;
static GLFitInfo	*CY_FCN_fitinfo;

/* prototypes of private utilities */

static void CY_covary(int n, double *r, int ldr, int *ipvt, double tol,
                      double *wa);
static GPCycleType CY_cycle_nosplit(const GLSpectrum *spectrum,
                                    GPFitVary *fitvary, double xtol,
                                    double ftol, GLFitRecord *fit,
                                    GLCurve *curve);
static GPCycleType CY_cycle_split(const GLSpectrum *spectrum,
                                  GPFitVary *fitvary, double xtol,
                                  double ftol, GLFitRecord *fit,
                                  GLCurve *curve);
static void CY_fcn_(int *m, int *n, double *x, double *fvec, double *fjac,
                    int *ldfjac, int *iflag);
static void CY_init_fitvary(GLPkwdMode pkwd_mode, GLFitInfo *fitinfo,
                            GPFitVary *fitvary);
static void CY_init_fitvary_split(int split_pass, GLPkwdMode pkwd_mode,
                                  GLFitInfo *fitinfo, GPFitVary *fitvary);
static GLRtnCode CY_nllsqs(const GLSpectrum *spectrum, GPFitVary *fitvary,
                           double xtol, double ftol, GLFitRecord *fit);
static void CY_store_covar(GPFitVary *fitvary, double *fjac, int ldfjac,
                           double (*covarr) [GL_COVARR_DIM]);

/*
 * GP_cycle	global private routine to produce one fit.
 */

GPCycleType GP_cycle(const GLSpectrum *spectrum, GPFitVary *fitvary,
                     double xtol, double ftol, GLFitRecord *fit, GLCurve *curve)
{
switch(fit->used_parms.split_mode)
   {
   case GL_SPLIT_ALLOWED:
      return(CY_cycle_split(spectrum, fitvary, xtol, ftol, fit, curve));
      /* break; */
   case GL_SPLIT_NONE:
   default:
      return(CY_cycle_nosplit(spectrum, fitvary, xtol, ftol, fit, curve));
      /* break; */
   }
}

/*
 * CY_covary	private routine
 *
 * Given an M x N matrix, the problem is to determine the covariance
 * matrix corresponding to A, defined as
 *
 *    inverse(A**T * A)
 *
 * This subroutine completes the solution of the problem if it is
 * provided with the necessary information from the QR factorization,
 * with column pivoting, of A.  That is, if
 *
 *    A*P = Q*R
 *
 * where P is a permutation matrix, Q has orthogonal columns, and R is
 * an upper triangular matrix with diagonal elements of nonincreasing
 * magnitude, then covary expects the full upper triangle of R and the
 * permutation matrix P.  The covariance matrix is computed as
 *
 *    P * inverse(R**T * R) * P**T
 *
 * If A is nearly rank deficient, it may be desirable to compute the
 * covariance matrix corresponding to the linearly independent columns
 * of A.  To define the numerical rank of A, covary uses the tolerance
 * 'tol'.  If L is the largest integer such that
 *
 *    abs(R[L][L]) > tol * abs(R[0][0])
 *
 * then covary computes the covarience matrix corresponding to the first
 * L columns of R.  For K greater than L, column and row ipvt[K] of the
 * covariance matrix are set to zero.
 *
 * Parameters:
 *
 * n       - order of r
 * r[][]   - contains the full upper triangle of the matrix r
 * ldr     - leading dimension of r
 * ipvt[n] - pivot vector defining the permutation matrix P such that A*P = Q*R
 * tol     - defines the numerical rank of A
 * wa[n]   - work array
 * r[][]   - contains the square symmetric covariance matrix - the answer.
 */

static void CY_covary(int n, double *r, int ldr, int *ipvt, double tol,
                      double *wa)
{
int	i, j, k, l, jj, kk, jk, ik, ij, ji, index;
double	tolr, temp;

/*
 * r is really the 2 dimensional fortran array 'fjac'
 * where 'ldr' is size of first dimension.
 *
 * so "r[i][j]" is addressed as r[ldr*j + i]
 */

/*
 * form the inverse of r in the full upper triangle of r
 */

   tolr = tol * fabs(r[0]);

   for (k = kk = 0, l = -1; k < n; k++, kk += ldr + 1)
      {
      if (fabs(r[kk]) <= tolr)
         break;

      r[kk] = 1.0 / r[kk];

      if (k >= 1)
         for (j = 0, jk = ldr * k; j < k; j++, jk++)
            {
            temp = r[kk] * r[jk];
            r[jk] = 0.0;

            for (i = 0, ij = ldr * j, ik = ldr * k; i <= j; i++, ij++, ik++)
               r[ik] = r[ik] - temp * r[ij];
            }

      l = k;
      }

/*
 * form the full upper triangle of the inverse of r**t * r
 * in the full upper triangle of r.
 */

   if (l >= 0)
      for (k = kk = 0; k <= l; k++, kk += ldr + 1)
         {
         if (k >= 1)
            for (j = 0, jk = ldr * k; j < k; j++, jk++)
               {
               temp = r[jk];

               for (i = 0, ij = ldr * j, ik = ldr * k; i <= j; i++, ij++, ik++)
                  r[ij] = r[ij] + temp * r[ik];
               }

         temp = r[kk];

         for (i = 0, ik = ldr * k; i <= k; i++, ik++)
            r[ik] = temp * r[ik];
         }

/*
 * form the full lower triangle of the covariance matrix in the
 * strict lower triangle of r and in wa.
 */

   for (j = 0, jj = 0; j < n; j++, jj += ldr + 1)
      {
      for (i = 0, ij = ldr * j; i <= j; i++, ij++)
         {
         if (j > l)
            r[ij] = 0.0;

         if (ipvt[i] > ipvt[j])
            {
            index = ldr * (ipvt[j] - 1) + (ipvt[i] - 1);
            r[index] = r[ij];
            }
         else if (ipvt[i] < ipvt[j])
            {
            index = ldr * (ipvt[i] - 1) + (ipvt[j] - 1);
            r[index] = r[ij];
            }
         }

      wa[ipvt[j] - 1] = r[jj];
      }

/*
 * symmetrize the covariance matrix in r
 */

   for (j = jj = 0; j < n; j++, jj += ldr + 1)
      {
      for (i = 0, ij = ldr * j, ji = j; i <= j; i++, ij++, ji += ldr)
         r[ij] = r[ji];

      r[jj] = wa[j];
      }

return;
}

/*
 * dpmpar_	global private routine called from MINPACK fortran routine
 *		lmder() and associated subroutines.
 */

double dpmpar_(int *choice)
{
double	answer;

switch(*choice)
   {
   case 1:
      answer = CY_DPMPAR_ONE;
      break;
   case 2:
      answer = CY_DPMPAR_TWO;
      break;
   default:
      answer = 0;
      break;
   }

return(answer);
}

/*
 * CY_cycle_nosplit	private routine to fit with only one call
 *			to the non-linear least squares fit routine.
 */

static GPCycleType CY_cycle_nosplit(const GLSpectrum *spectrum,
                                    GPFitVary *fitvary, double xtol,
                                    double ftol, GLFitRecord *fit,
                                    GLCurve *curve)
{
int		i, vary_count, regn_width, denominator;
GLRtnCode	ret_code;

regn_width = fit->used_chanrange.last - fit->used_chanrange.first + 1;

if (fit->cycle_number == 1)
   CY_init_fitvary(fit->used_parms.pkwd_mode, &fit->fitinfo, fitvary);

   /* determine how many parameters are allowed to vary. */

      vary_count = GP_info_count(fitvary);

   denominator = regn_width - vary_count;

if (denominator <= 1)
   return(GP_CYCLE_DONE);

/* do non-linear lt-sq fit */

   ret_code = CY_nllsqs(spectrum, fitvary, xtol, ftol, fit);

   if (ret_code == GL_BADMALLOC)
      printf("malloc error in nllsqs\n");

/*
 * Using computed curve residuals, calculate the chi-sq.
 * The curve is saved for later use.
 */

   GL_curve(spectrum, fit, 1, curve);

   for (i = 0, fit->chi_sq = 0.0; i < regn_width; i++)
      fit->chi_sq += curve->resid[i] * curve->resid[i];

   fit->chi_sq = fit->chi_sq / denominator;

if (fit->lmder_info == CY_TERMINATE)
   return(GP_CYCLE_DONE);

return(GP_CYCLE_CONTINUE);
}

/*
 * CY_cycle_split	private routine to fit by calling the nonlinear
 *			least squares routine repeatedly, with either
 *			the linear or nonlinear parameters, and finally
 *			with all parameters.
 */

static GPCycleType CY_cycle_split(const GLSpectrum *spectrum,
                                  GPFitVary *fitvary, double xtol,
                                  double ftol, GLFitRecord *fit,
                                  GLCurve *curve)
{
int		split_pass, i, j, vary_count, regn_width, denominator;
GLRtnCode	ret_code;

regn_width = fit->used_chanrange.last - fit->used_chanrange.first + 1;

for (split_pass = 1; split_pass <= 5; split_pass++)
   {
   CY_init_fitvary_split(split_pass, fit->used_parms.pkwd_mode, &fit->fitinfo,
                         fitvary);

   /* determine how many parameters are allowed to vary. */

      vary_count = GP_info_count(fitvary);

   denominator = regn_width - vary_count;

   if (denominator <= 1)
      return(GP_CYCLE_DONE);

   /* do non-linear lt-sq fit */

      ret_code = CY_nllsqs(spectrum, fitvary, xtol, ftol, fit);

      if (ret_code == GL_BADMALLOC)
         printf("malloc error in nllsqs\n");

   /*
    * Using computed curve residuals, calculate the chi-sq.
    * The curve is saved for later use.
    */

      GL_curve(spectrum, fit, 1, curve);

      for (i = 0, fit->chi_sq = 0.0; i < regn_width; i++)
         fit->chi_sq += curve->resid[i] * curve->resid[i];

      fit->chi_sq = fit->chi_sq / denominator;

   if (fit->lmder_info == CY_TERMINATE)
      return(GP_CYCLE_DONE);

   /*
    * If this isn't the last pass, constrain the peak height.
    */

   /* 11/11/96 - aee - add check for deleted peaks and don't adjust them */

      if (split_pass < 5)
         for (j = 0; j < fit->fitinfo.npeaks; j++)
            if ((fabs(fit->fitinfo.peakinfo[j].height) < 2.0) &&
                (fit->fitinfo.peakinfo[j].height != 0.0))
               fit->fitinfo.peakinfo[j].height = 2.0;
   }

return(GP_CYCLE_CONTINUE);
}

/*
 * CY_fcn_	private routine whose address is passed to the MINPACK
 *		fortran routine lmder().
 */

static void CY_fcn_(int *m, int *n, double *x, double *fvec, double *fjac,
                    int *ldfjac, int *iflag)
{
GP_vector_to_info(CY_FCN_fitvary, x, CY_FCN_fitinfo);

switch(*iflag)
   {
   case 1:
      GP_calc_functions(&CY_FCN_chanrange, &CY_FCN_spectrum, CY_FCN_fitinfo,
                        fvec);
      break;
   case 2:
      if (GP_calc_jacobian(CY_FCN_chanrange.first, *m, &CY_FCN_spectrum,
                           CY_FCN_fitvary, CY_FCN_fitinfo, fjac, *ldfjac)
          != GL_SUCCESS)
         *iflag = CY_TERMINATE;
      break;
   case 0:
   default:
      return;
      /* break; */
   }
return;
}

/*
 * CY_init_fitvary	private routine to initialize the fit parameters
 *			for CY_cycle_nosplit().
 */

static void CY_init_fitvary(GLPkwdMode pkwd_mode, GLFitInfo *fitinfo,
                            GPFitVary *fitvary)
{
int	j;

/*
 * Initialize variable parameter structure
 * so that nothing is allowed to vary.
 */

   fitvary->intercept = GL_FALSE;
   fitvary->slope = GL_FALSE;
   fitvary->step_height = GL_FALSE;
   fitvary->avg_width = GL_FALSE;
   for (j = 0; j < fitvary->listlength; j++)
      {
      fitvary->peakvary[j].height = GL_FALSE;
      fitvary->peakvary[j].centroid = GL_FALSE;
      fitvary->peakvary[j].addwidth_511 = GL_FALSE;
      }

/* define variable parameters */

/*
 * nonpeak parameters (allow intercept and slope to vary)
 */

   fitvary->intercept = GL_TRUE;

   /*
    * Historically, the parameter "ish" controlled whether to allow
    * slope or step_height to vary.  Marie and Dick hardwired
    * slope to vary, and step_height to not vary, during Gauss IX
    * development.
    */

      fitvary->slope = GL_TRUE;

   /*
    * allow average peak width to vary if pkwd_mode != GL_PKWD_FIXED
    * and atleast one normal peak.
    */

      if (pkwd_mode != GL_PKWD_FIXED)
         for (j = 0; j < fitinfo->npeaks; j++)
            if (fitinfo->peakinfo[j].fixed_centroid == GL_FALSE)
               {
               fitvary->avg_width = GL_TRUE;
               break;
               }

/*
 * peak parameters:
 *   allow peak height to vary
 *   allow peak location to vary if peak not required
 *   allow extra width for 511 peaks to vary if there is more
 *     than one peak or the width mode is fixed.
 */

   for (j = 0, fitvary->npeaks = fitinfo->npeaks; j < fitvary->npeaks; j++)
      {
      fitvary->peakvary[j].height = GL_TRUE;

      if (fitinfo->peakinfo[j].fixed_centroid == GL_FALSE)
         fitvary->peakvary[j].centroid = GL_TRUE;

      if ((fitinfo->peakinfo[j].addwidth_511 != 0.0) &&
          ((fitinfo->npeaks > 1) || (pkwd_mode == GL_PKWD_FIXED)))
         fitvary->peakvary[j].addwidth_511 = GL_TRUE;
      }
}

/*
 * CY_init_fitvary_split	private routine to initialize the fit
 *				parameters for CY_cycle_split().
 */

static void CY_init_fitvary_split(int split_pass, GLPkwdMode pkwd_mode,
                                  GLFitInfo *fitinfo, GPFitVary *fitvary)
{
int	j;

/*
 * Initialize variable parameter structure
 * so that nothing is allowed to vary.
 */

   fitvary->intercept = GL_FALSE;
   fitvary->slope = GL_FALSE;
   fitvary->step_height = GL_FALSE;
   fitvary->avg_width = GL_FALSE;
   for (j = 0; j < fitvary->listlength; j++)
      {
      fitvary->peakvary[j].height = GL_FALSE;
      fitvary->peakvary[j].centroid = GL_FALSE;
      fitvary->peakvary[j].addwidth_511 = GL_FALSE;
      }

/* define variable parameters */

/*
 * nonpeak parameters
 *   on first, third, and fifth pass, allow intercept, slope, and possibly
 *   average peak width to vary.
 */

   if ((split_pass != 2) && (split_pass != 4))
      {
      fitvary->intercept = GL_TRUE;

      /*
       * Historically, the parameter "ish" controlled whether to allow
       * slope or step_height to vary.  Marie and Dick hardwired
       * slope to vary, and step_height to not vary, during Gauss IX
       * development.
       */

         fitvary->slope = GL_TRUE;

      /*
       * allow average peak width to vary if pkwd_mode != GL_PKWD_FIXED
       * and atleast one normal peak.
       */

         if (pkwd_mode != GL_PKWD_FIXED)
            for (j = 0; j < fitinfo->npeaks; j++)
               if (fitinfo->peakinfo[j].fixed_centroid == GL_FALSE)
                  {
                  fitvary->avg_width = GL_TRUE;
                  break;
                  }
      }

/* peak parameters */

   for (j = 0, fitvary->npeaks = fitinfo->npeaks; j < fitvary->npeaks; j++)
      {
      /* on first, third, and fifth pass, allow peak height to vary */

         if ((split_pass != 2) && (split_pass != 4))
            fitvary->peakvary[j].height = GL_TRUE;

      /* on second, fourth, and fifth pass
       *   allow peak location to vary if peak not required
       *   allow extra width for 511 peaks to vary if there is more
       *     than one peak or the width mode is fixed.
       */

         if ((split_pass != 1) && (split_pass != 3))
            {
            if (fitinfo->peakinfo[j].fixed_centroid == GL_FALSE)
               fitvary->peakvary[j].centroid = GL_TRUE;

            if ((fitinfo->peakinfo[j].addwidth_511 != 0.0) &&
                ((fitinfo->npeaks > 1) || (pkwd_mode == GL_PKWD_FIXED)))
               fitvary->peakvary[j].addwidth_511 = GL_TRUE;
            }
      }
}

/*
 * CY_nllsqs	private routine to prepare the input to lmder_(),
 *		invoke lmder(), and process the answer.
 */

static GLRtnCode CY_nllsqs(const GLSpectrum *spectrum, GPFitVary *fitvary,
                           double xtol, double ftol, GLFitRecord *fit)
{
int		vary_count;
int		i, ii, m, n, mode, nprint, maxfev;
double		factor, gtol;
double		*fvec, *x;
double		*wa, *diag, *qtf, *wa1, *wa2, *wa3, *wa4;
int		nfev, njev, *ipvt;
int		ldfjac; /* region width */

/* fjac is really a 2d array in fortran */
/* fjac[ldfjac][n] */

double	*fjac;

vary_count = GP_info_count(fitvary);
ldfjac = fit->used_chanrange.last - fit->used_chanrange.first + 1;
m = ldfjac;
n = vary_count;

/* allocate space */

   if ((fvec = (double *) calloc(2 * m, sizeof(double))) == NULL)
      {
      fit->lmder_info = CY_TERMINATE;
      return(GL_BADMALLOC);
      }
   wa4 = &fvec[m];

   if ((wa = (double *) calloc((6 * n), sizeof(double))) == NULL)
      {
      fit->lmder_info = CY_TERMINATE;
      free(fvec);
      return(GL_BADMALLOC);
      }
   x = wa;
   diag = &wa[n];
   qtf = &wa[2 * n];
   wa1 = &wa[3 * n];
   wa2 = &wa[4 * n];
   wa3 = &wa[5 * n];

   if ((ipvt = (int *) calloc(n, sizeof(int))) == NULL)
      {
      fit->lmder_info = CY_TERMINATE;
      free(fvec);
      free(wa);
      return(GL_BADMALLOC);
      }

   if ((fjac = (double *) calloc((ldfjac * n), sizeof(double))) == NULL)
      {
      fit->lmder_info = CY_TERMINATE;
      free(fvec);
      free(wa);
      free(ipvt);
      return(GL_BADMALLOC);
      }

/* initialize */

   /* copy fitinfo into the x vector using the fitvary flags as a guide. */

      GP_info_to_vector(fitvary, &fit->fitinfo, x);

   factor = 100.0;
   maxfev = 100 * (n + 1);
   gtol = 0.0;
   mode = 1;
   nprint = 2;
   nfev = 0;
   njev = 0;

   fit->uncertainty.intercept = 0;
   fit->uncertainty.slope = 0;
   fit->uncertainty.step_height = 0;
   fit->uncertainty.avg_width = 0;
   for (i = 0; i < fit->uncertainty.listlength; i++)
      {
      fit->uncertainty.peakinfo[i].height = 0;
      fit->uncertainty.peakinfo[i].centroid = 0;
      fit->uncertainty.peakinfo[i].addwidth_511 = 0;
      fit->uncertainty.peakinfo[i].fixed_centroid = GL_FALSE;
      }
   fit->uncertainty.npeaks = 0;

   /* set global variables for use in fcn_, as called from lmder */

      CY_FCN_chanrange.first = fit->used_chanrange.first;
      CY_FCN_chanrange.last = fit->used_chanrange.last;

      CY_FCN_spectrum.listlength = spectrum->listlength;
      CY_FCN_spectrum.nchannels = spectrum->nchannels;
      CY_FCN_spectrum.firstchannel = spectrum->firstchannel;
      CY_FCN_spectrum.count = spectrum->count;
      CY_FCN_spectrum.sigcount = spectrum->sigcount;

      CY_FCN_fitvary = fitvary;
      CY_FCN_fitinfo = &fit->fitinfo;

/* lmder_ uses and modifies the x vector. */

   lmder_(CY_fcn_, &m, &n, x, fvec, fjac, &ldfjac, &ftol, &xtol, &gtol,
          &maxfev, diag, &mode, &factor, &nprint, &fit->lmder_info, &nfev,
          &njev, ipvt, qtf, wa1, wa2, wa3, wa4);

if ((fit->lmder_info == 9) || (fit->lmder_info == 10) ||
    (fit->lmder_info == 0))
   {
   fit->lmder_info = CY_TERMINATE;
   free(fvec);
   free(wa);
   free(ipvt);
   free(fjac);
   return(GL_SUCCESS);
   }

if (fit->lmder_info == 8)
   fit->lmder_info = 4;

/* copy x into fitinfo using the fitvary flags as a guide. */

   GP_vector_to_info(fitvary, x, &fit->fitinfo);

/* wa is being reused as a workarea.  values from lmder() not used. */

   CY_covary(vary_count, fjac, ldfjac, ipvt, 0.0, wa);

/* Put sqrt of fjac's diagonal elements into fit uncertainty. */
/* wa is being reused as a workarea. */
/* The fitvary flags are used as a guide. */

   for (i = ii = 0; i < vary_count; i++, ii+= ldfjac + 1)
      wa[i] = sqrt(fjac[ii]);

   GP_vector_to_info(fitvary, wa, &fit->uncertainty);

/* 1/3/2003 - Egger enhancement for Blackwood */
/* store covariance of background terms */
/* this value is found in either location 1,0 or 0,1 (same value) */
/* whether the slope varies or the step_height varies. */
/* row and column 0 are for the background intercept. */
/* row and column 1 are for the background slope or step_height. */

   fit->back_covarr = fjac[1];

CY_store_covar(fitvary, fjac, ldfjac, fit->covarr);

free(fvec);
free(wa);
free(ipvt);
free(fjac);
return(GL_SUCCESS);
}

/*
 * CY_store_covar	private routine to finish computing the covariance
 *			matrix of the fit.
 */

static void CY_store_covar(GPFitVary *fitvary, double *fjac, int ldfjac,
                           double (*covarr) [GL_COVARR_DIM])
{
int	i1, i2, i3, k, counter;

/*
 * fjac is used and set as a 2 dimensional fortran array in lmder()
 * fjac[ldfjac][vary_count]
 *
 * where vary_count is count of TRUE items in fitvary
 *
 * so "fjac[i][j]" is addressed as fjac[(ldfjac * j) + i]
 */

i1 = counter = -1;

if (fitvary->intercept == GL_TRUE)
   counter++;
if (fitvary->slope == GL_TRUE)
   counter++;
if (fitvary->step_height == GL_TRUE)
   counter++;

if (fitvary->avg_width == GL_TRUE)
   {
   counter++;
   i1 = counter;
   }

for (k = 0; k < fitvary->npeaks; k++)
   {
   i2 = i3 = -1;

   covarr[k][0] = covarr[k][1] = covarr[k][2] = 0;

   if (fitvary->peakvary[k].height == GL_TRUE)
      {
      counter++;
      i2 = counter;
      }

   if (fitvary->peakvary[k].centroid == GL_TRUE)
      counter++;

   if (fitvary->peakvary[k].addwidth_511 == GL_TRUE)
      {
      counter++;
      i3 = counter;
      }

   if ((i1 != -1) && (i2 != -1))
         covarr[k][0] = fjac[(ldfjac * i2) + i1];

   if ((i1 != -1) && (i3 != -1))
         covarr[k][1] = fjac[(ldfjac * i3) + i1];

   if ((i2 != -1) && (i3 != -1))
         covarr[k][2] = fjac[(ldfjac * i3) + i2];
   }

for (; k < fitvary->listlength; k++)
   covarr[k][0] = covarr[k][1] = covarr[k][2] = 0;
}

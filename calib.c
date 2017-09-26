/*
 * Copyright 2017 Battelle Energy Alliance
 */

/*
 *  calib.c - contains the Gauss Algorithms for energy and width
 *		Calibration.
 */

#include <stdlib.h>		/* for calloc() and NULL */
#include <math.h>		/* for sqrt */
#include "GaussLib.h"

/* prototypes for private utilities */

static GLRtnCode CA_find_max(int order, int k, double array[25][25],
                             double *amax, int *ik, int *jk);
static GLRtnCode CA_linf(int x_count, const double *x, const double *y,
                         double *sigy, int coeff_count, double *polycoeffs,
                         double *chi2);
static void CA_matinv(int order, double array[25][25], double *det);
static GLRtnCode CA_swap_columns(int order, int k, double array[25][25],
                                 double *amax, int *ik, int *jk);

/*
 * GL_ecalib	public routine to calibrate energy
 */

GLRtnCode GL_ecalib(int count, const double *channel, const double *energy,
                    const double *sige, GLEgyEqnMode mode, GLboolean weighted,
                    GLEnergyEqn *ex)
{
double	*error;
double	polycoeffs[3];
int	i, coeff_count, ret_code;

if ((error = (double *) calloc(count, sizeof(double))) == NULL)
   return(GL_BADMALLOC);

/* prepare energy errors based on weight */

   switch(weighted)
      {
      case GL_FALSE:
         for (i = 0; i < count; i++)
            error[i] = 1.0;
         break;
      case GL_TRUE:
      default:
         for (i = 0; i < count; i++)
            error[i] = sige[i];
         break;
      }

switch(mode)
   {
   case GL_EGY_LINEAR:
      coeff_count = 2;
      break;
   case GL_EGY_QUADRATIC:
   default:
      coeff_count = 3;
      break;
   }

ret_code = CA_linf(count, channel, energy, error, coeff_count, polycoeffs,
                   &ex->chi_sq);

free(error);

if (ret_code != GL_SUCCESS)
   return(GL_BADCALBLINF);

ex->a = polycoeffs[0];
ex->b = polycoeffs[1];

switch(mode)
   {
   case GL_EGY_LINEAR:
      ex->c = 0.0;
      break;
   case GL_EGY_QUADRATIC:
   default:
      ex->c = polycoeffs[2];
      break;
   }

ex->mode = mode;

return(GL_SUCCESS);
}

/*
 * GL_wcalib	public routine to calibrate peak width
 */

GLRtnCode GL_wcalib(int count, const double *channel, const double *wid,
                    const double *sigw, GLWidEqnMode mode, GLboolean weighted,
                    GLWidthEqn *wid_eqn)
{
double	*data, *error;
double	polycoeffs[2];
int	i, coeff_count, ret_code;

if ((data = (double *) calloc(count, sizeof(double))) == NULL)
   return(GL_BADMALLOC);

if ((error = (double *) calloc(count, sizeof(double))) == NULL)
   {
   free(data);
   return(GL_BADMALLOC);
   }

/* prepare data to be fit */

   switch(mode)
      {
      case GL_WID_LINEAR:
         for (i = 0; i < count; i++)
            data[i] = wid[i];
         break;
      case GL_WID_SQRT:
      default:
         for (i = 0; i < count; i++)
            data[i] = wid[i] * wid[i];
         break;
      }

/* prepare width errors based on weight */

   switch(weighted)
      {
      case GL_FALSE:
         for (i = 0; i < count; i++)
            error[i] = 1.0;
         break;
      case GL_TRUE:
      default:
         for (i = 0; i < count; i++)
            error[i] = sigw[i] * 2.0 * wid[i];
         break;
      }

coeff_count = 2;

ret_code = CA_linf(count, channel, data, error, coeff_count, polycoeffs,
                   &wid_eqn->chi_sq);

free(data);
free(error);

if (ret_code != GL_SUCCESS)
   return(GL_BADCALBLINF);

wid_eqn->alpha = polycoeffs[0];
wid_eqn->beta = polycoeffs[1];
wid_eqn->mode = mode;

return(GL_SUCCESS);
}

/*
 * CA_find_max	private routine to find and return the largest value
 *		in a 2 dimensional array beyond the first k rows and
 *		columns that is also greater than or equal to the amax
 *		value passed in.
 *		Also return the location of the found value in ik,jk.
 */

static GLRtnCode CA_find_max(int order, int k, double array[25][25],
                             double *amax, int *ik, int *jk)
{
int	i, j;

for (i = k; i < order; i++)
   for (j = k; j < order; j++)
      if (fabs(*amax) <= fabs(array[i][j]))
         {
         *amax = array[i][j];
         *ik = i;
         *jk = j;
         }

if (*amax == 0)
   return(GL_FAILURE);

return(GL_SUCCESS);
}

/*
 * CA_linf	private routine to do a linear least squares fit of
 *		a list of (x,y) coordinates.  (fitting a polynomial
 *		to a list of points)
 */

static GLRtnCode CA_linf(int x_count, const double *x, const double *y,
                         double *sigy, int coeff_count, double *polycoeffs,
                         double *chi2)
{
int	i, j, k, mp;
double	sp[5], a[5][5], ai[25][25];
double	temp, det, xnm, calc_y;

/* set up matrix and vector */

   for (i = 0; i < coeff_count; i++)
      for (j = 0; j < coeff_count; j++)
         for (k = 0, mp = i + j, sp[j] = a[i][j] = 0.0; k < x_count; k++)
            {
            sp[j] = sp[j] + (pow(x[k], j) * y[k]) / (sigy[k] * sigy[k]);
            a[i][j] = a[i][j] + pow(x[k], mp) / (sigy[k] * sigy[k]);
            }

/* solve matrix equations */

   for (i = 0; i < coeff_count; i++)
      for (j = 0; j < coeff_count; j++)
         {
         temp = a[i][i] * a[j][j];
         if (temp <= 0)
            return(GL_BADCALBLINF);

         xnm = sqrt(temp);
         ai[i][j] = a[i][j] / xnm;
         }

   CA_matinv(coeff_count, ai, &det);

   /* if singular matrix, return error */

      if (det == 0)
         return(GL_BADCALBLINF);

   for (i = 0; i < coeff_count; i++)
      for (j = 0; j < coeff_count; j++)
         {
         temp = a[i][i] * a[j][j];
         if (temp <= 0)
            return(GL_BADCALBLINF);

         xnm = sqrt(temp);
         ai[i][j] = ai[i][j] / xnm;
         }

/* find solution */

   for (i = 0; i < coeff_count; i++)
      for (j = 0, polycoeffs[i] = 0.0; j < coeff_count; j++)
         polycoeffs[i] = polycoeffs[i] + ai[i][j] * sp[j];

/* don't calculate chi2 if xcount == coeff_count */

   *chi2 = 0.0;

   if (x_count != coeff_count)
      {
      for (j = 0; j < x_count; j++)
         {
         for (i = 0, calc_y = 0.0; i < coeff_count; i++)
            calc_y = calc_y + polycoeffs[i] * pow(x[j], i);

         temp = (y[j] - calc_y) / sigy[j];
         *chi2 = *chi2 + temp * temp;
         }

      *chi2 = *chi2 / (x_count - coeff_count);
      }

return(GL_SUCCESS);
}

/*
 * CA_matinv	private routine to invert a symmetric matrix and
 *		calculate its determinant.
 *
 * order   - size of matrix (<= 25)
 * array   - input matrix which is replaced by its inverse
 * det     - determinant of input matrix
 *
 * routine taken from Bevington's "Data Reduction and Error Analysis
 * for the Physical Sciences" -- see Gauss-Jordan elimination
 */

static void CA_matinv(int order, double array[25][25], double *det)
{
int	i, j, k;
int	ik[25], jk[25];
double	amax, save;

/* change to malloc 'order' space for ik & jk */
/* and free when done */

for (k = 0, *det = 1.0; k < order; k++)
   {
   amax = 0.0;

   if (CA_swap_columns(order, k, array, &amax, &ik[k], &jk[k])
       != GL_SUCCESS)
      {
      *det = 0.0;
      return;
      }

   /* interchange rows to put amax in array [k][k] */

      while (jk[k] < k)
         if (CA_swap_columns(order, k, array, &amax, &ik[k], &jk[k])
             != GL_SUCCESS)
            {
            *det = 0.0;
            return;
            }

      if (jk[k] > k)
         for (i = 0, j = jk[k]; i < order; i++)
            {
            save = array[i][k];
            array[i][k] = array[i][j];
            array[i][j] = -save;
            }

   /* accumulate elements of inverse matrix */

      for (i = 0; i < order; i++)
         if (i != k)
            array[i][k] = - array[i][k] / amax;

      for (i = 0; i < order; i++)
         for (j = 0; j < order; j++)
            if ((i != k) && (j != k))
               array[i][j] = array[i][j] + array[i][k] * array[k][j];

      for (j = 0; j < order; j++)
         if (j != k)
            array[k][j] = array[k][j] / amax;

      array[k][k] = 1.0 / amax;
      *det = *det * amax;
   }

/* restore ordering of matrix */

   for (k = order - 1; k >= 0; k--)
      {
      j = ik[k];

      if (j > k)
         for (i = 0; i < order; i++)
            {
            save = array[i][k];
            array[i][k] = - array[i][j];
            array[i][j] = save;
            }

      i = jk[k];

      if (i > k)
         for (j = 0; j < order; j++)
            {
            save = array[k][j];
            array[k][j] = - array[i][j];
            array[i][j] = save;
            }
      }
}

/*
 * CA_swap_columns	private routine called from CA_matinv() to
 *			rearrange matrix entries.
 */

static GLRtnCode CA_swap_columns(int order, int k, double array[25][25],
                                 double *amax, int *ik, int *jk)
{
double	save;
int	i, j;

if (CA_find_max(order, k, array, amax, ik, jk) != GL_SUCCESS)
   return(GL_FAILURE);

/* interchange columns to put amax in array[k][k] */

   /*
      observation:  if order, k, and array are not changed, and
                    amax was the value returned from the last call,
                    I expect CA_find_max() to return the same values
                    for amax, ik, and jk again...
      question:     so why is it being called again?
    */

   while (*ik < k)
      if (CA_find_max(order, k, array, amax, ik, jk) != GL_SUCCESS)
         return(GL_FAILURE);

   if (*ik > k)
      for (j = 0, i = *ik; j < order; j++)
         {
         /* swap array[k][j] with array[ik][j] */

            save = array[k][j];
            array[k][j] = array[i][j];
            array[i][j] = -save;
         }

return(GL_SUCCESS);
}

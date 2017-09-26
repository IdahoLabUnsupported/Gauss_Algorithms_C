/*
 * Copyright 2017 Battelle Energy Alliance
 */

/*
 * pprocess.c -  contains the Gauss postprocessing algorithms.
 */

#include "GaussLib.h"
#include "GLprivate.h"		/* for GPmin() */

/* prototypes for private utilities */

static void PP_seval(int npeaks, const double *energy,
                     const GLEfficiency *efficy, double *effic);

/*
 * GL_postprocess	public routine to process the summary using
 *			efficiencies.
 */

GLRtnCode GL_postprocess(const GLEfficiency *efficy, const GLSummary *summary,
                         GLPostInfo *info)
{
int		j, use_pkcount;
GLRtnCode	ret_code;

use_pkcount = GPmin(info->listlength, info->npeaks);

if (info->listlength < info->npeaks)
   ret_code = GL_OVRLMT;
else
   ret_code = GL_SUCCESS;

PP_seval(use_pkcount, summary->energy, efficy, info->efficiency);

for (j = 0; j < use_pkcount; j++)
   {
   if (info->efficiency[j] != 0.0)
      {
      info->intensity[j] = summary->area[j] / info->efficiency[j];

      info->sigi[j] = summary->siga[j] / info->efficiency[j];
      }
   else
      info->intensity[j] = info->sigi[j] = (double) 0.0;
   }

return(ret_code);
}

/*
 * GL_spline	routine to use the efficiency->en and efficiency->eff
 *		arrays to calculate the efficiency->coeff array.
 *
 * The coefficients b[i][2] i = 0, 1, ..., n-1 are computed for
 * a cubic interpolating spline
 *
 * s(x) = y[i] +
 *        (b[i][0] * (x - x[i])) +
 *        (b[i][1] * (x - x[i]) **2) +
 *        (b[i][2] * (x - x[i]) **3)
 *
 * for x[i] <= x <= x[i+1]
 *
 * input:
 *
 * n = number of data points or knots (n >= 4)
 * x = abscissas of the knots in strictly increasing order
 * y = ordinates of the knots
 *
 * output:
 *
 * b = array of spline coefficients as defined above
 *
 * using p to denote differentiation
 *
 * y[i] = s[x[i]]
 * b[i][0] = sp[x[i]]
 * b[i][1] = spp[x[i]]/2
 * b[i][2] = sppp[x[i]]/6 (derivative from the right)
 *
 * the accompanying function subprogram seval can be used
 * to evaluate the spline.
 */

void GL_spline(GLEfficiency *efficiency)
{
int	i, j;
double	*x, *y, tmp;
double	(*b)[GL_NUM_COEFFS];

x = efficiency->en;
y = efficiency->eff;
b = efficiency->coeff;

j = efficiency->neff - 1;

/*
 * set up tridiagonal system
 *
 * b[i][0] = diagonal, b[i][2] = offdiagonal, b[i][1] = right hand side
 */

   b[0][2] = x[1] - x[0];
   b[1][1] = (y[1] - y[0]) / b[0][2];

   for (i = 1; i < j; i++)
      {
      b[i][2] = x[i+1] - x[i];
      b[i][0] = 2.0 * (b[i-1][2] + b[i][2]);
      b[i+1][1] = (y[i+1] - y[i]) / b[i][2];
      b[i][1] = b[i+1][1] - b[i][1];
      }

/*
 * end conditions, third derivative at x[0] and x[j]
 * obtained from divided differences
 */

   b[0][0] = - b[0][2];
   b[j][0] = - b[j-1][2];
   b[0][1] = (b[2][1] / (x[3] - x[1])) - (b[1][1] / (x[2] - x[0]));
   b[j][1] = (b[j-1][1] / (x[j] - x[j-2])) -
               (b[j-2][1] / (x[j-1] - x[j-3]));
   b[0][1] = (b[0][1] * b[0][2] * b[0][2]) / (x[3] - x[0]);
   b[j][1] = - (b[j][1] * b[j-1][2] * b[j-1][2]) / (x[j] - x[j-3]);

/*
 * forward elimination
 */

   for (i = 1; i < efficiency->neff; i++)
      {
      tmp = b[i-1][2] / b[i-1][0];
      b[i][0] = b[i][0] - tmp * b[i-1][2];
      b[i][1] = b[i][1] - tmp * b[i-1][1];
      }

/*
 * back substitution
 */

   b[j][1] = b[j][1] / b[j][0];

   for (i = j - 1; i >= 0; i--)
      b[i][1] = (b[i][1] - (b[i][2] * b[i+1][1])) / b[i][0];

/*
 * b[i][1] is now the sigma[i] of the text
 *
 * compute the polynomial coefficients
 */

   b[j][0] = ((y[j] - y[j-1]) / b[j-1][2]) +
               b[j-1][2] * (b[j-1][1] + 2.0 * b[j][1]);

   for (i = 0; i < j; i++)
      {
      b[i][0] = ((y[i+1] - y[i]) / b[i][2]) -
                b[i][2] * (b[i+1][1] + 2.0 * b[i][1]);
      b[i][2] = (b[i+1][1] - b[i][1]) / b[i][2];
      b[i][1] = 3.0 * b[i][1];
      }

   b[j][1] = 3.0 * b[j][1];
   b[j][2] = b[j-1][2];
}


/*
 * PP_seval	private routine that evaluates the cubic spline function
 *		and returns peak efficiencies.
 *
 * cubic spline function is:
 *   effic[ict] = efficy->eff[i] +
 *                efficy->coeff[i,0] * (energy[ict] - efficy->en[i]) +
 *                efficy->coeff[i,1] * (energy[ict] - efficy->en[i]) **2 +
 *                efficy->coeff[i,2] * (energy[ict] - efficy->en[i]) **3
 *
 * where efficy->en[i] < energy[ict] < efficy->en[i+1] using Horner's rule.
 *
 * if energy[ict] <  efficy->en[0], then i = 0 is used.
 * if energy[ict] >= efficy->en[neff-1], then i = neff-1 is used.
 *
 * energy  is the array of abscissa at which the spline is to be evaluated
 * effic   is the array of efficiencies interpolated by the spline
 */

static void PP_seval(int npeaks, const double *energy,
                     const GLEfficiency *efficy, double *effic)
{
int	i, ict;
double	delta;


for (ict = 0; ict < npeaks; ict++)
   {
   for (i = 0; i < efficy->neff - 1; i++)
      if ((energy[ict] >= efficy->en[i]) &&
          (energy[ict] < efficy->en[i+1]))
         break;

   /* evaluate spline */

      delta = energy[ict] - efficy->en[i];

      effic[ict] = efficy->eff[i] + delta *
                     (efficy->coeff[i][0] + delta *
                      (efficy->coeff[i][1] + delta * efficy->coeff[i][2]));
   }
}

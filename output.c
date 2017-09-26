/*
 * Copyright 2017 Battelle Energy Alliance
 */

/*
 * output.c - contains the Gauss code to summarize or curve the fit results.
 */

#include <math.h>		/* for fabs(), sqrt(), log(), erf() */
#include <limits.h>		/* for max integer value */
#include <stdlib.h>		/* for malloc() and NULL */

#include "GaussLib.h"
#include "GLprivate.h"		/* for GPmin() */

/* should user access be added to this constraint? */

#define	OU_CONSTRAINT	((double) 10.0)

/* prototypes of private utilities */

static void OU_area(int npeaks, const GLSpectrum *spectrum,
                    const GLChanRange *chanrange, const GLFitInfo *fitinfo,
                    int *peak_order, double constant, double *a,
                    double *background, double *ratio);
static double OU_calc_constant(void);
static double OU_calc_miuu(double regn_index, int peak_index, double constant,
                           const GLFitInfo *fitinfo);
static double OU_calc_miuu_constrained(double miuu);
static double OU_chan_to_back(double regn_index, const GLFitInfo *fitinfo);
static double OU_chan_to_fitpeak(double regn_index, int peak_index,
                                 double constant, const GLFitInfo *fitinfo,
                                 double *miuu, double *miuu_constrained,
                                 double *exp_miuu);
static double OU_fit_to_resid(double fit, int count, double sigcount);
static int OU_npoints(int regn_width, int nplots_per_chan);
static GLRtnCode OU_sortindex(double *value, int count, int *order);

/*
**    comments from Gauss VII:
**
**	evaluate peak shape function - as needed by nlin2
**         gaussian + gamanal tail + hypermet tail
**
**      assume - all y vary; x varies except for required fits;
**               peak height must be greater than -50 counts;
**               all widths equal except 511 keV;
**               width can only vary 20 from initial value;
**               widths all fixed if all x fixed;
**
**      will return large residuals if constraints are violated
**
**    6/3/2003 - studying old code, see that the version of GAUSS VII
**               that Marie Putnam gave me had some constraints commented
**               out, which is why my code does not have the one constraining
**               the peakheight >= -50 and why my code does not add more
**               to residual if centroid was too close to either end of
**               fit region.
*/

/*
 * GL_curve	public routine to compute the coordinates of a fit,
 *		its background, and its residual.
 */

GLRtnCode GL_curve(const GLSpectrum *spectrum, const GLFitRecord *fit,
                   int nplots_per_chan, GLCurve *curve)
{
int		i, j, k, regn_width, regn_index, fit_index, spec_index;
double		constant, unused_miuu, unused_miuu_constr, unused_exp_miuu;
double		x, step;
GLRtnCode	ret_code;

curve->chanrange.first = fit->used_chanrange.first;
curve->chanrange.last = fit->used_chanrange.last;

curve->nplots_per_chan = nplots_per_chan;
regn_width = curve->chanrange.last - curve->chanrange.first + 1;

curve->npoints = OU_npoints(regn_width, curve->nplots_per_chan);
if (curve->npoints > curve->listlength)
   {
   curve->npoints = curve->listlength;
   ret_code = GL_OVRLMT;
   }
else
   ret_code = GL_SUCCESS;

curve->npeaks = fit->fitinfo.npeaks;

constant = OU_calc_constant();
step = 1.0 / curve->nplots_per_chan;

for (i = regn_index = 0; i < curve->npoints; regn_index++)
   {
   for (j = 0, x = regn_index;
        (i < curve->npoints) && (j < curve->nplots_per_chan);
        j++, i++, x += step)
      {
      curve->x_offset[i] = x;

      curve->back[i] = OU_chan_to_back(x, &fit->fitinfo);

      for (k = 0, curve->fitcurve[i] = 0; k < curve->npeaks; k++)
         {
         curve->fitpeak[k][i] = OU_chan_to_fitpeak(x, k, constant,
                                                   &fit->fitinfo, &unused_miuu,
                                                   &unused_miuu_constr,
                                                   &unused_exp_miuu);
         curve->fitcurve[i] += curve->fitpeak[k][i];
         curve->fitpeak[k][i] += curve->back[i];
         }

      curve->fitcurve[i] += curve->back[i];
      }

   fit_index = regn_index * curve->nplots_per_chan;
   spec_index = curve->chanrange.first - spectrum->firstchannel + regn_index;

   curve->resid[regn_index] = OU_fit_to_resid(curve->fitcurve[fit_index],
                                              spectrum->count[spec_index],
                                              spectrum->sigcount[spec_index]);
   }

return(ret_code);
}

/*
 * GL_curve_alloc	public routine to allocate memory for a GLCurve
 *			structure.
 */

GLCurve *GL_curve_alloc(int regn_width, int nplots_per_chan, int npeaks)
{
GLCurve	*curve;
int	listlength, stg_count, i, j;
double	*storage;

listlength = OU_npoints(regn_width, nplots_per_chan);
stg_count = (listlength * (npeaks + 3)) + regn_width;

storage = NULL;

if (((curve = (GLCurve *) malloc(sizeof(GLCurve))) == NULL) ||
    ((storage = (double *) malloc(stg_count * sizeof(double))) == NULL) ||
    ((curve->fitpeak = (double **) calloc(npeaks, sizeof(double))) == NULL))
   {
   if (curve != NULL)
      free(curve);
   if (storage != NULL)
      free(storage);
   return(NULL);
   }

for (i = j = 0; i < npeaks; i++, j+= listlength)
   curve->fitpeak[i] = &storage[j];

curve->x_offset = &storage[j];
j += listlength;
curve->fitcurve = &storage[j];
j += listlength;
curve->back = &storage[j];
j += listlength;
curve->resid = &storage[j];

curve->listlength = listlength;

return(curve);
}

/*
 * GL_curve_free	public routine to free the memory previously
 *			allocated for a GLCurve structure.
 */

void GL_curve_free(GLCurve *curve)
{
if (curve->npeaks > 0)
   free(curve->fitpeak[0]);

free(curve->fitpeak);
free(curve);
}

/*
 * GL_fitrec_alloc	public routine to allocate memory for a GLFitRecord
 *			structure.
 */

GLFitRecord *GL_fitrec_alloc(int listlength)
{
GLFitRecord	*fitrec;

if ((fitrec = (GLFitRecord *) malloc(sizeof(GLFitRecord))) == NULL)
   return(NULL);

if ((fitrec->fitinfo.peakinfo = (GLPeakInfo *)
     calloc(listlength, sizeof (GLPeakInfo))) == NULL)
   {
   free(fitrec);
   return(NULL);
   }
fitrec->fitinfo.listlength = listlength;

if ((fitrec->uncertainty.peakinfo = (GLPeakInfo *)
     calloc(listlength, sizeof (GLPeakInfo))) == NULL)
   {
   free(fitrec->fitinfo.peakinfo);
   free(fitrec);
   return(NULL);
   }
fitrec->uncertainty.listlength = listlength;

if ((fitrec->covarr = (double (*)[GL_COVARR_DIM])
     calloc(listlength, GL_COVARR_DIM * sizeof(double))) == NULL)
   {
   free(fitrec->fitinfo.peakinfo);
   free(fitrec->uncertainty.peakinfo);
   free(fitrec);
   return(NULL);
   }

return(fitrec);
}

/*
 * GL_fitrec_free	public routine to free the memory previously
 *			allocated for a GLFitRecord structure.
 */

void GL_fitrec_free(GLFitRecord *fitrec)
{
free(fitrec->fitinfo.peakinfo);
free(fitrec->uncertainty.peakinfo);
free(fitrec->covarr);
free(fitrec);
}

/*
 * GL_summ_alloc	public routine to allocate memory for a GLSummary
 *			structure.
 */

GLSummary *GL_summ_alloc(int listlength)
{
GLSummary	*summ;
double		*storage;
int		i;

if ((summ = (GLSummary *) malloc(sizeof(GLSummary))) == NULL)
   return(NULL);

/* allocate storage in minimum of calls */

   if ((storage = (double *) calloc(10 * listlength, sizeof(double))) == NULL)
      {
      free(summ);
      return(NULL);
      }

/* set pointers for each double array in summary structure */

   i = 0;
   summ->channel = &storage[i];
   i += listlength;
   summ->sigc = &storage[i];
   i += listlength;
   summ->height = &storage[i];
   i += listlength;
   summ->sigh = &storage[i];
   i += listlength;
   summ->wid = &storage[i];
   i += listlength;
   summ->sigw = &storage[i];
   i += listlength;
   summ->area = &storage[i];
   i += listlength;
   summ->siga = &storage[i];
   i += listlength;
   summ->energy = &storage[i];
   i += listlength;
   summ->sige = &storage[i];

summ->listlength = listlength;

return(summ);
}

/*
 * GL_summ_free		public routine to free previously allocated
 *			memory for the GLSummary structure.
 */

void GL_summ_free(GLSummary *summary)
{
free(summary->channel);
free(summary);
}

/*
 * GL_summarize		public routine to summarize the fit information.
 */

GLRtnCode GL_summarize(const GLSpectrum *spectrum, const GLFitRecord *fit,
                       GLSummary *summary)
{
int		j, k;
int		*peak_order;
double		constant, dbl_back, t1, t2, t3, sta;
double		*dbl_chan;
GLRtnCode	ret_code;

/*
 * should I replace this someday with:
 *   constant = OU_calc_constant() * sqrt(3.14159);   ??
 */

   constant = sqrt((double) 3.14159 / ((double) 4.0 * log((double) 2.0)));

summary->npeaks = GPmin(summary->listlength, fit->fitinfo.npeaks);

if (summary->listlength < fit->fitinfo.npeaks)
   ret_code = GL_OVRLMT;
else
   ret_code = GL_SUCCESS;

/* set up work areas */

   peak_order = NULL;
   dbl_chan = NULL;

   if (((peak_order = (int *) malloc(summary->npeaks * sizeof(int)))
        == NULL) ||
       ((dbl_chan = (double *) malloc(summary->npeaks * sizeof(double)))
        == NULL))
      {
      if (peak_order != NULL)
         free(peak_order);
      return(GL_BADMALLOC);
      }

/* determine peak order */

   for (j = 0; j < summary->npeaks; j++)
      dbl_chan[j] = fit->fitinfo.peakinfo[j].centroid +
                      fit->used_chanrange.first;

   if (OU_sortindex(dbl_chan, summary->npeaks, peak_order) != GL_SUCCESS)
      {
      peak_order[0] = 0;

      for (k = 1; k < summary->npeaks; k++)
         peak_order[k] = k;
      }

/* compute peak areas */

   OU_area(summary->npeaks, spectrum, &fit->used_chanrange, &fit->fitinfo,
           peak_order, constant, summary->area, &dbl_back, &summary->ratio);

   for (k = 0; k < summary->npeaks; k++)
      {
      j = peak_order[k];

	  /* 10/27/2010
	   * Egger found example of negative peakwidth, but fit is good.
	   * realize that in constraint code, peakwidth always squared,
	   * so negative is not prevented. Notice also that area() takes
	   * the absolute value of the width plus additions. Do same here.
	   */
      summary->wid[k] = fabs(fit->fitinfo.avg_width +
                        fit->fitinfo.peakinfo[j].addwidth_511);

      t1 = fit->uncertainty.peakinfo[j].addwidth_511 *
           fit->uncertainty.peakinfo[j].addwidth_511 +
           fit->uncertainty.avg_width * fit->uncertainty.avg_width +
           2.0 * fit->covarr[j][1];

      t2 = fit->uncertainty.peakinfo[j].height *
           fit->uncertainty.peakinfo[j].height;

      t3 = fit->covarr[j][2] + fit->covarr[j][0];

      sta = constant * constant *
            (fit->fitinfo.peakinfo[j].height *
             fit->fitinfo.peakinfo[j].height * t1 +
             summary->wid[k] * summary->wid[k] * t2 +
             2.0 * summary->wid[k] * fit->fitinfo.peakinfo[j].height * t3);

      if (sta >= 0.0)
         summary->siga[k] = sqrt(sta);
      else
         summary->siga[k] = (double) 0.0;

      if (t1 >= 0.0)
         {
         summary->sigw[k] = sqrt(t1);
         }
      else
         summary->sigw[k] = (double) 0.0;
      }

/* store channel info and compute energies */

   for (k = 0; k < summary->npeaks; k++)
      {
      j = peak_order[k];

      summary->channel[k] = dbl_chan[j];
      summary->sigc[k] = fit->uncertainty.peakinfo[j].centroid;

      GL_chan_to_e(&fit->used_ex, summary->channel[k], &summary->energy[k]);
      summary->sige[k] = (fit->used_ex.b +
                          2 * summary->channel[k] * fit->used_ex.c) *
                         fit->uncertainty.peakinfo[j].centroid;
      }

/* finish copying into summary */

   for (k = 0; k < summary->npeaks; k++)
      {
      j = peak_order[k];

      summary->height[k] = fit->fitinfo.peakinfo[j].height;
      summary->sigh[k] = fit->uncertainty.peakinfo[j].height;
      }

free(peak_order);
free(dbl_chan);

return(ret_code);
}

/*
 * GP_calc_functions	global private routine for use only in fcn()
 */

void GP_calc_functions(const GLChanRange *chanrange, const GLSpectrum *spectrum,
                       const GLFitInfo *fitinfo, double *resid)
{
int	i, k, top, spec_index;
double	constant, fit, x, unused_miuu, unused_miuu_constr, unused_exp_miuu;

constant = OU_calc_constant();

for (i = 0, top = chanrange->last - chanrange->first + 1, x = i;
     i < top; i++, x = i)
   {
   for (k = 0, fit = 0; k < fitinfo->npeaks; k++)
      fit += OU_chan_to_fitpeak(x, k, constant, fitinfo, &unused_miuu,
                                &unused_miuu_constr, &unused_exp_miuu);

   fit += OU_chan_to_back(x, fitinfo);

   spec_index = chanrange->first - spectrum->firstchannel + i;

   resid[i] = OU_fit_to_resid(fit, spectrum->count[spec_index],
                              spectrum->sigcount[spec_index]);
   }
}

/*
 * GP_calc_jacobian	global private routine for use only in fcn()
 */

GLRtnCode GP_calc_jacobian(int chanrange_first, int regn_width,
                           const GLSpectrum *spectrum, GPFitVary *fitvary,
                           const GLFitInfo *fitinfo, double *fjac, int ldfjac)
{
GLFitInfo	workarea;
double		*row;
double		constant, miuu, miuu_constrained, exp_miuu, fitpeak,
		denominator;
int		vary_count, i, ii, j, spec_index;
double		x;

/*
 * fjac is assumed by lmder() to be a two-dimensional fortran array 
 * fjac[ldfjac][vary_count]
 *
 * where ldfjac is regn_width, and vary_count is number of trues in fitvary
 *
 * so "fjac[i][j]" is addressed as fjac[(ldfjac * j) + i]
 */

/* allocate workspace */

   vary_count = GP_info_count(fitvary);

   if ((row = (double *) calloc(vary_count, sizeof(double))) == NULL)
      return(GL_BADMALLOC);

   if ((workarea.peakinfo = (GLPeakInfo *)
        calloc(fitinfo->npeaks, sizeof(GLPeakInfo))) == NULL)
      return(GL_BADMALLOC);

   workarea.listlength = workarea.npeaks = fitinfo->npeaks;

constant = OU_calc_constant();

for (i = ii = 0, x = i; i < regn_width; i++, ii = i, x = i)
   {
   workarea.intercept = 1.0;
   workarea.slope = i;
   workarea.step_height = 0.0;
   workarea.avg_width = 0.0;

   for (j = 0; j < workarea.npeaks; j++)
      {
      fitpeak = OU_chan_to_fitpeak(x, j, constant, fitinfo, &miuu,
                                   &miuu_constrained, &exp_miuu);

      workarea.peakinfo[j].height = exp_miuu;

      if ((denominator =
           fitinfo->avg_width + fitinfo->peakinfo[j].addwidth_511) != 0.0)
         {
         workarea.peakinfo[j].addwidth_511 = (2.0 * fitpeak * miuu * miuu) /
                                             denominator;

         if (constant != 0)
            workarea.peakinfo[j].centroid = (2.0 * fitpeak * miuu) /
                                            (denominator * constant);
         else
            workarea.peakinfo[j].centroid = 0.0;
         }
      else
         {
         workarea.peakinfo[j].centroid = 0.0;
         workarea.peakinfo[j].addwidth_511 = 0.0;
         }

      workarea.avg_width += workarea.peakinfo[j].addwidth_511;
      }

   spec_index = chanrange_first - spectrum->firstchannel + i;

   GP_info_to_vector(fitvary, &workarea, row);

   /* now copy "row" into row of fjac[][] */

      for (j = 0; j < vary_count; j++, ii += ldfjac)
         fjac[ii] = - row[j] / spectrum->sigcount[spec_index];
   }

free(row);
free(workarea.peakinfo);
return(GL_SUCCESS);
}

/*
 * FCODE_add_gamanal	private routine to add a gamanal tail
 *			to the fit curve.
 *
 *			It is commented out because it is currently
 *			not used.
 */

/*
static void FCODE_add_gamanal()
{
 add gamanal tail
 current fixed peak shape parms gta1 & gta2 make p[0] == 0

 this has never been verified under Gauss IX

   if (p[1] != 0.0)
      e2 = xu / p[1];
   else
      e2 = 0.0;

   e3 = p[2] * xu * xu;

   if (fabs(e2) > 100.0)
      e2 = (100.0 * e2) / fabs(e2);

   if (fabs(e3) > 100.0)
      e3 = (100.0 * e3) / fabs(e3);

   *fit = *fit + fitinfo->peakinfo[j].height * p[0] * exp(e2) *
                 (1.0 - exp(-e3));
}
 */

/*
 * FCODE_add_hienergy	private routine to add a hi energy component
 *			to the fit curve.
 *
 *			It is commented out because it is currently
 *			not used.
 */

/*
void static FCODE_add_hienergy()
{
 add high energy tail
 current fixed peak shape parms hsa1 & hsa2 make p[5] == 0

 this has never been verified under Gauss IX

   if (p[6] != 0.0)
      e2 = xu / p[6];
   else
      e2 = 0.0;

   e3 = p[2] * xu * xu;

   if (fabs(e2) > 100.0)
      e2 = (100.0 * e2) / fabs(e2);

   if (fabs(e3) > 100.0)
      e3 = (100.0 * e3) / fabs(e3);

   *fit = *fit + (fitinfo->peakinfo[j].height * p[5] * exp(-e2) *
                  (1.0 - exp(-e3)));
}
 */

/*
 * FCODE_add_hypermet	private routine to add a hyperment tail to
 *			the fit curve.
 *
 *			It is commented out because it is currently
 *			not used.
 */

/*
static void FCODE_add_hypermet()
{
 add hypermet tail
 current fixed peak shape parms hta1 & hta2 make p[3] == 0

 comment this out currently because not using...

   if (p[4] != 0.0)
      {
      e2 = xu / p[4];
      e3 = xu + (.5 / p[4]);
      }
   else
      {
      e2 = 0.0;
      e3 = xu;
      }

   if (fabs(e2) > 100.0)
      e2 = (100.0 * e2) / fabs(e2);

   if (fabs(e3) > 6.0)
      e3 = (6.0 * e3) / fabs(e3);

   *fit = *fit + (.5 * fitinfo->peakinfo[j].height * p[3] * exp(e2) * erfc(e3));
}
 */

/*
 * OU_add_more		private routine to add more to the fit area.
 *
 *			It is commented out because it currently is
 *			not used.
 */

/*
static void OU_add_more()
{
 add more to the area - what is this being added?
 current fixed peak shape parms gta1 & gta2 make p[0] == 0

 this has never been verified under Gauss IX

   for (i = 0; i < regn_width; i++)
      {
       set up xu.  what is it? need better name?

         if ((denominator = fitinfo->avg_width +
                            fitinfo->peakinfo[j].addwidth_511) != 0)

             I think the 'chanrange.first - 1' is wrong.
             I added it because nolonger added in FitRegn/store

             8/14/2000 - study fcode.for from original GAUSS VII code
                         to determine what should be used.

            xu = (i + 1 -
                  (fitinfo->peakinfo[j].centroid + chanrange.first - 1))
                 / denominator;
         else
            xu = 0.0;

      if (xu > 0.0)
         break;

       set up e2.

         if (p[1] != 0)
            e2 = xu / p[1];
         else
            e2 = (double) 0.0;

         if (fabs(e2) > 100.0)
            e2 = (double) 100.0 * (e2 / fabs(e2));

       set up e3.

         e3 = p[2] * xu * xu;

         if (fabs(e3) > 100.0)
            e3 = (double) 100.0 * (e3 / fabs(e3));

       is aa the gamanal tail referred to in first comment?

         aa = fitinfo->peakinfo[j].height * p[0] * exp(e2) * (1.0 - exp(-e3));

         if (aa >= .5)
            a[k] = a[k] + aa;
      }
}
 */

/*
 * FCODE_add_step	private routine to add step to the background.
 *
 *			It is commented out because it currently is
 *			not used.
 */

/*
static void FCODE_add_step()
{
 add background step
 current fixed peak shape parms sh1 & sh2 make fitinfo->step_height == 0

 comment this out currently because not using...

   e2 = xu;

   if (fabs(e2) > 6.0)
      e2 = (6.0 * e2) / fabs(e2);

   backstep = .5 * fitinfo->peakinfo[j].height * fitinfo->step_height *
              erfc(e2);
   *back = *back + backstep;
   *fit = *fit + backstep;
   d[2] = d[2] + (backstep / fitinfo->step_height);
}
 */

/*
 * OU_area	private routine to calculate the area of each peak
 *
 * a          is array to return areas of each peak in the region
 * ratio      is the ratio of netsum to gaussian area (want near 1.0)
 */

static void OU_area(int npeaks, const GLSpectrum *spectrum,
                    const GLChanRange *chanrange, const GLFitInfo *fitinfo,
                    int *peak_order, double constant, double *a,
                    double *background, double *ratio)
{
int	j, k, regn_width;
int	regn_i, spec_i;
double	areasum, integralsum;

/*
 * since this is private utility - assume chanrange is legal.
 */

regn_width = chanrange->last - chanrange->first + 1;

/* determine areas with integration */

   for (k = 0; k < npeaks; k++)
      {
      j = peak_order[k];

      a[k] = fitinfo->peakinfo[j].height *
             fabs(fitinfo->peakinfo[j].addwidth_511 + fitinfo->avg_width) *
             constant;
      }

/* determine areas with summation */
/* this assumes fitinfo->step_height == 0, for linear background */

   /* midpoint of background curve */

      *background = fitinfo->intercept +
                    (fitinfo->slope * (regn_width - 1) / 2.0);

   for (regn_i = 1, areasum = 0.0,
        spec_i = chanrange->first - spectrum->firstchannel;
        regn_i <= regn_width; regn_i++, spec_i++)
      areasum = areasum + spectrum->count[spec_i] - *background;

/* comparison of integral and summation areas for region */
/* this is the peak to background ratio                  */

   for (k = 0, integralsum = 0.0; k < npeaks; k++)
      integralsum = integralsum + a[k];

   if (integralsum != 0.0)
      *ratio = areasum / integralsum;
   else
      *ratio = 0.0;
}

/*
 * OU_calc_constant	private routine to calculate a constant that is
 *			frequently used in computing gaussian curves.
 */

static double OU_calc_constant(void)
{
double	constant;

   constant = (double) 1.0 / sqrt((double) 4.0 * log((double) 2.0));

return(constant);
}

/*
 * OU_calc_miuu		private routine to calculate intermediate value
 *			needed by several other routines.
 */

static double OU_calc_miuu(double regn_index, int peak_index, double constant,
                           const GLFitInfo *fitinfo)
{
int	j;
double	denominator, miuu;

j = peak_index;

if ((denominator =
     (fitinfo->avg_width + fitinfo->peakinfo[j].addwidth_511) * constant)
    != 0.0)
   {
   /*
    * 8/14/2000.  Both region index and peakinfo[j].centroid changed
    *             to start at 0.
    */

      miuu = (regn_index - fitinfo->peakinfo[j].centroid) / denominator;
   }
else
   miuu = 0.0;

return(miuu);
}

/*
 * OU_calc_miuu_constrained	private routine to calculate intermediate
 *				value needed by several other routines.
 */

static double OU_calc_miuu_constrained(double miuu)
{
double	miuu_constrained;

miuu_constrained = miuu;

if (miuu_constrained > OU_CONSTRAINT)
   miuu_constrained = OU_CONSTRAINT;
else if (miuu_constrained < - OU_CONSTRAINT)
   miuu_constrained = - OU_CONSTRAINT;

return(miuu_constrained);
}

/*
 * OU_chan_to_back	private routine to calculate the background
 *			at any channel in the region.
 */

static double OU_chan_to_back(double regn_index, const GLFitInfo *fitinfo)
{
double	background;

/*
 * 8/14/2000.  Saw in original GAUSS VII code that region indexing
 *             starts at 0.   See line 76 of fit2.for.
 *             Now, when one plot per channel, regn_index here ranges
 *             from 0 to region width - 1
 */

background = fitinfo->intercept + (fitinfo->slope * regn_index);

return(background);
}

/*
 * OU_chan_to_fitpeak	private routine to calculate the fit of
 *			a given peak in the region.
 */

static double OU_chan_to_fitpeak(double regn_index, int peak_index,
                                 double constant, const GLFitInfo *fitinfo,
                                 double *miuu, double *miuu_constrained,
                                 double *exp_miuu)
{
int	j;
double	fitpeak;

/*
 * 8/14/2000.  Saw in original GAUSS VII code that region indexing
 *             starts at 0.   See line 76 of fit2.for.
 *             Now, when one plot per channel, regn_index here ranges
 *             from 0 to region width - 1
 */

j = peak_index;

*miuu = OU_calc_miuu(regn_index, j, constant, fitinfo);
*miuu_constrained = OU_calc_miuu_constrained(*miuu);
*exp_miuu = exp(-(*miuu_constrained * *miuu_constrained));
fitpeak = fitinfo->peakinfo[j].height * *exp_miuu;

return(fitpeak);
}

/*
 * OU_fit_to_resid	private routine to calculate the residual of
 *			the fit at a given channel.
 */

static double OU_fit_to_resid(double fit, int count, double sigcount)
{
double	residual;

if (sigcount != 0.0)
   residual = (count - fit) / sigcount;
else
   residual = 0.0;

return(residual);
}

/*
 * OU_subtr_step	private routine to subtract the step background
 *			from the peak area.
 *
 *			It is commented out because it currently is not used.
 *

static void OU_subtr_step()
{
 subtract step background from area when fitinfo->step_height is non-zero
 current fixed peak shape parms sh1 & sh2 make fitinfo->step_height == 0

 this has not been verified....

   constant = (double) 1.0 / sqrt((double) 4.0 * log((double) 2.0));

   for (regn_i = 1, spec_i = chanrange.first - spectrum.firstchannel;
        regn_i <= regn_width; regn_i++, spec_i++)
      {
      for (k = 0, background = 0.0; k < npeaks; k++)
         {
         j = peak_order[k];

         if ((denominator =
             (fitinfo->avg_width + fitinfo->peakinfo[j].addwidth_511) *
             constant) != 0)

             I think the 'chanrange.first - 1' is wrong.
             I added it because nolonger added in FitRegn/store
            
             8/14/2000 - study fcode.for from original GAUSS VII code
                         to determine what should be used.

            xu = (regn_i -
                  (fitinfo->peakinfo[j].centroid + chanrange.first - 1)) /
                 denominator;
         else
            xu = 0.0;

         if (fabs(xu) > 6.0)
            xu = 6.0 * (xu / fabs(xu));

         background = background +
           (.5 * fitinfo->peakinfo[j].height * fitinfo->step_height * erfc(xu));
         }

      background = background + fitinfo->intercept +
                   (fitinfo->slope * (regn_i - chanrange.first));
      areasum = areasum + spectrum.count[spec_i] - background;
      }
}
 */

/*
 * OU_npoints	private routine to determine how many coordinates
 *		will be needed to plot the fit for the given region
 *		and resolution.
 */

static int OU_npoints(int regn_width, int nplots_per_chan)
{
int	npoints;

npoints = ((regn_width - 1) * nplots_per_chan) + 1;

return(npoints);
}

/*
 * OU_sortindex		private routine that takes a list of double
 *			precision numbers, and determines which is smallest.
 *			It's index is stored in order[0]...  and so on.
 */

static GLRtnCode OU_sortindex(double *value, int count, int *order)
{
int	i, j, min_index;
double	min_value, *list;

if (count < 1)
   return(GL_FAILURE);

if (count == 1)
   {
   order[0] = 0;
   return(GL_SUCCESS);
   }

if ((list = (double *) malloc(sizeof(double) * count)) == NULL)
   return(GL_BADMALLOC);

for (i = 0; i < count; i++)
   {
   list[i] = value[i];
   order[i] = 0;
   }

/* identify smallest remaining value */

   for (i = 0; i < count; i++)
      {
      min_index = 0;
      min_value = list[min_index];

      for (j = 1; j < count; j++)
         if (list[j] <= min_value)
            {
            min_value = list[j];
            min_index = j;
            }

      order[i] = min_index;
      list[min_index] = INT_MAX;
      }

free(list);
return(GL_SUCCESS);
}

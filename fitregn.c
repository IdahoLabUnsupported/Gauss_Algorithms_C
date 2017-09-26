/*
 * Copyright 2017 Battelle Energy Alliance
 */

/*
 * fitregn.c - contains starting point of Gauss fit algorithm
 *             and methods for evaluating fit results.
 */

#include <math.h>		/* for fabs() */
#include <stdlib.h>		/* for malloc() and NULL */
#include <stdio.h>		/* for printf() */
#include "GaussLib.h"
#include "GLprivate.h"		/* for checkfit(), cycle() */

/* prototypes of private utilities */

static GLFitRecList *FR_alloc_fitreclist_item(int max_npeaks);
static void FR_copy_ex(const GLEnergyEqn *source, GLEnergyEqn *destination);
static void FR_copy_fitrecord(const GLFitRecord *source,
                              GLFitRecord *destination);
static void FR_copy_info(const GLFitInfo *source, GLFitInfo *destination);
static void FR_copy_parms(const GLFitParms *source, GLFitParms *destination);
static void FR_copy_range(const GLChanRange *source, GLChanRange *destination);
static void FR_copy_wx(const GLWidthEqn *source, GLWidthEqn *destination);
static GLPeak *FR_create_peak(double channel, GLboolean fixed_centroid,
							  const GLEnergyEqn *ex);
static void FR_freelist(GLFitRecList **top);
static void FR_init_fitinfo(const GLChanRange *chanrange,
                            const GLSpectrum *spectrum, const GLEnergyEqn *ex,
                            const GLWidthEqn *wx, const GLPeakList *peaks,
                            GLboolean *arg_peak511, double *init_peakwidth,
                            GLFitInfo *fitinfo);
static void FR_orderlist(GLFitRecList **top);
static void FR_prunelist(int nout, GLFitRecList **top);
static void FR_set_cnvg(GLPkwdMode pkwd_mode, GLCCType cc_type,
                        double *ftol, double *xtol);

/*
 * GL_fitreclist_free	public routine to release results
 */

void GL_fitreclist_free(GLFitRecList *fitreclist)
{
if (fitreclist == NULL)
   return;

if (fitreclist->next != NULL)
   GL_fitreclist_free(fitreclist->next);

GL_fitrec_free(fitreclist->record);
free(fitreclist);

return;
}

/*
 * GL_fitregn	public routine to fit a region
 */

GLRtnCode GL_fitregn(const GLChanRange *chanrange, const GLSpectrum *spectrum,
                     const GLPeakList *peaks, const GLFitParms *fitparms,
                     const GLEnergyEqn *ex, const GLWidthEqn *wx,
                     GLFitRecList **fitlist)
{
GLRtnCode	ret_code;
int		i, j;
int		ncycles, index, *valid_npeaks;
GLCCType	cc_type;
GLFitRecord	*curr_record;
GLFitRecList	*curr_list, *last_list;
GLCurve		*curr_curve;
int		regn_width;
GLboolean	peak511;
double		init_peakwidth;
GPFitVary	fitvary;
double		ftol, xtol;

/* check input */

   if (peaks->npeaks > fitparms->max_npeaks)
      return(GL_BADNPKS);

/* initialize */

   ret_code = GL_SUCCESS;

   if (*fitlist != NULL) 
      FR_freelist(fitlist);

   if ((fitvary.peakvary = (GPPeakVary *)
        calloc(fitparms->max_npeaks, sizeof(GPFitVary))) == NULL)
      return(GL_BADMALLOC);
   fitvary.listlength = fitparms->max_npeaks;

   if ((valid_npeaks = (int *) malloc(sizeof(int) * fitparms->ncycle)) == NULL)
      {
      free(fitvary.peakvary);
      return(GL_BADMALLOC);
      }

   regn_width = chanrange->last - chanrange->first + 1;
   if ((curr_curve = GL_curve_alloc(regn_width, 1, fitparms->max_npeaks))
       == NULL)
      {
      free(fitvary.peakvary);
      free(valid_npeaks);
      return(GL_BADMALLOC);
      }

   if ((*fitlist = FR_alloc_fitreclist_item(fitparms->max_npeaks)) == NULL)
      {
      free(fitvary.peakvary);
      free(valid_npeaks);
      GL_curve_free(curr_curve);
      return(GL_BADMALLOC);
      }

      cc_type = fitparms->cc_type;
      FR_set_cnvg(fitparms->pkwd_mode, cc_type, &ftol, &xtol);

   for (i = 0; i < fitparms->ncycle; i++)
      valid_npeaks[i] = 0;

/* set up fit record storage */

   ncycles = 1;
   last_list = NULL;
   curr_list = *fitlist;
   curr_record = curr_list->record;
   curr_record->cycle_number = ncycles;
   FR_copy_range(chanrange, &curr_record->used_chanrange);
   FR_copy_parms(fitparms, &curr_record->used_parms);
   FR_copy_ex(ex, &curr_record->used_ex);
   FR_copy_wx(wx, &curr_record->used_wx);
   curr_record->lmder_info = 0;
   curr_record->chi_sq = 0.0;
   FR_init_fitinfo(chanrange, spectrum, ex, wx, peaks, &peak511,
                   &init_peakwidth, &curr_record->fitinfo);
   curr_record->uncertainty.npeaks = 0;

if (GP_cycle(spectrum, &fitvary, xtol, ftol, curr_record, curr_curve)
    == GP_CYCLE_CONTINUE)
   {
   /* save count of undeleted peaks */

      valid_npeaks[0] = curr_record->fitinfo.npeaks;
   }
else
   {
   FR_freelist(fitlist);
   free(fitvary.peakvary);
   free(valid_npeaks);
   GL_curve_free(curr_curve);
   return(GL_FAILURE);
   }

for (ncycles = 2; ncycles <= fitparms->ncycle; ncycles++)
   {
   index = ncycles - 1;

   /* set up next record */

      if ((curr_list->next = FR_alloc_fitreclist_item(fitparms->max_npeaks))
          == NULL)
         {
         free(fitvary.peakvary);
         free(valid_npeaks);
         GL_curve_free(curr_curve);

         ret_code = GL_BADMALLOC;
         break;
         }

      /* carry over peaks, peakshapes, and used_* into subsequent cycles */

         FR_copy_fitrecord(curr_record, curr_list->next->record);

      last_list = curr_list;
      curr_list = curr_list->next;
      curr_record = curr_list->record;

      curr_record->cycle_number = ncycles;

   /* add/delete component for next cycle */

      if (GP_checkfit(spectrum, curr_record, &fitvary, curr_curve,
                      valid_npeaks, &peak511, init_peakwidth) == GP_CYCLE_DONE)
         {
         FR_freelist(&last_list->next);
         break;
         }

   if (GP_cycle(spectrum, &fitvary, xtol, ftol, curr_record, curr_curve)
       == GP_CYCLE_CONTINUE)
      {
      /*
       * Save count of undeleted peaks.
       * (they are deleted if peakheight == 0)
       */

         for (j = 0, valid_npeaks[index] = curr_record->fitinfo.npeaks;
              j < curr_record->fitinfo.npeaks; j++)
            if (curr_record->fitinfo.peakinfo[j].height == 0.0)
               valid_npeaks[index]--;
      }
   else
      {
      FR_freelist(&last_list->next);
      break;
      }
   }   /* endfor ncycles */

/* order the fits by chi_sq, smallest first, and discard extra cycles */

   FR_orderlist(fitlist);
   FR_prunelist(fitparms->nout, fitlist);

free(fitvary.peakvary);
free(valid_npeaks);
GL_curve_free(curr_curve);
return(ret_code);
}

/*
 * GL_get_big_resid_peak	public routine to return first peak
 *                          in the fit with big residual.
 *                          big == 10 times "count+1"
 */

GLPeak *GL_get_big_resid_peak(const GLSpectrum *spectrum,
							  const GLFitRecord *fit)
{
GLChanRange chan_range;
int regn_width;
int npeaks;
int i, j, k;
GLCurve *curve;
double *fitpeak;
double count, diff;
GLPeak *answer;

answer = NULL;
chan_range = fit->used_chanrange;
regn_width = chan_range.last - chan_range.first + 1;
npeaks = fit->fitinfo.npeaks;

if ((curve = GL_curve_alloc(regn_width, 1, npeaks)) == NULL)
   {
   return(answer);
   }

if (GL_curve(spectrum, fit, 1, curve) != GL_SUCCESS)
   {
   GL_curve_free(curve);
   return(answer);
   }

k = chan_range.first - spectrum->firstchannel;

for (i = 0; i < npeaks; i++)
   {
   fitpeak = curve->fitpeak[i];

   // loop over region, comparing fit to spectrum
   for (j = 0; j < regn_width; j++)
      {
      count = spectrum->count[k + j];
      diff = fabs(fitpeak[j] - count);
      if (diff > (10 * (count + 1)))
         {
         answer = FR_create_peak(fit->fitinfo.peakinfo[i].centroid,
			 fit->fitinfo.peakinfo[i].fixed_centroid, &(fit->used_ex));
         GL_curve_free(curve);
         return(answer);
         }
      }
   }

GL_curve_free(curve);
return(answer);
}

/*
 * GL_get_neg_peak			public routine to return first peak
 *                          in the fit with negative height
 */

GLPeak *GL_get_neg_peak(const GLFitRecord *fit)
{
int i;
int npeaks;
double height;
GLPeak *answer;

answer = NULL;
npeaks = fit->fitinfo.npeaks;

for (i = 0; i < npeaks; i++)
   {
   height = fit->fitinfo.peakinfo[i].height;

   if (0 > height)
      {
      answer = FR_create_peak(fit->fitinfo.peakinfo[i].centroid,
		  fit->fitinfo.peakinfo[i].fixed_centroid, &(fit->used_ex));

      break;
      }
   }

return(answer);
}

/*
 * GL_get_outside_peak		public routine to return first peak in the fit
 *                          with centroid located outside region
 */

GLPeak *GL_get_outside_peak(const GLFitRecord *fit)
{
GLChanRange used_chanrange;
int i;
int npeaks;
double centroid;
GLPeak *answer;

answer = NULL;
used_chanrange = fit->used_chanrange;
npeaks = fit->fitinfo.npeaks;

for (i = 0; i < npeaks; i++)
   {
   centroid = fit->fitinfo.peakinfo[i].centroid;

   if ((centroid < used_chanrange.first) ||
	   (centroid > used_chanrange.last))
      {
      answer = FR_create_peak(centroid,
		  fit->fitinfo.peakinfo[i].fixed_centroid, &(fit->used_ex));

      break;
      }
   }

return(answer);
}

/*
 * GL_get_pos_neg_peak		public routine to return first "positive-negative"
 *                          peak pair found in the fit
 */

GLPeakList *GL_get_posneg_peakpair(const GLFitRecord *fit)
{
int npeaks;
int i, j;
double neg_height, neg_centroid;
double other_height, other_centroid;
double half_width;
GLPeak *pos_peak, *neg_peak;
GLPeakList *answer;

answer = NULL;
half_width = fit->fitinfo.avg_width / 2.0;
npeaks = fit->fitinfo.npeaks;

for (i = 0; i < npeaks; i++)
   {
   neg_height = fit->fitinfo.peakinfo[i].height;
   neg_centroid = fit->fitinfo.peakinfo[i].centroid;

   if (neg_height < 0)
      {
      /*
	   * check to see whether there is a positive height peak located
       * within half of the peak width of the negative height peak
	   */
      for (j = 0; j < npeaks; j++)
         {
         other_height = fit->fitinfo.peakinfo[j].height;
		 if (other_height > 0)
            {
            other_centroid = fit->fitinfo.peakinfo[j].centroid;
            if (fabs(neg_centroid - other_centroid) < half_width)
               {
               pos_peak = FR_create_peak(other_centroid,
				   fit->fitinfo.peakinfo[j].fixed_centroid, &(fit->used_ex));
               neg_peak = FR_create_peak(neg_centroid,
				   fit->fitinfo.peakinfo[i].fixed_centroid, &(fit->used_ex));

               if ((answer = GL_peaks_alloc(2)) != NULL)
                  {
                  GL_add_peak(pos_peak, answer);
                  GL_add_peak(neg_peak, answer);
                  }

               break;
               }
            }
         }
      }
   }

return(answer);
}

/*
 * FR_alloc_fitreclist_item	private routine to allocate memory
 *				for the GLFitRecList structure.
 */

static GLFitRecList *FR_alloc_fitreclist_item(int max_npeaks)
{
GLFitRecList	*fitreclist;

if ((fitreclist = (GLFitRecList *) malloc(sizeof(GLFitRecList))) == NULL)
   return(NULL);

fitreclist->next = NULL;

if ((fitreclist->record = GL_fitrec_alloc(max_npeaks)) == NULL)
   {
   free(fitreclist);
   return(NULL);
   }

return(fitreclist);
}

/*
 * FR_copy_ex	private routine to copy the GLEnergyEqn structure.
 */

static void FR_copy_ex(const GLEnergyEqn *source, GLEnergyEqn *destination)
{
destination->a = source->a;
destination->b = source->b;
destination->c = source->c;
destination->chi_sq = source->chi_sq;
destination->mode = source->mode;
}

/*
 * FR_copy_fitrecord	private routine to copy the GLFitRecord structure.
 */

static void FR_copy_fitrecord(const GLFitRecord *source,
                              GLFitRecord *destination)
{
int	i, j;

destination->cycle_number = source->cycle_number;

FR_copy_range(&source->used_chanrange, &destination->used_chanrange);

FR_copy_parms(&source->used_parms, &destination->used_parms);

FR_copy_ex(&source->used_ex, &destination->used_ex);

FR_copy_wx(&source->used_wx, &destination->used_wx);

destination->lmder_info = source->lmder_info;
destination->chi_sq = source->chi_sq;

FR_copy_info(&source->fitinfo, &destination->fitinfo);
FR_copy_info(&source->uncertainty, &destination->uncertainty);

for (i = 0; i < source->fitinfo.npeaks; i++)
   for (j = 0; j < GL_COVARR_DIM; j++)
      destination->covarr[i][j] = source->covarr[i][j];
}

/*
 * FR_copy_info		private routine to copy the GLFitInfo structure.
 */

static void FR_copy_info(const GLFitInfo *source, GLFitInfo *destination)
{
int	j;

destination->intercept = source->intercept;
destination->slope = source->slope;
destination->step_height = source->step_height;
destination->avg_width = source->avg_width;
destination->npeaks = source->npeaks;

for (j = 0; j < destination->npeaks; j++)
   {
   destination->peakinfo[j].height = source->peakinfo[j].height;
   destination->peakinfo[j].centroid = source->peakinfo[j].centroid;
   destination->peakinfo[j].addwidth_511 = source->peakinfo[j].addwidth_511;
   destination->peakinfo[j].fixed_centroid = source->peakinfo[j].fixed_centroid;
   }
}

/*
 * FR_copy_parms	private routine to copy the GLFitParms structure.
 */

static void FR_copy_parms(const GLFitParms *source, GLFitParms *destination)
{
destination->ncycle = source->ncycle;
destination->nout = source->nout;
destination->max_npeaks = source->max_npeaks;
destination->pkwd_mode = source->pkwd_mode;
destination->split_mode = source->split_mode;
destination->cc_type = source->cc_type;
destination->max_resid = source->max_resid;
}

/*
 * FR_copy_range	private routine to copy the GLChanRange structure.
 */

static void FR_copy_range(const GLChanRange *source, GLChanRange *destination)
{
destination->first = source->first;
destination->last = source->last;
}

/*
 * FR_copy_wx	private routine to copy the GLWidthEqn structure.
 */

static void FR_copy_wx(const GLWidthEqn *source, GLWidthEqn *destination)
{
destination->alpha = source->alpha;
destination->beta = source->beta;
destination->chi_sq = source->chi_sq;
destination->mode = source->mode;
}

/*
 * FR_create_peak	private routine to allocate and initialize a peak
 */

static GLPeak *FR_create_peak(double channel, GLboolean fixed_centroid,
							  const GLEnergyEqn *ex)
{
GLPeak *answer;

if ((answer = (GLPeak *) malloc(sizeof(GLPeak))) != NULL)
   {
   answer->channel = channel;
   answer->channel_valid = GL_TRUE;
   answer->type = GL_PEAK_CHANNEL;
   answer->fixed_centroid = fixed_centroid;
   GL_chan_to_e(ex, channel, &(answer->energy));
   answer->energy_valid = GL_TRUE;
   }

return(answer);
}

/*
 * FR_freelist	private routine to free the memory allocation of a
 *		GLFitRecList structure.
 */

static void FR_freelist(GLFitRecList **top)
{
GL_fitreclist_free(*top);
*top = NULL;

return;
}

/*
 * FR_init_fitinfo	private routine to initialize the data in a
 *			GLFitInfo structure.
 */

static void FR_init_fitinfo(const GLChanRange *chanrange,
                            const GLSpectrum *spectrum, const GLEnergyEqn *ex,
                            const GLWidthEqn *wx, const GLPeakList *peaks,
                            GLboolean *arg_peak511, double *init_peakwidth,
                            GLFitInfo *fitinfo)
{
GLboolean	peak511;
int		i, j, top, rounded_chan, index, chan_of_highest, regn_width;
double		xmid, energy, highest, channel;

peak511 = *arg_peak511;

/*
 * set initial parameters in fitinfo with
 * spectrum, region, and user peaklist
 */

/* nonpeak items */

   peak511 = GL_FALSE;

   index = chanrange->last - spectrum->firstchannel;

   fitinfo->intercept = 0.5 * (spectrum->count[index - 1] +
                               spectrum->count[index]);
   fitinfo->slope = 0.0;
   fitinfo->step_height = 0.0;

/*
 * 8/18/2000 - Larry Blackwood noted that as originally calculated,
 *             xmid was rounded to the closest channel.  He said that
 *             floating point was the "right" thing to calculate here.
 *             So I changed it.
 */

   xmid = (chanrange->first + chanrange->last + 1) / 2;

   /* initialize initial peakwidth for future use */

      if (GL_chan_to_w(wx, xmid, &fitinfo->avg_width) != GL_SUCCESS)
         fitinfo->avg_width = 1;	/* default setting */
      *init_peakwidth = fitinfo->avg_width;

/* peak parameters */

   fitinfo->npeaks = 0;

   for (i = j = 0; (i < peaks->npeaks) && (j < peaks->npeaks); i++)
      {
      if (peaks->peak[i].channel_valid == GL_TRUE)
         {
         rounded_chan = (int) (peaks->peak[i].channel + 0.5);

         if ((rounded_chan >= chanrange->first) &&
             (rounded_chan <= chanrange->last))
            {
            fitinfo->peakinfo[j].centroid = peaks->peak[i].channel -
                                            chanrange->first;

            index = rounded_chan - spectrum->firstchannel;
            fitinfo->peakinfo[j].height = spectrum->count[index] -
                                          fitinfo->intercept;

            fitinfo->peakinfo[j].addwidth_511 = 0.0;
            fitinfo->peakinfo[j].fixed_centroid = peaks->peak[i].fixed_centroid;

            if (peaks->peak[i].energy_valid == GL_TRUE)
               energy = peaks->peak[i].energy;
            else
               GL_chan_to_e(ex, peaks->peak[i].channel, &energy);

            /*
             * if any user peak is close to 511KeV, flag it as such.
             */

               if (peak511 == GL_FALSE)
                  if (fabs(energy - 511.0) <= GP_511PK_THRESHOLD)
                     {
                     peak511 = GL_TRUE;
                     fitinfo->peakinfo[j].addwidth_511 = fitinfo->avg_width;
                     }

            j++;
            }
         }
      }

   fitinfo->npeaks = j;

   /*
    * If no user peaks to add, then add one at the highest count
    * in the region.  But if this is too close to region ends, then
    * add a peak at the region midpoint.
    */

      if (fitinfo->npeaks <= 0)
         {
         for (i = chanrange->first - spectrum->firstchannel,
              top = chanrange->last - spectrum->firstchannel,
              j = 0, chan_of_highest = 0, highest = 0.0;
              i <= top;
              i++, j++)
            if (spectrum->count[i] > highest)
               {
               chan_of_highest = j;
               highest = spectrum->count[i];
               }

         regn_width = chanrange->last - chanrange->first + 1;

         if ((chan_of_highest < 2) || (chan_of_highest > regn_width - 3))
            chan_of_highest = regn_width / 2;

         fitinfo->peakinfo[0].centroid = chan_of_highest;
         channel = chanrange->first - spectrum->firstchannel + chan_of_highest;
         index = (int) channel;
         fitinfo->peakinfo[0].height = spectrum->count[index] -
                                         fitinfo->intercept;
         fitinfo->peakinfo[0].addwidth_511 = 0.0;
         fitinfo->peakinfo[0].fixed_centroid = GL_FALSE;

         GL_chan_to_e(ex, channel, &energy);

         /*
          * if peak is close to 511KeV, flag it as such.
          */

            if (fabs(energy - 511.0) <= 4*GP_511PK_THRESHOLD)
               {
               peak511 = GL_TRUE;
               fitinfo->peakinfo[0].addwidth_511 = fitinfo->avg_width;
               }

         fitinfo->npeaks = 1;
      }

*arg_peak511 = peak511;
}

/*
 * FR_orderlist		private routine to order the fit records by
 *			the value of chi squared in increasing order.
 */

static void FR_orderlist(GLFitRecList **top)
{
GLFitRecList	*curr_list, *prev_next_list, *next_list;
GLFitRecord	*temp_record;

for (curr_list = *top; curr_list != NULL; curr_list = curr_list->next)
   {
   for (prev_next_list = curr_list, next_list = curr_list->next;
        next_list != NULL;
        prev_next_list = next_list, next_list = next_list->next)
      {
      if (curr_list->record->chi_sq > next_list->record->chi_sq)
         {
         /* want next_list->record to come before curr_list->record */

            /* pull next_list from fitlist */

               prev_next_list->next = next_list->next;

            /* insert next_list just after curr_list */

               next_list->next = curr_list->next;
               curr_list->next = next_list;

            /* swap records between curr_list and next_list */

               temp_record = next_list->record;
               next_list->record = curr_list->record;
               curr_list->record = temp_record;

         /* reset next_list */

            next_list = prev_next_list;
         }
      }
   }

return;
}

/*
 * FR_prunelist		private routine to shorten a list of fit records
 *			to the desired length and free any remaining
 *			records.
 */

static void FR_prunelist(int nout, GLFitRecList **top)
{
int		count;
GLFitRecList	*curr_list, *previous;

if (*top == NULL)
   return;

for (count = 0, previous = NULL, curr_list = *top; count < nout;
     count++, previous = curr_list, curr_list = curr_list->next)
   if (curr_list->next == NULL)
      return;

if (previous == NULL)
   FR_freelist(top);
else
   FR_freelist(&previous->next);
}

/*
 * FR_set_cnvg	private routine to set the convergence criteria for use
 *		in the nonlinear least squares fitting routine.
 */

static void FR_set_cnvg(GLPkwdMode pkwd_mode, GLCCType cc_type,
                        double *ftol, double *xtol)
{
static double	ftolm[4] = {1.0e-4, 1.0e-3, 1.0e-3, 1.0e-2};
static double	xtolm[4] = {1.0e-4, 1.0e-3, 3.0e-3, 1.0e-2};
int		index;

/*
 * set up convergence criteria for minpack routines
 *
 * values used are a function of cc_type (formerly iconv)
 *   and/or pkwd_mode (formerly 'iw')
 */

/*
 * In Gauss VII, if used "peak width varies; constrained" mode,
 * and cc_type was GL_CC_LARGER, the resulting index was zero.
 * Now the resulting index is never zero.
 */

   *ftol = *xtol = 0.0;

   index = (int) cc_type;

/* define xtol and ftol values */

   if (cc_type == GL_CC_LARGER)
      switch(pkwd_mode)
         {
         case GL_PKWD_VARIES:
            index = GL_CC_LARGER_INC;
            break;
         case GL_PKWD_FIXED:
            index = 3;
            break;
         default:
            break;
         }

   *ftol = ftolm[index];
   *xtol = xtolm[index];
}

/*
 * FR_init_pkshap	private routine to initialize the fixed peak
 *			shape parameters.
 *
 *			It is commented out because these parameters
 *			are currently not used.
 */

/*
static void FR_init_pkshap(GLChanRange chanrange, GLEnergyEqn ex,
                           GLPkShape *pkshapes, double *p)
{
int	regn_width;
double	xmid, energy;

 *
 * 8/18/2000 - Larry Blackwood noted that as originally calculated,
 *             xmid was rounded to the closest channel.  He said that
 *             floating point was the "right" thing to calculate here.
 *             So it should be changed to the following:
 *
 *             xmid = (chanrange.first + chanrange.last + 1) / 2;
 *

 * compute fixed peak shape parameters

   regn_width = chanrange.last - chanrange.first + 1;
   xmid = chanrange.first + (regn_width / 2);
   GL_chan_to_e(&ex, xmid, &energy);

   p[0] = pkshapes->gta1 + pkshapes->gta2 * energy;
   p[1] = pkshapes->gts1 + pkshapes->gts2 * energy;
   p[2] = pkshapes->gtc;
   p[3] = pkshapes->hta1 + pkshapes->hta2 * energy;
   p[4] = pkshapes->hts1 + pkshapes->hts2 * energy;
   p[5] = pkshapes->hsa1 + pkshapes->hsa2 * energy;
   p[6] = pkshapes->hss1 + pkshapes->hss2 * energy;
}
 */

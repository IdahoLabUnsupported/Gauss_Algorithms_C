/*
 * Copyright 2017 Battelle Energy Alliance
 */

/*
 * checkfit.c - is part of the Gauss fit algorithm.
 */

/*
 * determines whether to add or delete a peak and recycle
 */

#include <stdio.h>		/* for printf() */
#include <math.h>		/* for fabs() */
#include "GaussLib.h"
#include "GLprivate.h"		/* for GPCycleType */

/* prototypes of private utilities */

static GLboolean CK_add_peak(GLChanRange chanrange, double max_resid,
                             const GLCurve *curve, const GLSpectrum *spectrum,
                             GLEnergyEqn *ex, GLboolean *arg_peak511,
                             int max_npeaks, GLPkwdMode pkwd_mode,
                             int current_valid_npeaks,
                             GPFitVary *fitvary, GLFitInfo *fitinfo);
static GPCycleType CK_dlt_peak(int regn_width, int current_valid_npeaks,
                               int *saved_valid_npeaks, int ncycles,
                               GLPkwdMode pkwd_mode, double init_peakwidth,
                               GLboolean *arg_peak511,
                               GLFitInfo *fitinfo, GPFitVary *fitvary);

/*
 * GP_checkfit		global private routine to add or delete a peak
 *			in preparation for another fit.
 *
 * saved_valid_npeaks - array containing count of undeleted peaks for each cycle
 * peak511	- indicates whether one of existing peaks is 511keV
 * init_peakwidth - value of fitinfo.avg_width from first cycle before varied
 */

GPCycleType GP_checkfit(const GLSpectrum *spectrum, GLFitRecord *fit,
                        GPFitVary *fitvary, const GLCurve *curve,
                        int *saved_valid_npeaks, GLboolean *arg_peak511,
                        double init_peakwidth)
{
int		j, ncycles, regn_width, reset_count;
GPCycleType	cycle_flag;
GLboolean	peak511;

peak511 = *arg_peak511;
ncycles = fit->cycle_number - 1;	/* count of completed fit cycles */
regn_width = fit->used_chanrange.last - fit->used_chanrange.first + 1;
reset_count = fit->fitinfo.npeaks;

cycle_flag = CK_dlt_peak(regn_width, saved_valid_npeaks[ncycles-1],
                         saved_valid_npeaks, ncycles,
                         fit->used_parms.pkwd_mode, init_peakwidth,
                         &peak511, &fit->fitinfo, fitvary);

switch(cycle_flag)
   {
   case GP_CYCLE_CONTINUE:
      if (CK_add_peak(fit->used_chanrange, fit->used_parms.max_resid,
                      curve, spectrum, &fit->used_ex, &peak511,
                      fit->used_parms.max_npeaks, fit->used_parms.pkwd_mode,
                      saved_valid_npeaks[ncycles-1], fitvary,
                      &fit->fitinfo) == GL_FALSE)
         return(GP_CYCLE_DONE);
      break;
   case GP_CYCLE_DLT:
      break;
   case GP_CYCLE_DONE:
   default:
      return(GP_CYCLE_DONE);
      /* break; */
   }

/* constrain the other parameters */

   for (j = 0; j < reset_count; j++)
      {
      if ((fit->fitinfo.peakinfo[j].height < 10.0) &&
          (fit->fitinfo.peakinfo[j].height != 0.0))
         fit->fitinfo.peakinfo[j].height = 10.0;

      if (fit->fitinfo.peakinfo[j].addwidth_511 < 0.0)
         fit->fitinfo.peakinfo[j].addwidth_511 = init_peakwidth;

/*
 * 8/14/2000.  looking at parg2.for from GAUSS VII original code.
 *             it constrains centroid to be between ISTART+2
 *             and IEND-2.
 *
 *             and my peakinfo[j].centroid is index into region starting at 0.
 */

/*
      if (fit->fitinfo.peakinfo[j].centroid < 2.0)
         fit->fitinfo.peakinfo[j].centroid = 2.0;

      if (fit->fitinfo.peakinfo[j].centroid > (regn_width - 2))
         fit->fitinfo.peakinfo[j].centroid = regn_width - 2;
 */
      if (fit->fitinfo.peakinfo[j].centroid < 2.0)
         fit->fitinfo.peakinfo[j].centroid = 2.0;

      if (fit->fitinfo.peakinfo[j].centroid > (regn_width - 3))
         fit->fitinfo.peakinfo[j].centroid = regn_width - 3;
      }

return(GP_CYCLE_CONTINUE);
}

/*
 * CK_add_peak	private routine to add a peak for the next fit
 *		if conditions warrant it.
 */

static GLboolean CK_add_peak(GLChanRange chanrange, double max_resid,
                             const GLCurve *curve, const GLSpectrum *spectrum,
                             GLEnergyEqn *ex, GLboolean *arg_peak511,
                             int max_npeaks, GLPkwdMode pkwd_mode,
                             int current_valid_npeaks,
                             GPFitVary *fitvary, GLFitInfo *fitinfo)
{
int		i, j, index, regn_width, chan_of_maxresid;
double		maxresid_found, fit_at_maxresid, channel, energy;
GLboolean	peak511;

peak511 = *arg_peak511;

regn_width = chanrange.last - chanrange.first + 1;

/* find the maximum residual and its corresponding channel and fit */

   for (i = chan_of_maxresid = 0, fit_at_maxresid = maxresid_found = 0.0;
        i < regn_width; i++)
      if (curve->resid[i] >= maxresid_found)
         {
         maxresid_found = curve->resid[i];
         index = i * curve->nplots_per_chan;
         fit_at_maxresid = curve->fitcurve[index];
         chan_of_maxresid = i;	/* same range as peakinfo[j].centroid */
         }

/*
 * If the maximum residual found is less than the threshold,
 * don't add a peak.
 */

   if (maxresid_found < max_resid)
      return(GL_FALSE);

/*
 * If the maximum residual is too close to an existing peak,
 * don't add a peak.
 */

   for (j = 0; j < fitinfo->npeaks; j++)
      if (fabs(chan_of_maxresid - fitinfo->peakinfo[j].centroid) < 1.0)
         return(GL_FALSE);

/*
 * If the maximum residual's neighboring residuals are not greater
 * than zero, don't add a peak.
 */

   /* fit curve residuals are indexed from zero. */
   /* channels of maxresiduals and centroids are indexed from zero now... */

      index = chan_of_maxresid;

   if ((index < 1) || (index > regn_width-2))
      return(GL_FALSE);

   if ((curve->resid[index-1] <= 0.0) && (curve->resid[index+1] < 0.0))
      return(GL_FALSE);

/* If there is no room for another peak, don't add one. */

   if ((fitinfo->npeaks + 1 > max_npeaks) ||
       (fitinfo->npeaks + 1 > fitinfo->listlength))
      return(GL_FALSE);

/* add a component */

   channel = chanrange.first - spectrum->firstchannel + chan_of_maxresid;
   index = (int) channel;

   j = fitinfo->npeaks;
   fitinfo->peakinfo[j].height = spectrum->count[index] - fit_at_maxresid;
   fitinfo->peakinfo[j].centroid = chan_of_maxresid;
   fitinfo->peakinfo[j].addwidth_511 = 0.0;
   fitinfo->peakinfo[j].fixed_centroid = GL_FALSE;

/* update which parameters are allowed to vary */

   if (pkwd_mode != GL_PKWD_FIXED)
      fitvary->avg_width = GL_TRUE;

   if (peak511 == GL_TRUE)
      {
      /*
       * If there was only one peak before, and it was the 511KeV peak,
       * now we can let its additional width vary.
       */

      if (current_valid_npeaks == 1)
         fitvary->peakvary[0].addwidth_511 = GL_TRUE;
      }
   else
      {
      GL_chan_to_e(ex, channel, &energy);

      if (fabs(energy - 511.0) <= GP_511PK_THRESHOLD)
         {
         peak511 = GL_TRUE;
         fitinfo->peakinfo[j].addwidth_511 = fitinfo->avg_width;
         fitvary->peakvary[j].addwidth_511 = GL_TRUE;
         }
      }

   fitvary->peakvary[j].height = GL_TRUE;
   fitvary->peakvary[j].centroid = GL_TRUE;

fitinfo->npeaks++;
fitvary->npeaks++;

*arg_peak511 = peak511;
return(GL_TRUE);
}

/*
 * CK_dlt_peak	private routine to delete a peak for the next fit
 *		if conditions warrant it.
 */

static GPCycleType CK_dlt_peak(int regn_width, int current_valid_npeaks,
                               int *saved_valid_npeaks, int ncycles,
                               GLPkwdMode pkwd_mode, double init_peakwidth,
                               GLboolean *arg_peak511,
                               GLFitInfo *fitinfo, GPFitVary *fitvary)
{
int		i, j, k, top, dlt_index;
double		threshold;
GLboolean	delete_flag, peak511;

peak511 = *arg_peak511;

/*
 * Loop over undeleted peak pairs and check for peaks closer than
 * 1/5 of peakwidth.
 */

   threshold = 0.2 * init_peakwidth;

   for (delete_flag = GL_FALSE, j = 0, top = fitinfo->npeaks - 1;
        (delete_flag != GL_TRUE) && (j < top); j++)
      {
      if (fitinfo->peakinfo[j].height != 0.0)

         for (k = j + 1; k < fitinfo->npeaks; k++)

            /* if peaks are too close */

            if ((fitinfo->peakinfo[k].height != 0.0) &&
                (fabs(fitinfo->peakinfo[j].centroid -
                      fitinfo->peakinfo[k].centroid) < threshold))
               {
               /* check of npeaks to prevent repeat of fit */
               /* with same number of components       */

                  if (ncycles > 1)
                     for (i = 0; i < ncycles; i++)
                        if ((current_valid_npeaks - 1) == saved_valid_npeaks[i])
                           return(GP_CYCLE_DONE);

               dlt_index = j;

               /* if height of kth peak is less than that of jth */
               /* then delete the kth peak instead.              */
               /* See Gauss VII manual, section 2.2.5, page 22   */

                  if (fitinfo->peakinfo[k].height <
                      fitinfo->peakinfo[j].height)
                     dlt_index = k;

               delete_flag = GL_TRUE;
               break;
               }
      }

if (delete_flag == GL_TRUE)
   {
   j = dlt_index;

   /*
    * if peak being deleted is the 511KeV peak, update the flag
    */

      if (fitinfo->peakinfo[j].addwidth_511 != 0)
         peak511 = GL_FALSE;

   /*
    * Take peak out of list of parameters, and zero the height.
    * Peak still takes up space in fitinfo, but its height,
    * centroid and additional width are no longer allowed to vary.
    */

      fitvary->peakvary[j].height = GL_FALSE;
      fitvary->peakvary[j].centroid = GL_FALSE;

      fitinfo->peakinfo[j].height = 0.0;

      /*
       * If previously two peaks and width mode not fixed,
       * don't let the addwidth_511 vary for the remaining peak.
       */

         if (current_valid_npeaks == 2)
            {
            if (pkwd_mode != GL_PKWD_FIXED)
               for (i = 0; i < fitinfo->npeaks; i++)
                  fitvary->peakvary[i].addwidth_511 = GL_FALSE;
            }
         else
            fitvary->peakvary[j].addwidth_511 = GL_FALSE;

   /*
    * Leave fitvary->avg_width alone.
    * If this deletes last peak, would need to update fitvary->avg_width.
    * But probably not possible to delete last peak.
    */

   return(GP_CYCLE_DLT);
   }

*arg_peak511 = peak511;
return(GP_CYCLE_CONTINUE);
}

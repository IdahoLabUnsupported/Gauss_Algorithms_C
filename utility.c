/*
 * Copyright 2017 Battelle Energy Alliance
 */

/*
 *  utility.c - contains miscellaneous public utilities for the Gauss
 *		Algorithms library.
 */

#include <math.h>		/* for log, sqrt */
#include <stdlib.h>		/* for malloc() and NULL */
#include "GaussLib.h"
#include "GLprivate.h"		/* for GPmin(), GPmax() */


/*
 * GL_add_chanpeak	public routine to add a peak to the list
 *			in terms of it's channel location.
 */

GLRtnCode GL_add_chanpeak(double channel, GLPeakList *peaks)
{
int	i;

i = peaks->npeaks;

if ((i+1) > peaks->listlength)
   return(GL_OVRLMT);

peaks->peak[i].type = GL_PEAK_CHANNEL;
peaks->peak[i].channel_valid = GL_TRUE;
peaks->peak[i].channel = channel;
peaks->peak[i].energy_valid = GL_FALSE;
peaks->peak[i].fixed_centroid = GL_FALSE;
peaks->npeaks++;

return(GL_SUCCESS);
}

/*
 * GL_add_egypeak	public routine to add a peak to the list
 *			in terms of it's energy location.
 */

GLRtnCode GL_add_egypeak(double energy, double sige, GLPeakList *peaks)
{
int	i;

i = peaks->npeaks;

if ((i+1) > peaks->listlength)
   return(GL_OVRLMT);

peaks->peak[i].type = GL_PEAK_ENERGY;
peaks->peak[i].channel_valid = GL_FALSE;
peaks->peak[i].energy_valid = GL_TRUE;
peaks->peak[i].energy = energy;
peaks->peak[i].sige = sige;
peaks->peak[i].fixed_centroid = GL_TRUE;
peaks->npeaks++;

return(GL_SUCCESS);
}

/*
 * GL_add_peak		public routine to add a peak to the list.
 */

GLRtnCode GL_add_peak(const GLPeak *peak, GLPeakList *peaks)
{
int	i;

i = peaks->npeaks;

if ((i+1) > peaks->listlength)
   return(GL_OVRLMT);

peaks->peak[i].type = peak->type;
peaks->peak[i].channel_valid = peak->channel_valid;
peaks->peak[i].channel = peak->channel;
peaks->peak[i].energy_valid = peak->energy_valid;
peaks->peak[i].energy = peak->energy;
peaks->peak[i].sige = peak->sige;
peaks->peak[i].fixed_centroid = peak->fixed_centroid;
peaks->npeaks++;

return(GL_SUCCESS);
}

/*
 * GL_chan_to_e		public routine to calculate the energy of
 *			a given channel.
 */

void GL_chan_to_e(const GLEnergyEqn *ex, double channel, double *energy)
{
switch(ex->mode)
   {
   case GL_EGY_LINEAR:
      *energy = ex->a + (ex->b * channel);
      break;
   case GL_EGY_QUADRATIC:
   default:
      *energy = ex->a + (ex->b * channel) + (ex->c * channel * channel);
      break;
   }
}

/*
 * GL_chan_to_w		public routine to calculate the peakwidth
 *			of a peak at a given channel.
 */

GLRtnCode GL_chan_to_w(const GLWidthEqn *wx, double channel, double *width)
{
double	temp;

temp = wx->alpha + (wx->beta * channel);

switch(wx->mode)
   {
   case GL_WID_LINEAR:
      *width = temp;
      return(GL_SUCCESS);
      /* break; */
   case GL_WID_SQRT:
   default:
      if (temp < 0)
         return(GL_FAILURE);
      else
         {
         *width = sqrt(temp);
         return(GL_SUCCESS);
         }
      /* break; */
   }
}

/*
 * GL_e_to_chan		public routine to calculate the channel for
 *			a given energy.
 */

GLRtnCode GL_e_to_chan(const GLEnergyEqn *ex, double energy, double *channel)
{
double	bsqr_4ac; 	/* from quadratic formula for ax**2 + bx + c */

if ((ex->mode == GL_EGY_LINEAR) || (ex->c == 0.0))
   {
   if (ex->b == 0.0)
      return(GL_FAILURE);
   else
      *channel = (energy - ex->a) / ex->b;
   }
else
   {
   bsqr_4ac = (ex->b * ex->b) - (4.0 * ex->c * (ex->a - energy));

   if (bsqr_4ac < 0.0)
      return(GL_FAILURE);
   else
      *channel = GPmax(0.0, (-ex->b + sqrt(bsqr_4ac)) / (2.0 * ex->c));
   }

return(GL_SUCCESS);
}

/*
 * GL_find_regnpks	public routine that returns a sublist of all
 *			the peaks that are found in the indicated region.
 */

GLRtnCode GL_find_regnpks(const GLChanRange *region, const GLPeakList *peaks,
                          GLPeakList *pks_in_rgn)
{
int	i;

if (pks_in_rgn->listlength <= 0)
   return(GL_OVRLMT);

pks_in_rgn->npeaks = 0;

for (i = 0; i < peaks->npeaks; i++)
   {
   if ((peaks->peak[i].channel_valid == GL_TRUE) &&
       (peaks->peak[i].channel >= region->first) &&
       (peaks->peak[i].channel <= region->last))
      {
      if (GL_add_peak(&peaks->peak[i], pks_in_rgn) == GL_OVRLMT)
         {
         return(GL_OVRLMT);
         }
      }
   }

return(GL_SUCCESS);
}

/*
 * GL_get_version	public routine to return the version string
 *                  for Gauss Algorithms.
 */

char *GL_get_version()
{
	return GAUSS_ALGORITHMS_VERSION;
}

/*
 * GL_peak_results_alloc  public routine to allocate memory
 *                        for a GLPeakSearchResults structure.
 */

GLPeakSearchResults *GL_peak_results_alloc(int peak_listlength,
										   int spectrum_nchannels)
{
GLPeakSearchResults	*results;

if ((results = (GLPeakSearchResults *)
	malloc(sizeof(GLPeakSearchResults))) == NULL)
   return(NULL);

if ((results->peaklist = GL_peaks_alloc(peak_listlength)) == NULL)
   {
   free(results);
   return(NULL);
   }

if ((results->refinements = (GLPeakRefinement *)
	calloc(peak_listlength, sizeof(GLPeakRefinement))) == NULL)
   {
   GL_peaks_free(results->peaklist);
   free(results);
   return(NULL);
   }

if ((results->crossproducts = (int *)
	calloc(spectrum_nchannels, sizeof(int))) == NULL)
   {
   free(results->refinements);
   GL_peaks_free(results->peaklist);
   free(results);
   return(NULL);
   }
results->listlength = spectrum_nchannels;

return(results);
}

/*
 * GL_peak_results_free		public routine to free memory for a
 *                          GLPeakSearchResults structure.
 */

void GL_peak_results_free(GLPeakSearchResults *results)
{
free(results->crossproducts);
free(results->refinements);
GL_peaks_free(results->peaklist);
free(results);
}

/*
 * GL_peaks_alloc	public routine to allocate memory for a GLPeakList
 *			structure.
 */

GLPeakList *GL_peaks_alloc(int listlength)
{
GLPeakList	*peaks;

if ((peaks = (GLPeakList *) malloc(sizeof(GLPeakList))) == NULL)
   return(NULL);

if ((peaks->peak = (GLPeak *) calloc(listlength, sizeof (GLPeak))) == NULL)
   {
   free(peaks);
   return(NULL);
   }
peaks->listlength = listlength;

return(peaks);
}

/*
 * GL_peaks_free	public routine to free memory for a GLPeakList
 *			structure.
 */

void GL_peaks_free(GLPeakList *peaks)
{
free(peaks->peak);
free(peaks);
}

/*
 * GL_prune_rqdpks	public routine that deletes required peaks from
 *			the list if they are too close to a regular
 *			(search) peak.
 */

GLRtnCode GL_prune_rqdpks(const GLWidthEqn *wx, const GLPeakList *searchpks,
                          const GLPeakList *curr_rqd, GLPeakList *new_rqd)
{
double		peakwid, threshold;
GLboolean	save_flag;
int		i, j;

/* compare required channels to search channels */

   for (i = 0; i < curr_rqd->npeaks; i++)
      if (curr_rqd->peak[i].channel_valid == GL_TRUE)
         {
         for (j = 0, save_flag = GL_TRUE; j < searchpks->npeaks; j++)
            if (searchpks->peak[j].channel_valid == GL_TRUE)
               {
               if (GL_chan_to_w(wx, searchpks->peak[j].channel, &peakwid)
                   != GL_SUCCESS)
                  peakwid = 0;
               if (peakwid <= 0.0)
                  peakwid = 3;
               threshold = .2 * peakwid;

               if (fabs(curr_rqd->peak[i].channel -
                        searchpks->peak[j].channel) <
                   threshold)
                  {
                  save_flag = GL_FALSE;
                  break;
                  }
               }

         /* if required channel not too close to a search channel, save it */

            if (save_flag == GL_TRUE)
               {
               if (GL_add_peak(&curr_rqd->peak[i], new_rqd) == GL_OVRLMT)
                  return(GL_OVRLMT);
               }
         }

return(GL_SUCCESS);
}

/*
 * GL_regions_alloc	public routine to allocate memory for a GLRegions
 *			structure.
 */

GLRegions *GL_regions_alloc(int listlength)
{
GLRegions	*regions;

if ((regions = (GLRegions *) malloc(sizeof(GLRegions))) == NULL)
   return(NULL);

if ((regions->chanrange =
       (GLChanRange *) calloc(listlength, sizeof (GLChanRange))) == NULL)
   {
   free(regions);
   return(NULL);
   }
regions->listlength = listlength;

return(regions);
}

/*
 * GL_regions_free	public routine to free memory for a GLRegions
 *			structure.
 */

void GL_regions_free(GLRegions *regions)
{
free(regions->chanrange);
free(regions);
}

/*
 * GL_set_sigcount	public routine to calculate the spectrum
 *			uncertainty.
 */

void GL_set_sigcount(GLSpectrum *spec)
{
int	i, top;
double	temp;

top = GPmin(spec->listlength, spec->nchannels);

for (i = 0; i < top; i++)
   {
   temp = GPmax(0.0, spec->count[i]);
   spec->sigcount[i] = sqrt(temp);

   if (spec->sigcount[i] <= 0.0)
      spec->sigcount[i] = .3;
   }

/* do small count correction for counts <= 10 based on GW Phillips, */
/* NIM 153 (1978), p. 449.                                          */

   for (i = 0; i < 2; i++)
      if (spec->count[i] <= 10)
         {
         temp = (spec->count[i] + spec->count[i+1] + spec->count[i+2]) / 3.0;
         spec->sigcount[i] = sqrt(temp);

         if (spec->sigcount[i] <= 0.0)
            spec->sigcount[i] = .5773503;
         }

   for (; i < top - 2; i++)
      if (spec->count[i] <= 10)
         {
         temp = spec->count[i-2] + spec->count[i+2] +
                (2 * (spec->count[i-1] + spec->count[i+1])) +
                (3 * spec->count[i]);
         temp = GPmax(0, temp / 9.0);
         spec->sigcount[i] = sqrt(temp);

         if (spec->sigcount[i] <= 0.0)
            spec->sigcount[i] = .3333333;
         }

   for (; i < top; i++)
      if (spec->count[i] <= 10)
         {
         temp = (spec->count[i-2] + spec->count[i-1] + spec->count[i]) / 3.0;
         spec->sigcount[i] = sqrt(temp);

         if (spec->sigcount[i] <= 0.0)
            spec->sigcount[i] = .5773503;
         }
}

/*
 * GL_update_peaklist	public routine to update a peaklist.
 */

void GL_update_peaklist(const GLEnergyEqn *ex, GLPeakList *peaks)
{
int	i;

if (ex != NULL)
   {
   for (i = 0; i < peaks->npeaks; i++)
      switch(peaks->peak[i].type)
         {
         case GL_PEAK_CHANNEL:
            GL_chan_to_e(ex, peaks->peak[i].channel, &peaks->peak[i].energy);
            peaks->peak[i].sige = 0;
            peaks->peak[i].energy_valid = GL_TRUE;
            break;
         case GL_PEAK_ENERGY:
            GL_e_to_chan(ex, peaks->peak[i].energy, &peaks->peak[i].channel);
            peaks->peak[i].channel_valid = GL_TRUE;
            break;
         default:
            break;
         }
   }
else
   {
   for (i = 0; i < peaks->npeaks; i++)
      switch(peaks->peak[i].type)
         {
         case GL_PEAK_CHANNEL:
            peaks->peak[i].energy_valid = GL_FALSE;
            break;
         case GL_PEAK_ENERGY:
            peaks->peak[i].channel_valid = GL_FALSE;
            break;
         default:
            break;
         }
   }
}

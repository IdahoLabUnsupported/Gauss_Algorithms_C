/*
 * Copyright 2017 Battelle Energy Alliance
 */

/*
 *  peaksrch.c - contains the Gauss peak search algorithm
 */

#include <math.h>		/* for log, sqrt, fabs */
#include <stdlib.h>		/* for malloc() and NULL */
#include <stdio.h>      /* for printf() */
#include "GaussLib.h"
#include "GLprivate.h"		/* for GPmin(), GPmax() */

#define	PS_UPDATE_INTERVAL	10
#define PS_MIN_PEAKWIDTH	1
#define	PS_MAX_PEAKWIDTH	10
#define	PS_MIN_SQWAV		(PS_MIN_PEAKWIDTH * 3)	/* for peak search */
#define PS_MAX_FITWIDTH_ODD	1001 /* for refining location of peak */

/* prototypes for utilities private to this file */

static double PS_fitpeak(int peak, double peakwid, int *counts,
						 int nchannels, int firstchannel,
						 GLPeakRefinement *answer);
static GLboolean PS_markpeak(int *int_pklist, int *int_pkcount, int newPeak,
							 double peakWidth, int max_allowed);
static int PS_squarewave_width(const GLWidthEqn *wx, double channel);

/*
 * GL_peaksearch - public routine to search for peaks in a region
 *
 *
 * First, construct a square wave with zero area, and width calculated
 * using the width equation for each 1K range of the spectrum.  The wave
 * looks approximately like:
 *
 *                 -----
 *                 |   |
 *                 |   |
 *     ***************************** <-- y=0
 *             |   |   |   |
 *             -----   -----
 *
 * Then take the cross product of the square wave curve with the uncertainties
 * of the spectrum counts, multiplying the square wave to the uncertainties,
 * a waves' width at a time.  Save the result of each multiplication.
 * (The uncertainties of the counts should be previously calculated with
 *  the subroutine "GL_set_sigcount" in utilities.c. The uncertainties
 *  are approximately the square root of the counts.)
 *
 * This process effectively reinforces the peak shapes and washes out
 * the noise and background.
 *
 * When the resulting shape is examined, and a local maximum is greater
 * than the threshold passed into the search function, the location is
 * tagged as a peak.
 *
 * Taking the cross product of the square wave with the counts is discussed in:
 *   K. Debertin & R. G. Helmer, Gamma- and X-Ray Spectrometry with
 *   Semiconductor Detectors, (Amsterdam: Elsevier Science B.V., 1988),
 *   p. 172-175. 
 *
 * Using the count uncertainties in the cross product is discussed in:
 *   C. M. McCullagh & R. G. Helmer, GAUSS VII A Computer Program for the
 *   Analysis of Gamma-Ray Spectra from Ge Semiconductor Spectrometers,
 *   (EGG-PHYS-5890, October 1982), p. 5-7.
 * This reference states that the count uncertainties are used as an
 * attempt to make the sensitivity more independent of the magnitude
 * of the counts.
 *
 */

GLRtnCode GL_peaksearch(const GLChanRange *chanrange, const GLWidthEqn *wx,
                        int threshold, const GLSpectrum *spectrum,
                        GLboolean refine, GLPeakSearchResults *results)
{
/*
 * Arguments:
 *
 * chanrange	defines first and last spectral channels for searching
 * wx			peak width calibration
 * threshold	used to determine when a peak is found
 *		(use 20 for low, 10 for medium, and 5 for high sensitivity)
 * spectrum		contains spectral counts & uncertainties at each channel
 * refine		indicates whether to refine peak locations with lls fitting
 * results		holds the found peaks, array of cross products, and
 * 				  intermediate values from fitting/refining the location
 */

/* Variables used in this subroutine: */

int		sqwav_wid;	/* width of middle part of square wave */
int		sqwav_y;	/* y coord of square wave at current channel */
double	peakwid;	/* peakwidth at begin of each 1K interval */
double	centroid;	/* location of each peak found */
int		hi_chan;	/* last channel to use making cross products */
int		mult_interval;	/* current multiple of update interval */
int     *int_sigcount; /* round uncertainties for use in loop */
int		crossprod;	/* cross product of current channel */
int		*int_pklist;	/* array to hold the found peak centroids */
int		pass_count; /* count of consecutive cross products above threshold */
int		int_pkcount;	/* current count of found peaks */
int		chan;		/* current channel of cross product */
int		i;		/* used to index over channels of square wave */
int		top;		/* last channel in square wave */
GLboolean	exceed_flag;	/* set true when find too many peaks */

if ((chanrange->last > (spectrum->nchannels + spectrum->firstchannel - 1)) ||
    (chanrange->first > chanrange->last))
   return(GL_BADCHNRNG);

if (threshold < 0)
   return(GL_BADTHRESH);

if (results->peaklist->listlength <= 0)
   return(GL_FAILURE);

/* calculate hi_chan to leave room for cross product */

   if (GL_chan_to_w(wx, (double) chanrange->last, &peakwid) != GL_SUCCESS)
      return(GL_BADPEAKWD);
   hi_chan = chanrange->last - (int) (3.0 * peakwid);

/* malloc space for working lists */

   if ((int_pklist = (int *) calloc(results->peaklist->listlength,
                                    sizeof(int))) == NULL)
      return(GL_BADMALLOC);

/* initialize list of cross products */

   for (i = 0; i < spectrum->nchannels; i++)
      results->crossproducts[i] = 0;

/* initialize sqwav_wid */

   sqwav_wid = PS_squarewave_width(wx, (double) chanrange->first);

/*
 * Determine first multiple of update interval after the start of the channels
 * to be searched over.
 */

   /* increment i until the interval multiples exceed the starting point */

      for (i = 1; chanrange->first >= (PS_UPDATE_INTERVAL * i); i++) ;

   /* now initialize the multiple of interval */

      mult_interval = PS_UPDATE_INTERVAL * i;

int_pkcount = 0;

/*
 * Make integer array of spectrum uncertainties for efficient
 * use inside the cross product loop.
 */
   if ((int_sigcount = (int *) calloc(spectrum->nchannels,
	                                  sizeof(int))) == NULL)
      {
      free(int_pklist);
      return(GL_BADMALLOC);
      }

   for (i = 0; i < spectrum->nchannels; i++)
      int_sigcount[i] = (int) spectrum->sigcount[i];

/*
 * Make first guess at peak locations using cross product with
 * zero area squarewave.
 */

   for (chan = chanrange->first; chan <= hi_chan; chan++)
      {
      /* calc new sqwav_wid when chan is multiple of interval */

         if (chan == mult_interval)
            {
            sqwav_wid = PS_squarewave_width(wx, (double) chan);
            mult_interval += PS_UPDATE_INTERVAL;
            }

      /*
       * Calculate the cross product of square wave and count uncertainties
       * at "chan".
       *
       * NOTE - GAUSS VII only used integer part of sigcount
       */

         /* from "chan" to "chan + sqwav_wid - 1" */

            sqwav_y = -1;
            for (i = chan - spectrum->firstchannel, crossprod = 0,
                 top = i + sqwav_wid; i < top; i++)
               crossprod += sqwav_y * int_sigcount[i];

         /* from "chan + sqwav_wid" to "chan + (sqwav_wid * 2) - 1" */

            sqwav_y = 2;
            for (top += sqwav_wid; i < top; i++)
               crossprod += sqwav_y * int_sigcount[i];

         /* from "chan + (sqwav_wid * 2)" to "chan + (sqwav_wid * 3) - 1" */

            sqwav_y = -1;
            for (top += sqwav_wid; i < top; i++)
               crossprod += sqwav_y * int_sigcount[i];

         /* store cross product for this channel */

            i = chan - spectrum->firstchannel;
            results->crossproducts[i] = crossprod;

      }  /* end for (chan = chanrange->first... */

/* review cross products to find peaks */

   exceed_flag = GL_FALSE;
   pass_count = 0;
   for (chan = chanrange->first + 2, i = chan - spectrum->firstchannel;
	    chan < hi_chan; i++, chan++)
      {
      /* calculate peak width */
      peakwid = PS_MIN_PEAKWIDTH;
      GL_chan_to_w(wx, (double) chan, &peakwid);
      sqwav_wid = PS_squarewave_width(wx, (double) chan);

	  if (peakwid < PS_MAX_PEAKWIDTH) /* use old algorithm for marking */
	     {
		 /*
          * If previous cross product was greater than the threshold,
          * and this cross product is decreasing,
          * and previous two cross products were the same or increasing,
          * then record peak at 'previous chan + 3/2 of wave width'
          * because this is indexing from the left end of the square wave.
          */

            if ((results->crossproducts[i-1] > threshold) &&
                (results->crossproducts[i-1] > results->crossproducts[i]) &&
                (results->crossproducts[i-1] >= results->crossproducts[i-2]))
               {
               if (PS_markpeak(int_pklist, &int_pkcount,
                               chan - 1 + (int) (1.5 * sqwav_wid), peakwid,
                               results->peaklist->listlength) == GL_FALSE)
                  {
                  exceed_flag = GL_TRUE;
                  break;
                  }

               pass_count = 0;
               }
         }
      else  /* use new algorithm for wider peaks */
         {
         if (results->crossproducts[i] >= threshold)
            {
            pass_count++;
            }
		 else /* cross product below threshold */
            {
            if (pass_count > 0)
               {
               /* want to mark half-way back plus 3/2 of wave width */
               if (PS_markpeak(int_pklist, &int_pkcount,
                       (int) (chan - (.5 * pass_count) + (1.5 * sqwav_wid)),
							   peakwid, results->peaklist->listlength)
                   == GL_FALSE)
                  {
                  exceed_flag = GL_TRUE;
                  break;
                  }
               }

            pass_count = 0;
            }
         }
      }  /* end for loop */

/* fine-tune peak locations using linear least-square fit */

   for (i = 0; i < int_pkcount; i++)
      {
      if ((GL_chan_to_w(wx, (double) int_pklist[i], &peakwid) != GL_SUCCESS) ||
          (peakwid <= 0.0))
         {
         centroid = int_pklist[i];
         results->refinements[i].raw_channel = int_pklist[i];
         results->refinements[i].refine_region.first = int_pklist[i];
         results->refinements[i].refine_region.last = int_pklist[i];
         results->refinements[i].net_area = 0;
         results->refinements[i].background = 0;
         results->refinements[i].refined_channel = int_pklist[i];
         }
      else
         {
#ifndef GL_LINUX
         centroid = PS_fitpeak(int_pklist[i], peakwid, (int *) spectrum->count,
                    spectrum->nchannels, spectrum->firstchannel,
					&results->refinements[i]);
#else
         centroid = PS_fitpeak(int_pklist[i], peakwid, spectrum->count,
                    spectrum->nchannels, spectrum->firstchannel,
					&results->refinements[i]);
#endif
         }

      if (refine != GL_TRUE)
         {
         centroid = int_pklist[i];
         }

	  if (GL_add_chanpeak(centroid, results->peaklist) != GL_SUCCESS)
         break;

      }

results->peaklist->npeaks = i;
free(int_pklist);
free(int_sigcount);

switch(exceed_flag)
   {
   case GL_TRUE:
      return(GL_OVRLMT);
      /* break; */
   case GL_FALSE:
   default:
      return(GL_SUCCESS);
      /* break; */
   }
}

/*
 * PS_fitpeak - private routine to fine-tune a peak centroid
 *
 *
 * First, determine the peak width.
 *
 * Then take the log of the count in each channel of the peak. This
 * takes data that is expected to be gaussian shaped and produces
 * data that is expected to be shaped like a parabola.
 *
 * Then fit a parabola to the "logged" data.
 *
 * Finally, find the parabola's maximum which coincides with the centroid
 * of the gaussian shape.
 *
 * The reference:
 *   K. Debertin & R. G. Helmer, Gamma- and X-Ray Spectrometry with
 *   Semiconductor Detectors, (Amsterdam: Elsevier Science B.V., 1988),
 *   pg. 175. 
 * says that the centroid can be determined to within .1 channels using
 * this fine-tuning if the statistical quality of the data allows.
 *
 */

static double PS_fitpeak(int centroid, double peakwid, int *counts,
						 int nchannels, int firstchannel,
						 GLPeakRefinement *answer)
{
/*
 * Arguments:
 *
 * centroid		the centroid to be fitted and refined
 * peakwid		the width of the peak to be fitted and refined
 * counts		the spectral counts
 * nchannels		the number of channels in the spectrum
 * firstchannel		the first channel of the spectrum
 * answer		the intermediate values & new centroid
 */

/* Variables used in this subroutine: */

int		low;			/* first channel of peak */
int		hi;				/* last channel of peak */
int		pcw;			/* peak width in channels */
int		hpcw;			/* half of peak width in channels */
int		i, j;			/* loop indices */
int		top;			/* indicates last index value in loop */
int		pre_avg_back;	/* average background before peak */
int		post_avg_back;	/* average background after peak */
int		avg_back;		/* average background of peak */
double	net_counts;		/* counts minus background in current channel */
double	x[PS_MAX_FITWIDTH_ODD];	/* x coords of parabola */
double	y[PS_MAX_FITWIDTH_ODD];	/* y coords of parabola */
double	B1;				/* 0th moment of y[i] wrt x[i] */
double	B2;				/* 1st moment of y[i] wrt x[i] */
double	B3;				/* 2nd moment of y[i] wrt x[i] */
double	sumx0;			/* 0th moment of x[i] */
double	sumx1;			/* 1st moment of x[i] */
double	sumx2;			/* 2nd moment of x[i] */
double	sumx3;			/* 3rd moment of x[i] */
double	sumx4;			/* 4th moment of x[i] */
double	detAlpha_a2;	/* |Alpha| * parabola's linear coefficient */
double	detAlpha_a3;	/* |Alpha| * parabola's quadratic coefficient */
double	newcentroid;	/* contains the centroid location returned */
float	diff;			/* diff between centroid & new centroid */
double	net_area;		/* sum area of peak */
double	background;		/* sum background of peak */

/* If peak is too close to one of the spectrum's ends, it cannot be fitted. */

   if ((centroid < firstchannel + peakwid) ||
       (centroid > (nchannels + firstchannel - 1 - peakwid)))
      {
      answer->raw_channel = centroid;
      answer->refine_region.first = centroid;
      answer->refine_region.last = centroid;
      answer->net_area = 0;
      answer->background = 0;
      answer->refined_channel = centroid;

      return(centroid);
      }

/* establish peak boundaries */

   pcw = (int) peakwid;
   pcw = GPmin(PS_MAX_FITWIDTH_ODD, pcw);
   hpcw = pcw / 2;
   pcw = (hpcw * 2) + 1; /* force to be odd */

   low = centroid - hpcw;
   hi = centroid + hpcw + 1;

/* estimate average background */

   for (i = low - 5 - firstchannel, pre_avg_back = 0,
        top = low - 1 - firstchannel; i <= top; i++)
      pre_avg_back += counts[i];
   pre_avg_back = pre_avg_back/5;

   for (i = hi + 1 - firstchannel, post_avg_back = 0,
        top = hi + 5 - firstchannel; i <= top; i++)
      post_avg_back += counts[i];
   post_avg_back = post_avg_back/5;

   avg_back = GPmin(pre_avg_back, post_avg_back);

/*
 * Using the Gaussian shaped peak data, construct parabolic shaped peak data
 * by taking the natural log of the number of counts in each channel.
 */

   answer->raw_channel = centroid;
   answer->refine_region.first = centroid - hpcw;
   answer->refine_region.last = centroid + hpcw;

   net_area = 0;
   background = 0;
   for (i = 0, j = low - firstchannel; i < pcw; i++, j++)
      {
      x[i] = i;
      net_counts = GPmax(1.0, counts[j] - avg_back);
      y[i] = log(net_counts);
      net_area += counts[j] - avg_back;
      background += avg_back;
      }

   answer->net_area = net_area;
   answer->background = background;

/*
 * To fit the parabolic shaped peak data, follow Bevington's method for
 * solving a matrix equation for a least-squares fit to a polynomial.
 * (P.R. Bevington, Data Reduction and Error Analysis for the Physical
 *  Sciences (New York: McGraw-Hill, 2003), p. 116-123, 238-244.)
 *
 * Define the row matrix "Beta" to have as its components the 0th, 1st, and
 * 2nd moments of y[i] with respect to x[i].
 *
 *     Beta = (B1 B2 B3)
 *
 * Define the row matrix "a" to have as its components the coefficients of
 * the polynomial describing the fitted parabola y = a1 + a2*x + a3*x**2.
 *
 *     a = (a1 a2 a3)
 *
 * Define the symmetric matrix Alpha to have as its components the 0th
 * through 4th moments of x[i].
 *
 *             (sumx0 sumx1 sumx2)
 *     Alpha = (sumx1 sumx2 sumx3)
 *             (sumx2 sumx3 sumx4)
 *
 * Then the matrix equation is
 *
 *     Beta = a * Alpha    (see equation B.6 in Bevington)
 *
 * The formula for a[i] is shown in equation B.33, and an example of
 * computing 3 x 3 determinants is in equations B.20 and B.21.
 *
 *
 * The following reference shows the solution for a[i] in terms of
 * Cramer's rule: (same solution in different notation)
 *    D.C. Lay, Linear Algebra and Its Applications, 2nd Ed.
 *    (Reading, MA: Addison-Wesley, 1999), p. 195-199.
 *
 */

   /* calculate components for the matrices Beta and Alpha */

      B1 = B2 = B3 = 0;
      sumx0 = sumx1 = sumx2 = sumx3 = sumx4 = 0;
      for (i = 0; i < pcw; i++)
         {
         B1 += y[i];
         B2 += y[i] * x[i];
         B3 += y[i] * x[i] * x[i];
         sumx0 += 1;
         sumx1 += x[i];
         sumx2 += x[i] * x[i];
         sumx3 += x[i] * x[i] * x[i];
         sumx4 += x[i] * x[i] * x[i] * x[i];
         }

   /*
    * Calculate the product of the determinate of Alpha times a3 by
    * calculating the determinate on the right in the following equation.
    *
    *                |sumx0 sumx1 B1|
    * |Alpha| * a3 = |sumx1 sumx2 B2|  (see B.33 in Bevington)
    *                |sumx2 sumx3 B3|
    *
    *
    */

      detAlpha_a3 = (sumx0 * sumx2 * B3);
      detAlpha_a3 += - (sumx0 * sumx3 * B2);
      detAlpha_a3 += - (sumx1 * sumx1 * B3);
      detAlpha_a3 += (sumx1 * B2 * sumx2);
      detAlpha_a3 += (B1 * sumx1 * sumx3);
      detAlpha_a3 += - (B1 * sumx2 * sumx2);

   if (detAlpha_a3 == 0.0)
      {
      /* cannot improve the centroid location */

         newcentroid = centroid;
		 answer->refined_channel = centroid;
      }
   else
      {
      /*
       * Calculate the product of the determinate of Alpha times a2 by
       * calculating the determinate on the right in the following equation.
       *
       *                |sumx0 B1 sumx2|
       * |Alpha| * a2 = |sumx1 B2 sumx3|  (see B.33 in Bevington)
       *                |sumx2 B3 sumx4|
       */

         detAlpha_a2 = (sumx0 * B2 * sumx4);
         detAlpha_a2 += - (sumx0 * sumx3 * B3);
         detAlpha_a2 += - (B1 * sumx1 * sumx4);
         detAlpha_a2 += (B1 * sumx3 * sumx2);
         detAlpha_a2 += (sumx2 * sumx1 * B3);
         detAlpha_a2 += - (sumx2 * B2 * sumx2);

      /*
       * Find the location of the maximum of the fitted parabola which
       * is also the location of the Gaussian centroid. The maximum of
       * the parabola will be where the slope is zero.
       *
       * The formula for the slope is y' = a2 + (2 * a3 * x).
       *
       * Solving for x where y' is zero yields: x = -a2 / (2 * a3)
       *
       * The values "|Alpha| * a3" and "|Alpha| * a2" have already
       * been calculated. Taking a ratio of these values as follows:
       *
       *     x = (-|Alpha| * a2)/(2 * |Alpha| * a3)
       *
       * is a convenient way to calculate x because |Alpha| cancels out
       * from the numerator and denominator.
       */

         newcentroid = (- detAlpha_a2 / (2.0 * detAlpha_a3)) + low;
		 answer->refined_channel = newcentroid;

      /*
       * If the new centroid is within a half peakwidth of the original, use
       * the new centroid.
       *
       * Otherwise, use the original centroid.
       */

         if ((diff = (float) fabs(centroid - newcentroid)) > hpcw)
              newcentroid = centroid;
      }

return(newcentroid);
}

/*
 * PS_markpeak - add a peak to the integer list if it is not closer than
 *               a peakwidth to the previous peak
 */
static GLboolean PS_markpeak(int *int_pklist, int *int_pkcount, int newPeak,
							 double peakWidth, int max_allowed)
{
if (*int_pkcount >= max_allowed)
   return(GL_FALSE);

if ((*int_pkcount <= 0) ||
    ((int_pklist[*int_pkcount-1] + peakWidth) < newPeak))
   {
   int_pklist[*int_pkcount] = newPeak;
   *int_pkcount = *int_pkcount + 1;
   }

   return(GL_TRUE);
}

/*
 * PS_squarewave_width - calculate width of the square wave
 */
static int PS_squarewave_width(const GLWidthEqn *wx, double channel)
{
double peakWidthDouble;
int sqwav_wid;

   if ((GL_chan_to_w(wx, channel, &peakWidthDouble) != GL_SUCCESS) ||
       (peakWidthDouble <= 0.0))
      {
      sqwav_wid = PS_MIN_SQWAV;
      }
   else
      {
      sqwav_wid = (int) peakWidthDouble;
      sqwav_wid = ((sqwav_wid/2) * 2) + 1;	/* ensure it is odd integer */
      sqwav_wid = GPmax(PS_MIN_SQWAV, sqwav_wid);
      }

   return sqwav_wid;
}

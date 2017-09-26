/*
 * Copyright 2017 Battelle Energy Alliance
 */
 
/*
 *  regnsrch.c - contains the Gauss region search algorithm.
 */

#include <stdlib.h>		/* for malloc() and NULL */
#include "GaussLib.h"
#include "GLprivate.h"	/* for GPmin(), GPmax() */

#define	RS_MIN_REGN_WIDTH	4	/* check with Dick... */
#define	RS_MIN_PKWID		1	/* for region search */

/* prototypes for utilities private to this file */

static void RS_dlt_small_bkg(const GLChanRange *chanrange, int firstchannel,
                             int nchannels, int maxrgnwid,
							 GLboolean *regn_flag);
static void RS_dlt_small_rgn(const GLChanRange *chanrange, double threshold,
                             const GLSpectrum *spectrum, int *background,
                             int maxrgnwid, GLboolean *regn_flag);
static void RS_init_rgnbkg(const GLWidthEqn *wx, const GLChanRange *chanrange,
                           const GLSpectrum *spectrum, int *background,
                           GLboolean *regn_flag);
static double RS_peakwidth(const GLWidthEqn *wx, double dchannel);
static void RS_prune_rgns(const GLWidthEqn *wx, const GLSpectrum *spectrum,
                          const GLPeakList *peaks, int *background,
                          double threshold, int maxrgnwid, GLRegions *regions);
static void RS_rgn_pad(const GLWidthEqn *wx, const GLChanRange *chanrange,
                       int irw, int irch, int maxrgnwid, GLRegions *regions);
static void RS_rgn_rmv(GLRegions *regions, int regn_index);
static void RS_rgnsforpks(const GLWidthEqn *wx, const GLPeakList *peaks,
                          int firstchannel, int nchannels,
                          GLboolean *regn_flag);
static GLRtnCode RS_store_rgns(const GLChanRange *chanrange, int firstchannel,
                               int nchannels, GLboolean *regn_flag,
                               GLRegions *regions);

/*
 * GL_regnsearch	public routine to search for regions in a spectrum.
 *
 * A standard smoothing technique is used to generate the background
 * curve of the entire spectrum.  Then the spectrum is compared to
 * this background.  Wherever a significant part of the spectrum lies
 * above the background, a region is flagged.
 *
 */

GLRtnCode GL_regnsearch(const GLChanRange *chanrange, const GLWidthEqn *wx,
                        double threshold, int irw, int irch,
                        const GLSpectrum *spectrum, const GLPeakList *peaks,
                        GLRgnSrchMode mode, int maxrgnwid, GLRegions *regions)
{
int			i;
GLboolean	*regn_flag;
int			*background;
int			regn_width;
GLRtnCode	ret_code;

if ((chanrange->last > spectrum->nchannels + spectrum->firstchannel - 1) ||
    (chanrange->first > chanrange->last))
   return(GL_BADCHNRNG);

if (irch < 0)
   return(GL_BADIRCH);
if (irw < 0)
   return(GL_BADIRW);
if (threshold < 0)
   return(GL_BADTHRESH);

/* initialize */

   if ((regn_flag =
         (GLboolean *) malloc(sizeof(GLboolean) * spectrum->nchannels)) ==
       NULL)
      return(GL_BADMALLOC);

   if ((background =
          (int *) malloc(sizeof(int) * spectrum->nchannels)) == NULL)
      {
      free(regn_flag);
      return(GL_BADMALLOC);
      }

   for (i = 0; i < spectrum->nchannels; i++)
      {
      regn_flag[i] = GL_FALSE;
      background[i] = 0;
      }

if (mode != GL_RGNSRCH_FORPKS)
   RS_init_rgnbkg(wx, chanrange, spectrum, background, regn_flag);

RS_rgnsforpks(wx, peaks, spectrum->firstchannel, spectrum->nchannels,
              regn_flag);

if (mode != GL_RGNSRCH_FORPKS)
   RS_dlt_small_rgn(chanrange, threshold, spectrum, background, maxrgnwid,
                    regn_flag);

RS_dlt_small_bkg(chanrange, spectrum->firstchannel, spectrum->nchannels,
                 maxrgnwid, regn_flag);

ret_code = RS_store_rgns(chanrange, spectrum->firstchannel,
                         spectrum->nchannels, regn_flag, regions);

if (mode != GL_RGNSRCH_FORPKS)
   RS_prune_rgns(wx, spectrum, peaks, background, threshold, maxrgnwid,
                 regions);

RS_rgn_pad(wx, chanrange, irw, irch, maxrgnwid, regions);

free(regn_flag);
free(background);

/* check answer for regions that are too wide */

   if (ret_code == GL_SUCCESS)
      {
      for (i = 0; i < regions->nregions; i++)
         {
         regn_width = regions->chanrange[i].last -
                      regions->chanrange[i].first + 1;
         if (regn_width > maxrgnwid)
            {
            ret_code = GL_BADRGNWD;
            }
         }
      }

return(ret_code);
}

/*
 * RS_dlt_small_bkg	private routine that joins regions that are close.
 */

static void RS_dlt_small_bkg(const GLChanRange *chanrange, int firstchannel,
                             int nchannels, int maxrgnwid,
							 GLboolean *regn_flag)
{
int	i, j, bottom_chan, top_chan;
int first_region_count, second_region_count;

/* If regions separated by one channel, join them. */

   bottom_chan = chanrange->first + 6 - firstchannel;
   bottom_chan = GPmax(1, bottom_chan);
   top_chan = chanrange->last - 6 - firstchannel;
   top_chan = GPmin((nchannels - 2), top_chan);

   first_region_count = 0;
   for (i = bottom_chan; i <= top_chan; i++)
      {
      if ((regn_flag[i] == GL_FALSE) &&
          (regn_flag[i-1] == GL_TRUE) &&
          (regn_flag[i+1] == GL_TRUE))
         {
         /* Count second region and decide if we can join. */

            second_region_count = 0;
            for (j = i + 1; j < top_chan; j++)
               {
               if (regn_flag[j] == GL_TRUE)
                  second_region_count++;
               else
                  break;
               }
            if (first_region_count + second_region_count <= maxrgnwid)
               regn_flag[i] = GL_TRUE;
         }

      if (regn_flag[i] == GL_TRUE)
         first_region_count++;
      else
         first_region_count = 0;
      }
}

/*
 * RS_dlt_small_rgn	private routine that deletes or enlarges small
 *			regions.
 */

static void RS_dlt_small_rgn(const GLChanRange *chanrange, double threshold,
                             const GLSpectrum *spectrum, int *background,
                             int maxrgnwid, GLboolean *regn_flag)
{
int	i, j, bottom_chan, top_chan;
int first_region_count, second_region_count;

/* Discard single channel regions unless both adjacent counts are greater */
/* than background; and at center, y greater than background + threshold. */

   bottom_chan = chanrange->first + 6 - spectrum->firstchannel;
   bottom_chan = GPmax(1, bottom_chan);
   top_chan = chanrange->last - 6 - spectrum->firstchannel;
   top_chan = GPmin((spectrum->nchannels - 2), top_chan);

   for (i = bottom_chan; i <= top_chan; i++)
      if ((regn_flag[i] == GL_TRUE) &&
          (regn_flag[i-1] == GL_FALSE) && (regn_flag[i+1] == GL_FALSE))
         {
         if ((spectrum->count[i] >= (background[i] +
                             (int) (threshold * spectrum->sigcount[i]))) &&
             (spectrum->count[i-1] >= background[i-1]) &&
             (spectrum->count[i+1] >= background[i+1]))
            {
            /* Don't create new region that is bigger than max allowed. */

               first_region_count = 3;
               for (j = i - 2; j >= bottom_chan; j--)
                  {
                  if (regn_flag[j] == GL_TRUE)
                     first_region_count++;
                  else
                     break;
                  }
               second_region_count = 0;
               for (j = i + 2; j <= top_chan; j++)
                  {
                  if (regn_flag[j] == GL_TRUE)
                     second_region_count++;
                  else
                     break;
                  }
               if (first_region_count + second_region_count <= maxrgnwid)
                  {
                  regn_flag[i-1] = GL_TRUE;
                  regn_flag[i+1] = GL_TRUE;
                  }
               else
                  regn_flag[i] = GL_FALSE;
            }
         else
            regn_flag[i] = GL_FALSE;
         }

   /* Check very first and very last channels for small regions. */

      if ((regn_flag[0] == GL_TRUE) && (regn_flag[1] == GL_FALSE))
            regn_flag[0] = GL_FALSE;

      if ((regn_flag[top_chan+1] == GL_TRUE) &&
          (regn_flag[top_chan] == GL_FALSE))
            regn_flag[top_chan+1] = GL_FALSE;
}

/*
 * RS_init_rgnbkg	private routine that initializes the background
 *			used to find regions.
 */

static void RS_init_rgnbkg(const GLWidthEqn *wx, const GLChanRange *chanrange,
                           const GLSpectrum *spectrum, int *background,
                           GLboolean *regn_flag)
{
int			i, j, k;
int         bottom_chan, top_chan, peakwidth, sum_top_chan;
int			sum;
double		d_peakwidth;
GLboolean	change;

/* Initial region search - using background. */

   bottom_chan = chanrange->first + 5 - spectrum->firstchannel;
   bottom_chan = GPmax(0, bottom_chan);
   top_chan = chanrange->last - 5 - spectrum->firstchannel;
   top_chan = GPmin((spectrum->nchannels - 1), top_chan);

   for (i = 0; i < 30; i++)
      {
      /* Set values in background[]. */

         for (j = bottom_chan; j <= top_chan; j++)
            {
            d_peakwidth = RS_peakwidth(wx, j);
            peakwidth = (int) ((d_peakwidth + .1) * 1.5);

            k = j - peakwidth;
            k = GPmax(0, k);
            sum_top_chan = j + peakwidth;
            sum_top_chan = GPmin((spectrum->nchannels - 1), sum_top_chan);
            for (sum = 0; k <= sum_top_chan; k++)
               if (regn_flag[k] == GL_TRUE)
                  sum += background[k];
               else
                  sum += spectrum->count[k];

            background[j] = (sum + peakwidth) / ((2 * peakwidth) + 1);
            }

      /* Set values in regn_flag[]. */

         change = GL_FALSE;
         for (j = bottom_chan; j <= top_chan; j++)
            if ((regn_flag[j] == GL_FALSE) &&
                ((background[j] + (int) (2 * spectrum->sigcount[j])) <=
                 spectrum->count[j]) &&
                (spectrum->count[j] > 1))
               {
               regn_flag[j] = GL_TRUE;
               change = GL_TRUE;
               }

      if (change == GL_FALSE)
         break;
      } /* end i loop */
}

/*
 * RS_peakwidth		private routine to calculate the peakwidth
 *			at the indicated channel.
 */
static double RS_peakwidth(const GLWidthEqn *wx, double dchannel)
{
double   dpeakwidth;

if (GL_chan_to_w(wx, dchannel, &dpeakwidth) != GL_SUCCESS)
   dpeakwidth = RS_MIN_PKWID;

dpeakwidth = GPmax(RS_MIN_PKWID, dpeakwidth);

return(dpeakwidth);
}

/*
 * RS_prune_rgns	private routine that deletes regions that are
 *			not needed for an existing peak and are not
 *			wide enough above background.
 */

static void RS_prune_rgns(const GLWidthEqn *wx, const GLSpectrum *spectrum,
                          const GLPeakList *peaks, int *background,
                          double threshold, int maxrgnwid, GLRegions *regions)
{
int			i, j;
double		d_peakwidth;
double		d_channel;
int			peakwidth;
GLboolean	regn_w_pk;
int			diff, maxdiff;
int			temp_peak, bottom_chan, top_chan, below_count, pkrgncount;

/* If region contains peak from peaklist (search peaks), keep it. */
/* Else if region has < peakwidth points above background near */
/* the highest point in the region, then discard the region.   */

   for (i = 0; i < regions->nregions; i++)
      {
      regn_w_pk = GL_FALSE;
      for (j = 0; j < peaks->npeaks; j++)
         {
         if (peaks->peak[j].channel_valid == GL_TRUE)
            {
            if ((peaks->peak[j].channel >= regions->chanrange[i].first) &&
                (peaks->peak[j].channel <= regions->chanrange[i].last))
               {
               regn_w_pk = GL_TRUE;
               break;
               }
            }
         }

      if (! regn_w_pk)
         {
         maxdiff = 0;
         bottom_chan = regions->chanrange[i].first - spectrum->firstchannel;
         top_chan = regions->chanrange[i].last - spectrum->firstchannel;
         for (j = bottom_chan, temp_peak = j; j <= top_chan; j++)
            {
            if ((diff = spectrum->count[j] - background[j]) > maxdiff)
               {
               maxdiff = diff;
               temp_peak = j;
               }
            }
      
         d_peakwidth = RS_peakwidth(wx, temp_peak);
         peakwidth = (int) (((d_peakwidth + .5) / 2.0) - 1.0);
         peakwidth = GPmax(peakwidth, RS_MIN_PKWID);

         below_count = 0;
         bottom_chan = temp_peak - peakwidth - spectrum->firstchannel;
         top_chan = temp_peak + peakwidth + 1 - spectrum->firstchannel;
         for (j = bottom_chan; j <= top_chan; j++)
            {
            if (spectrum->count[j] <= background[j])
               below_count++;
            }

         if ((below_count > 1) ||
             (spectrum->count[temp_peak] <= background[temp_peak] +
              (int) (threshold * spectrum->sigcount[temp_peak])))
            {
            RS_rgn_rmv(regions, i);
            i--;
            }
         } /* endif no peaks in region */
      } /* end i loop over regions */

/* If regions close, data between them is above background, */
/* and combined regions have <= 4 peaks, then join regions. */
/* Don't join regions if it would exceed max size allowed.  */

   for (i = 1; i < regions->nregions; i++)
      {

      /* Through Gauss Algorithms 1.30 (1/16/2012), a peakwidth
       * calculated at channel 3000 was used inside the loop.
       * From version 1.31 through 2.1 (1/17/2012 - 3/11/2013),
       * a peakwidth at a very small channel was incorrectly used.
       * Fixed in version 2.2 (4/2013) to use peakwidth at midpoint
       * between regions. */

      /* calculate peakwidth at midpoint between regions */

      d_channel = ((double) (regions->chanrange[i].first +
			  regions->chanrange[i-1].last)) / 2.0;
      d_peakwidth = RS_peakwidth(wx, d_channel);
      peakwidth = (int) (d_peakwidth + .5);

      if (regions->chanrange[i].first - regions->chanrange[i-1].last
          <= peakwidth)
         {
         bottom_chan = regions->chanrange[i-1].last + 1 -
			 spectrum->firstchannel;
         top_chan = regions->chanrange[i].first - spectrum->firstchannel;
         for (j = bottom_chan; j < top_chan; j++)
            {
            if (spectrum->count[j] < background[j])
               break;
            }

         if (j >= top_chan)
            {
            pkrgncount = 0;
            for (j = 0; j < peaks->npeaks; j++)
               {
               if (peaks->peak[j].channel_valid == GL_TRUE)
                  {
                  if ((peaks->peak[j].channel >= regions->chanrange[i-1].first)
                      &&
                      (peaks->peak[j].channel <= regions->chanrange[i].last))
                     pkrgncount++;
                  }
               }

            if ((pkrgncount <= 4) &&
                (regions->chanrange[i].last - regions->chanrange[i].first
                 <= 90) &&
                (regions->chanrange[i].last - regions->chanrange[i-1].first
                 <= maxrgnwid))
               {
               /* Join regions. */

                  regions->chanrange[i-1].last = regions->chanrange[i].last;
                  RS_rgn_rmv(regions, i);
                  /* Before i loop counter is incremented, which would
                   * skip over the NEXT region, decrement it here.
                   */
                     i--;
                  /* decrement one more time so a check is made
                   * to combine new region with next one.
                   */
                     i--;
               }
            } /* endif j >= top_chan */
         } /* endif regions within a peakwidth apart */
      } /* end loop i over regions */
}

/*
 * RS_rgn_pad	private routine that pads the ends of regions
 *		using the irw & irch parameters.
 */

static void RS_rgn_pad(const GLWidthEqn *wx, const GLChanRange *chanrange,
                       int irw, int irch, int maxrgnwid, GLRegions *regions)
{
int		i;
double	d_peakwidth;
int		peakwidth;
int		regn_gap, regn_pad, prev_regnwidth, prev_regn_pad, curr_regnwidth;

if (regions->nregions <= 0)
   return;

/* Pad lower end of first region. */

   d_peakwidth = RS_peakwidth(wx, regions->chanrange[0].first);
   peakwidth = (int) (d_peakwidth + .5);

   if ((regions->chanrange[0].last - regions->chanrange[0].first +
        (2 * peakwidth) <= maxrgnwid) &&
       (regions->chanrange[0].first >= chanrange->first))
      regions->chanrange[0].first =
        GPmax(chanrange->first, regions->chanrange[0].first - (2 * peakwidth));

/* Loop through regions starting with second. */
/* Pad upper end of prev and lower end of curr regions. */

   for (i = 1; i < regions->nregions; i++)
      {
      d_peakwidth = RS_peakwidth(wx, regions->chanrange[i].first);
      peakwidth = (int) (d_peakwidth + .5);

      regn_gap = regn_pad = regions->chanrange[i].first -
                              regions->chanrange[i-1].last;

      /* irch is # chans to subtract from ends of regions */
      /* by removing channels from the padding. */

      if (regn_gap >= irch)
         regn_pad = regn_gap - irch;

      /* irw controls region expansion in multiples of peakwidth */
      /* if gap is big enough. */

      if ((regn_gap / peakwidth) > irw)
         regn_pad = irw * peakwidth;

      prev_regnwidth = regions->chanrange[i-1].last -
                         regions->chanrange[i-1].first;

      if ((prev_regnwidth + regn_pad) > maxrgnwid)
         {
         prev_regn_pad = maxrgnwid - prev_regnwidth - 1;
	     prev_regn_pad = GPmax(0, prev_regn_pad);
         }
      else
         prev_regn_pad = regn_pad;

      regions->chanrange[i-1].last += prev_regn_pad;

      curr_regnwidth = regions->chanrange[i].last -
                         regions->chanrange[i].first;

      if ((curr_regnwidth + regn_pad) > maxrgnwid)
         {
         regn_pad = maxrgnwid - curr_regnwidth - 1;
		 regn_pad = GPmax(0, regn_pad);
         }

      regions->chanrange[i].first -= regn_pad;
      }

/* Add pad to upper end of last region. */

   i = regions->nregions - 1;

   d_peakwidth = RS_peakwidth(wx, regions->chanrange[i].last);
   peakwidth = (int) (d_peakwidth + .5);

   regn_gap = (chanrange->last - 5) - regions->chanrange[i].last;
   regn_pad = GPmax(0, regn_gap - irch);

   if (regn_gap / peakwidth > irw)
      regn_pad = irw * peakwidth;

   curr_regnwidth = regions->chanrange[i].last - regions->chanrange[i].first;

   if (curr_regnwidth + regn_pad < RS_MIN_REGN_WIDTH)
      RS_rgn_rmv(regions, i);
   else
      {
      if (curr_regnwidth + regn_pad > maxrgnwid)
         {
         regn_pad = maxrgnwid - curr_regnwidth - 1;
		 regn_pad = GPmax(0, regn_pad);
         }

      regions->chanrange[i].last += regn_pad;
      }
}

/*
 * RS_rgn_rmv	private routine that removes a region from a list.
 */

static void RS_rgn_rmv(GLRegions *regions, int regn_index)
{
int	i, j;

for (i = regn_index, j = i + 1; j < regions->nregions; i++, j = i + 1)
   {
   regions->chanrange[i].first = regions->chanrange[j].first;
   regions->chanrange[i].last = regions->chanrange[j].last;
   }

regions->nregions = i;
}

/*
 * RS_rgnsforpks	private routine that forces regions for
 *			existing peaks.
 */

static void RS_rgnsforpks(const GLWidthEqn *wx, const GLPeakList *peaks,
                          int firstchannel, int nchannels,
                          GLboolean *regn_flag)
{
int		i, j;
double	d_peakwidth;
int		bottom_chan, top_chan;

/* Assure search peaks are included. */

   for (i = 0; i < peaks->npeaks; i++)
      {
      if (peaks->peak[i].channel_valid == GL_TRUE)
         {
         d_peakwidth = RS_peakwidth(wx, peaks->peak[i].channel);

         /* Set regn_flag TRUE in area around peak. */

            bottom_chan = (int) (peaks->peak[i].channel - firstchannel -
				d_peakwidth);
            bottom_chan = GPmax(0, bottom_chan);

            top_chan = (int) (peaks->peak[i].channel - firstchannel +
				d_peakwidth + .5);
            top_chan = GPmin((nchannels - 1), top_chan);

            for (j = bottom_chan; j <= top_chan; j++)
               regn_flag[j] = GL_TRUE;
         }
      }
}

/*
 * RS_store_rgns	private routine that translates an array
 *			of flags into a list of regions.
 */

static GLRtnCode RS_store_rgns(const GLChanRange *chanrange, int firstchannel,
                               int nchannels, GLboolean *regn_flag,
                               GLRegions *regions)
{
int			i, j;
int			bottom_chan, top_chan;
GLboolean	within_region;
GLboolean	exceed_flag;

/*
 * Comments and change by Egger on 9/10/99.
 *
 * In original Gauss VII fortran, the scan for regions to store ranges
 * from first channel in search range + 5
 * to last channel in search range - 5.
 *
 * This potentially misses a region that was flagged for an
 * existing peak by the RS_rgnsforpks() routine because this
 * region could start before the scanning starts.
 *
 * So start the scanning at the beginning of the search range.
 *
 * Egger change on 9/17/2001 - first region could start in first channel.
 */

/* Record start and end of each region. */

   exceed_flag = GL_FALSE;
   within_region = GL_FALSE;
   j = 0;
   bottom_chan = chanrange->first - firstchannel;
   bottom_chan = GPmax(0, bottom_chan);
   top_chan = chanrange->last - firstchannel;
   top_chan = GPmin((nchannels - 1), top_chan);
   for (i = bottom_chan; i <= top_chan; i++)
      {
      if ((within_region == GL_FALSE) && (regn_flag[i] == GL_TRUE))
	     {
         if (regions->listlength <= j)
		    {
			exceed_flag = GL_TRUE;
			break;
			}

         within_region = GL_TRUE;
         regions->chanrange[j].first = i + firstchannel;
	     }
      else if ((within_region == GL_TRUE) && (regn_flag[i] == GL_FALSE))
         {
         within_region = GL_FALSE;
         regions->chanrange[j].last = i - 1 + firstchannel;
         j++;
         }
      }
   regions->nregions = j;

   /* If the last region did not end before the search range,
    * then add one last region that ends with the search range.
    */
      if ((exceed_flag == GL_FALSE) && (within_region == GL_TRUE))
         {
         regions->chanrange[j].last = chanrange->last;
         regions->nregions++;
         }

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

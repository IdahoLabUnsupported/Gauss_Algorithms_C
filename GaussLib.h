/*
 * Copyright 2017 Battelle Energy Alliance
 */

/*
 *  GaussLib.h - contains typedefs and prototypes for the Gauss Algorithms
 *		library
 */

#ifndef GAUSSLIB_H
#define GAUSSLIB_H

/*
 * When building a DLL in Visual C++, declspec() is required.
 * Otherwise, leave this undefined.
 */

#ifdef _WINDOWS
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT
#endif

#define GAUSS_ALGORITHMS_VERSION "3.3"

/*
 * Typedefs for prototype arguments:
 */


/*
 * stuff for boolean logic
 */

#ifndef False
#define	False	0
#define True	1
#endif

   typedef enum
      {
      GL_FALSE = False,
      GL_TRUE  = True
      } GLboolean;


/*
 * definition of fitting range
 */

   typedef struct
      {
      int		first;
      int		last;
      } GLChanRange;


/*
 * Energy equation modes
 */

   typedef enum
      {
      GL_EGY_LINEAR,
      GL_EGY_QUADRATIC
      } GLEgyEqnMode;


/*
 * GLEnergyEqn holds coefficients for the energy equation
 * e(x) = a + b*x + c*x**2
 * and chi squared from corresponding calibration.
 */

   typedef struct
      {
      double		a;
      double		b;
      double		c;
      double		chi_sq;
      GLEgyEqnMode	mode;
      } GLEnergyEqn;


/*
 * Width equation modes
 */

   typedef enum
      {
      GL_WID_LINEAR,
      GL_WID_SQRT
      } GLWidEqnMode;


/*
 * GLWidthEqn holds coefficients for the width equation
 *
 * linear: w(x) = alpha + beta*x
 *
 * sqrt:   w(x) = (alpha + beta*x)**1/2
 *
 * and chi squared from corresponding calibration.
 */

   typedef struct
      {
      double		alpha;
      double		beta;
      double		chi_sq;
      GLWidEqnMode	mode;
      } GLWidthEqn;


/*
 * GLSpectrum holds the counts per channel of a spectrum
 */

   typedef struct
      {
      int		listlength;	/* array size of count & sigcount */
      int		nchannels;
      int		firstchannel;
#if (! defined(GL_LINUX)) && (! defined(GL_MACOSX))
      __int32		*count;		/* for Visual C++ */
#else
      int		*count;		/* list of counts per channel */
#endif
      double		*sigcount;	/* count uncertainty          */
      } GLSpectrum;


/*
 * GLPeakType indicates which data form was used to add/define the peak.
 */

   typedef enum
      {
      GL_PEAK_CHANNEL = 0,
      GL_PEAK_ENERGY = 1
      } GLPeakType;


/*
 * GLPeak is a structure to hold all information about a single peak.
 */

   typedef struct
      {
      GLPeakType	type;
      GLboolean		channel_valid;
      double		channel;
      GLboolean		energy_valid;
      double		energy;
      double		sige;
      GLboolean		fixed_centroid;
      } GLPeak;


/*
 * GLPeakList is a structure to hold a list of peaks.
 */

   typedef struct
      {
      int		listlength;	/* array size of peak */
      int		npeaks;
      GLPeak		*peak;		/* array of peaks */
      } GLPeakList;

/*
 * GLPeakRefinement is a structure to hold the results
 *                  of a peak refinement.
 */

   typedef struct
      {
      double		raw_channel;
      GLChanRange	refine_region;
      double		net_area;
      double		background;
      double		refined_channel;
      } GLPeakRefinement;


/*
 * GLPeakSearchResults is a structure to hold the results
 *                     of a peak search.
 */

   typedef struct
      {
	  GLPeakList	*peaklist;         /* found peaks */
	  GLPeakRefinement	*refinements;		/* array of peak refinements */
	  int			listlength;        /* array size of crossproducts */
	  int			*crossproducts;    /* cross-correlations */
      } GLPeakSearchResults;


/*
 * GLRegions is a structure to hold list of region coordinates.
 *
 * Gauss IX limits the region width to 150 channels.
 */

   typedef struct
      {
      int		listlength;	/* array size of chanrange */
      int		nregions;
      GLChanRange	*chanrange;	/* list of channel ranges */
      } GLRegions;


/*
 * Region search modes
 *
 * RGNSRCH_FORPKS is used when want only regions for existing peaks
 * and no regions without peaks.
 */

   typedef enum
      {
      GL_RGNSRCH_ALL,
      GL_RGNSRCH_FORPKS
      } GLRgnSrchMode;


/*
 * enumerations of legal fit parameter values
 */

   /*
    * legal values for the pkwd_mode parameter which is used to
    * set convergence criteria and determine when to vary peakwidth.
    *
    * See the description of GLCCType for a discussion of setting
    * the convergence criteria.
    *
    * For determining at beginning of each fit cycle whether
    * to vary peakwidth:
    *
    * if pkwd_mode == GL_PKWD_VARIES
    *      if atleast one peak centroid is not fixed, let average
    *      peakwidth vary.
    *      if jth peak is ~511keV and it is not only peak,
    *         let fitinfo.peakinfo[j].addwidth_511 vary.
    *
    * if pkwd_mode == GL_PKWD_FIXED
    *      if jth peak is ~511keV,
    *         let fitinfo.peakinfo[j].addwidth_511 vary.
    */

      typedef enum 
         {
         GL_PKWD_VARIES = 0,
         GL_PKWD_FIXED = 1
         } GLPkwdMode;


   /*
    * legal values for the split_mode parameter which is used to control
    * whether the linear and nonlinear parameters are fitted separately.
    *
    * if split_mode == GL_SPLIT_NONE
    *      a fit cycle calls lmder() only once with all
    *      parameters allowed to vary together.
    *
    * if split_mode == GL_SPLIT_ALLOWED
    *      a fit cycle calls lmder() 5 times (5 passes).
    *
    *      After each of first 4 passes, peakheight is constrained
    *      to be greater than 2.0.
    *
    *      On 1st & 3rd passes, the linear parameters (intercept, slope,
    *	   peak heights, avg_width) are allowed to vary.
    *
    *      On 2nd & 4th passes, the nonlinear parameters (centroids,
    *      peak addwidth_511s) are allowed to vary.
    *
    *      On 5th pass, all parameters are varied together.
    *
    *      The effect of splitting the parameters like this is to
    *      allow the fast linear least squares fitting method to
    *      improve these parameter estimates a great deal before
    *      the nonlinear fitting is done.
    *
    * GL_SPLIT_ALLOWED takes more time and so it is only recommended
    * to investigate the fitting process itself.
    *
    * Gauss IX discourages the use of GL_SPLIT_ALLOWED, forcing the user
    * to use a command line option to invoke this choice.
    */

      typedef enum
         {
         GL_SPLIT_NONE = 0,
         GL_SPLIT_ALLOWED = 1
         } GLSplitMode;


   /*
    * enumeration of types for the cc_type parameter
    *
    * cc_type is used in conjunction with pkwd_mode to choose values
    * for the convergence criteria ftol & xtol, which are inputs to
    * the marquardt algorithm lmder() which is used by fitregn().
    * ftol measures the relative error desired in the sum of squares.
    * xtol measures the relative error desired in the approximate
    * solution.
    *
    * setting ftol & xtol:
    *
    * if   cc_type == GL_CC_LARGER
    * and  pkwd_mode == GL_PKWD_VARIES
    * then ftol = .001
    * and  xtol = .003
    *
    * if   cc_type == GL_CC_LARGER
    * and  pkwd_mode == GL_PKWD_FIXED
    * then ftol = xtol = .01
    *
    * if   cc_type == GL_CC_SMALLER
    * then ftol = xtol = .001
    *
    * if   cc_type == GL_CC_LARGER_INC
    * then ftol = .001
    * and  xtol = .003
    */

      typedef enum
         {
         GL_CC_LARGER     = 0,
         GL_CC_SMALLER    = 1,
         GL_CC_LARGER_INC = 2
         } GLCCType;


/*
 * Fit Parameters
 *
 * ncycle		- maximum number of fit cycles allowed
 * nout			- maximum number of the best fits to be used
 * max_npeaks		- maximum number of peaks allowed in a fit
 *			  (Gauss IX uses '10')
 * pkwd_mode		- control peakwidth (see above)
 * split_mode		- controls splitting of linear and nonlinear
 *			  parameter fits (see above)
 * cc_type		- used to set convergence criteria (see above)
 * max_resid		- fitregn() recycles until the fit residual at
 *			  each channel in the region is less than the
 *			  value of max_resid. (if not stopped sooner)
 *			  (Gauss IX uses default value of 2.0 counts)
 */

   typedef struct
      {
      int		ncycle;		/* hardcoded to 5 in Gauss VII */
      int		nout;		/* hardcoded to 3 in Gauss VII */
      int		max_npeaks;	/* hardcoded to 10 in Gauss VII */
      GLPkwdMode	pkwd_mode;	/* known as iw in Gauss VII */
      GLSplitMode	split_mode;	/* known as ifit in Gauss VII */
      GLCCType		cc_type;	/* known as iconv in Gauss VII */
      float		max_resid;	/* known as sd in Gauss VII */
      } GLFitParms;


/*
 * A record of fit information corresponding to one fit cycle
 *
 * The 'used' information is the information that was passed
 * into the fitregn() algorithm for the fit.
 *
 * The rest of the information in this structure is passed out
 * by the fitregn() algorithm.
 *
 * lmder_info	- the info parameter returned from lmder()
 * chi_sq	- is reduced chi squared of the fit of this region
 * fitinfo	- used to compute summary and curve
 * uncertainty	- used to compute summary
 * back_covarr	- copied from jacobian in nllsqs() - it is the
 *		  covariance of the background's intercept and
 *		  other term (slope or step_height).
 * covarr	- computed from jacobian in nllsqs()
 *
 */

#define GL_COVARR_DIM		3

   typedef struct
      {
      double		height;		/* gaussian height in counts */
      double		centroid;	/* gaussian position in channels */
      double		addwidth_511;	/* width correction for 511keV */
      GLboolean		fixed_centroid;
      } GLPeakInfo;

   typedef struct
      {
      int		listlength;	/* storage length of peakinfo */
      double		intercept;	/* bckgrnd intercept @ chnrng.first */
      double		slope;		/* bckgrnd slope */
      double		step_height;	/* unused - fraction of gaussian ht */
      double		avg_width;	/* fwhm in channels */
      int		npeaks;		/* # of elements in peakinfo array */
      GLPeakInfo	*peakinfo;
      } GLFitInfo;

   /* assume that fitinfo.listlength == uncertainty.listlength */
   /* assume that first dimension of covarr[][] is same listlength */

   typedef struct
      {
      int		cycle_number;
      GLChanRange	used_chanrange;
      GLFitParms	used_parms;
      GLEnergyEqn	used_ex;
      GLWidthEqn	used_wx;
      int		lmder_info;
      double		chi_sq;
      GLFitInfo		fitinfo;
      GLFitInfo		uncertainty;
      double		back_covarr;
      double		(*covarr) [GL_COVARR_DIM]; /*covarr[i][GL_COVARR_DIM]*/
      } GLFitRecord;


/*
 * A list of fit records returned by fitregn()
 *
 * Each record in the list corresponds to one fit cycle.
 */

   typedef struct fitreclist
      {
      GLFitRecord	*record;
      struct fitreclist	*next;
      } GLFitRecList;


/*
 * GLSummary is output from summarize().
 */

   typedef struct
      {
      int	listlength;	/* array size of channel, etc. */
      int	npeaks;
      double	ratio;		/* of summation area to integral area */
      double	*channel;
      double	*sigc;
      double	*height;
      double	*sigh;
      double	*wid;
      double	*sigw;
      double	*area;
      double	*siga;
      double	*energy;
      double	*sige;
      } GLSummary;


/*
 * GLCurve contains coordinates for the fit, background, and residual
 * at each channel in the region.
 *
 * fitcurve is an array of component curves.  There can be up to 'npeaks'
 * of these curves.  and each curve has regn_width*nplots_per_chan coordinates.
 *
 * so listlength should be >= regn_width*nplots_per_chan
 */

   typedef struct
      {
      int		listlength;	/* array size - usually "npoints" */
      GLChanRange	chanrange;
      int		nplots_per_chan;	/* # plotpoints per channel */
      int		npoints;
      int		npeaks;
      double		*x_offset;	/* array of npoints offsets relative */
					/* to start of chanrange */
      double		**fitpeak;	/* array of npeak curve arrays */
      double		*fitcurve;	/* array of npoints curve y coords */
      double		*back;		/* array of npoints back y coords */
      double		*resid;		/* array of residual at each channel */
      } GLCurve;


/*
 * GLEfficiency holds information from the efficiency file and
 * related calculations from spline().
 *
 * en is the table of data abscissas
 * eff is the table of ordinates
 * coeff is the 2D array of spline coefficients computed by spline() (seval()?)
 */

#define GL_MAX_EFFICY		100
#define GL_NUM_COEFFS		3

   typedef struct
      {
      int	neff;				/* INTEGER*4    */
      double	en[GL_MAX_EFFICY];		/* REAL*8 x 100 */
      double	eff[GL_MAX_EFFICY];		/* REAL*8 x 100 */
      double	coeff[GL_MAX_EFFICY][GL_NUM_COEFFS]; /* REAL*8 x 300 */
      } GLEfficiency;


/*
 * GLPostInfo contains
 *
 * efficiency - detector efficiency at each peak
 * intensity  - area/efficiency at each peak
 * sigi       - the intensity uncertainty
 */

   typedef struct
      {
      int	listlength;	/* array size of efficiency, etc. */
      int 	npeaks;
      double	*efficiency;
      double	*intensity;
      double	*sigi;
      } GLPostInfo;


/*
 * Return codes from the procedures in the Gauss Library:
 *
 * GL_SUCCESS		no problems in execution of procedure
 * GL_FAILURE		unspecified problem in execution of procedure
 *
 * GL_BADMALLOC		failure to allocate temporary workspace in memory;
 *			procedure returned without completing.
 *
 * GL_BADCALBLINF	in calibration, error returned by linf_()
 * GL_BADPEAKWD		in any procedure, GL_chan_to_w() returned error
 *			and recovery was not possible.
 *
 * GL_BADCHNRNG		the channel range is illegal or
 *			inconsistent with spectrum
 *
 * GL_BADIRCH		the value of irch passed in is <= 0
 * GL_BADIRW		the value of irw passed in is <= 0
 * GL_BADTHRESH		the value of threshold passed in is < 0
 *
 * GL_OVRLMT		more peaks or regions were found than structure
 *			could hold; the extra peaks or regions are ignored;
 *			the returned list contains the items found up to
 *			the limit.
 * GL_BADNPKS		more peaks provided for fit than fitparms allow
 * GL_BADRGNWD		region search returned one or more regions
 *			wider than user limit
 */

   typedef enum
      {
      GL_SUCCESS,
      GL_FAILURE,
      GL_BADMALLOC,
      GL_BADCALBLINF,
      GL_BADPEAKWD,
      GL_BADCHNRNG,
      GL_BADIRCH,
      GL_BADIRW,
      GL_BADTHRESH,
      GL_OVRLMT,
      GL_BADNPKS,
      GL_BADRGNWD
      } GLRtnCode;


/*
 * Prototypes for procedures in the Gauss Library
 */

#ifdef __cplusplus
extern "C" {
#endif


/*
 * GL_add_chanpeak	adds a peak to the list in terms of channel.
 *			Type will be set to GL_PEAK_CHANNEL, centroid
 *			will be initialized to NOT fixed, and energy
 *			will be left undefined.  If the list does
 *			not have room, GL_OVRLMT is returned.
 *			Otherwise, GL_SUCCESS is returned.
 */

   DLLEXPORT GLRtnCode GL_add_chanpeak(double channel, GLPeakList *peaks);


/*
 * GL_add_egypeak	adds a peak to the list in terms of energy.
 *			Type will be set to GL_PEAK_ENERGY, centroid
 *			will be initialized to fixed, and channel
 *			will be left undefined.  If the list does
 *			not have room, GL_OVRLMT is returned.
 *			Otherwise, GL_SUCCESS is returned.
 */

   DLLEXPORT GLRtnCode GL_add_egypeak(double energy, double sige,
                                      GLPeakList *peaks);


/*
 * GL_add_peak		adds a peak to the list.  If the list does not
 *			have room, GL_OVRLMT is returned.  Otherwise,
 *			GL_SUCCESS is returned.
 */

   DLLEXPORT GLRtnCode GL_add_peak(const GLPeak *peak, GLPeakList *peaks);


/*
 * GL_chan_to_e		converts channel number to energy using the given
 *			energy equation.  The calling routine must provide
 *			space for the answer in 'energy'.
 */

   DLLEXPORT void GL_chan_to_e(const GLEnergyEqn *ex, double channel,
                               double *energy);


/*
 * GL_chan_to_w		calculates the peak width at the indicated channel.
 *			The calling routine must provide space for the
 *			answer in 'width'.
 */

   DLLEXPORT GLRtnCode GL_chan_to_w(const GLWidthEqn *wx, double channel,
                                    double *width);


/*
 * GL_curve		extract curve information from fitdata.
 *			The calling routine must provide space for
 *			the answer in 'curve'.
 *
 *			"nplots_per_chan - 1" is how many coordinates
 *			will be plotted BETWEEN channels. 
 *
 */

   DLLEXPORT GLRtnCode GL_curve(const GLSpectrum *spectrum,
                                const GLFitRecord *fit, int nplots_per_chan,
                                GLCurve *curve);


/*
 * GL_curve_alloc	allocate memory for the curve structure.
 *			listlength is set in the returned structure.
 *
 *			If routine fails, returns NULL.
 */

   DLLEXPORT GLCurve *GL_curve_alloc(int regn_width, int nplots_per_chan,
                                     int npeaks);


/*
 * GL_curve_free	free curve structure memory that was
 *			allocated with GL_curve_alloc().
 */

   DLLEXPORT void GL_curve_free(GLCurve *curve);


/*
 * GL_e_to_chan		converts energy to channel using the given energy
 *			equation.  The calling routine must provide space
 *			for the answer in 'channel'.
 */

   DLLEXPORT GLRtnCode GL_e_to_chan(const GLEnergyEqn *ex, double energy,
                                    double *channel);


/*
 * GL_ecalib	calibrates the energy equation.  channel, energy, and sige
 *		are assumed to be arrays of size count.  If 'weighted'
 *		is true, then 'sige' is used; otherwise, error for each
 *		energy is fixed to be '1'.  The calling routine must
 *		provide space for the answer in 'ex'.
 */

   DLLEXPORT GLRtnCode GL_ecalib(int count, const double *channel,
                                 const double *energy, const double *sige,
                                 GLEgyEqnMode mode, GLboolean weighted,
                                 GLEnergyEqn *ex);


/*
 *  GL_find_regnpks	return a list of peaks that are within the
 *			indicated region.  The peak list must provide
 *			the channel for each peak.  The calling routine
 *			must provide space for the answer in 'pks_in_rgn'.
 */

   DLLEXPORT GLRtnCode GL_find_regnpks(const GLChanRange *region,
                                       const GLPeakList *peaks,
                                       GLPeakList *pks_in_rgn);


/*
 * GL_fitrec_alloc	allocate memory for the fitrecord structure.
 *			Enter max_npeaks for the listlength.
 *			Listlength is set in the returned structure.
 *
 *			If routine fails, returns NULL.
 */

   DLLEXPORT GLFitRecord *GL_fitrec_alloc(int listlength);


/*
 * GL_fitrec_free	free fitrecord structure memory that was
 *			allocated with GL_fitrec_alloc().
 */

   DLLEXPORT void GL_fitrec_free(GLFitRecord *fitrec);


/*
 * GL_fitreclist_free	free structure memory that was allocated and
 *			returned by GL_fitregn() in "fitlist".
 */

   DLLEXPORT void GL_fitreclist_free(GLFitRecList *fitreclist);


/*
 * GL_fitregn		fit the region indicated by chanrange and store
 *			the answer in fitlist.  Currently, the peak list
 *			must contain only those peaks that lie within
 *			the region indicated by chanrange.  The order
 *			does not appear to matter.  The peak list must
 *			also provide the channel for each peak. If the
 *			number of peaks in the list is greater than
 *			that allowed by the fitparms, error code
 *			GL_BADNPKS is returned and no fitting is done.
 *
 *			All memory allocation for the fitlist is done
 *			inside.  If the GLFitRecList structure is
 *			previously allocated, it is freed first.
 *			Otherwise the calling routine must set
 *			*fitlist=NULL before invoking GL_fitregn().
 *
 */

   DLLEXPORT GLRtnCode GL_fitregn(const GLChanRange *chanrange,
                                  const GLSpectrum *spectrum,
                                  const GLPeakList *peaks,
                                  const GLFitParms *fitparms,
                                  const GLEnergyEqn *ex, const GLWidthEqn *wx,
                                  GLFitRecList **fitlist);


/*
 * GL_get_big_resid_peak	returns first peak with big residual
 *                          where big is ten times "count+1"
 */

   DLLEXPORT GLPeak *GL_get_big_resid_peak(const GLSpectrum *spectrum,
                                           const GLFitRecord *fit);


/*
 * GL_get_neg_peak			returns first peak with negative height
 */

   DLLEXPORT GLPeak *GL_get_neg_peak(const GLFitRecord *fit);


/*
 * GL_get_outside_peak		returns first peak whose centroid
 *                          is located outside the region
 */

   DLLEXPORT GLPeak *GL_get_outside_peak(const GLFitRecord *fit);


/*
 * GL_get_posneg_peakpair		returns first pos-neg peak pair
 */

   DLLEXPORT GLPeakList *GL_get_posneg_peakpair(const GLFitRecord *fit);


/*
 * GL_get_version	public routine to return the version string
 *                  for Gauss Algorithms.
 */

   DLLEXPORT char *GL_get_version();


/*
 * GL_peak_results_alloc  public routine to allocate memory
 *                        for a GLPeakSearchResults structure.
 */

   DLLEXPORT GLPeakSearchResults *GL_peak_results_alloc(int peak_listlength,
	                                                    int spectrum_nchannels);


/*
 * GL_peak_results_free		public routine to free memory for a
 *                          GLPeakSearchResults structure.
 */
   
   DLLEXPORT void GL_peak_results_free(GLPeakSearchResults *results);


/*
 * GL_peaks_alloc	allocate memory for a GLPeakList structure.
 *			Listlength is set in the returned structure.
 */

   DLLEXPORT GLPeakList *GL_peaks_alloc(int listlength);


/*
 * GL_peaks_free	free GLPeakList structure memory that was
 *			allocated with GL_peaks_alloc().
 */

   DLLEXPORT void GL_peaks_free(GLPeakList *peaks);


/*
 * GL_peaksearch	searches for peaks in the spectrum.  Threshold
 *			controls the pruning out of insignificant peaks.
 *			Large values increase the pruning.
 *
 *			for low search sensitivity, use  threshold = 20
 *			for med search sensitivity, use  threshold = 10
 *			for high search sensitivity, use threshold =  5
 *
 *			The refine argument indicates whether the peak locations
 *			are refined with a summing fit.
 *
 *			The channel of each peak is returned in 'peaks'.
 *			The type of each peak is GL_PEAK_CHANNEL, and
 *			'fixed_centroid' and energy_valid are both set to
 *			GL_FALSE for each peak.
 *
 *			The calling routine must provide space for the
 *			answer in 'results'.
 */

   DLLEXPORT GLRtnCode GL_peaksearch(const GLChanRange *chanrange,
                                     const GLWidthEqn *wx, int threshold,
                                     const GLSpectrum *spectrum,
                                     GLboolean refine,
                                     GLPeakSearchResults *results);


/*
 * GL_postprocess	using efficiency file and summary, compute
 *			efficiencies, intensities, and intensity
 *			uncertainties and store in postinfo.  The
 *			calling routine must provide space for the
 *			answer in 'info'.
 *
 *			info is array of nsummaries structures.
 */

   DLLEXPORT GLRtnCode GL_postprocess(const GLEfficiency *efficy,
                                      const GLSummary *summary,
                                      GLPostInfo *info);


/*
 *  GL_prune_rqdpks	remove any required peaks that are too close
 *			to any found peaks.  The peak channels must
 *			be provided in the peak lists.  The calling routine
 *			must provide space for the answer in 'new_rqd'.
 */

   DLLEXPORT GLRtnCode GL_prune_rqdpks(const GLWidthEqn *wx,
                                       const GLPeakList *searchpks,
                                       const GLPeakList *curr_rqd,
                                       GLPeakList *new_rqd);


/*
 * GL_regions_alloc	allocate memory for a GLRegions structure.
 *			Listlength is set in the returned structure.
 */

   DLLEXPORT GLRegions *GL_regions_alloc(int listlength);


/*
 * GL_regions_free	free GLRegions structure memory that was
 *			allocated with GL_regions_alloc().
 */

   DLLEXPORT void GL_regions_free(GLRegions *regions);


/*
 * GL_regnsearch	searches for regions in the spectrum.  Threshold
 *			controls the pruning out of insignificant regions.
 *			Large values increase the pruning.
 *
 *			for low search sensitivity, use  threshold = 3
 *			for med search sensitivity, use  threshold = 2
 *			for high search sensitivity, use threshold = 1
 *
 *			'irw' and 'irch' control padding the ends of
 *			the found regions.  The pad starts with
 *			the gap between regions, is decremented by irch,
 *			and is increased to 'irw * peakwidth'.
 *
 *			recommended starting value for  irw = 3
 *			recommended starting value for irch = 2
 *
 *			The peak list must include the channel for each peak.
 *
 *			maxrgnwid limits the size of each region (number
 *			of channels in a region).  A good value to use
 *			is maxrgnwid = 150.
 *
 *			The calling routine must provide space for the
 *			answer in 'regions'.
 */

   DLLEXPORT GLRtnCode GL_regnsearch(const GLChanRange *chanrange,
                                     const GLWidthEqn *wx, double threshold,
                                     int irw, int irch,
                                     const GLSpectrum *spectrum,
                                     const GLPeakList *peaks,
                                     GLRgnSrchMode mode, int maxrgnwid,
                                     GLRegions *regions);


/*
 * GL_set_sigcount	Using the channel counts provided in a spectral
 *			structure (spec->count), calculate the count
 *			uncertainty (spec->sigcount) and return it in
 *			the spectral structure.  The calling routine
 *			must provide space for the answer in
 *			'spec->sigcount'.
 *
 *			Use small count correction for counts <= 10
 *			based on GW Phillips, NIM 153 (1978), p. 449.
 */

   DLLEXPORT void GL_set_sigcount(GLSpectrum *spec);


/*
 * GL_spline		Use members en and eff of efficiency to
 *			compute and fill member coeff.
 */

   DLLEXPORT void GL_spline(GLEfficiency *efficiency);


/*
 * GL_summ_alloc	allocate memory for the summary structure.
 *			Enter npeaks for the listlength.
 *			Listlength is set in the returned structure.
 *
 *			If routine fails, returns NULL.
 */

   DLLEXPORT GLSummary *GL_summ_alloc(int listlength);


/*
 * GL_summ_free		free summary structure memory that was
 *			allocated with GL_summ_alloc().
 */

   DLLEXPORT void GL_summ_free(GLSummary *summary);


/*
 * GL_summarize		extract summary information from the fit data.
 *			The calling routine must provide the space for
 *			the answer in 'summary'.
 */

   DLLEXPORT GLRtnCode GL_summarize(const GLSpectrum *spectrum,
                                    const GLFitRecord *fit, GLSummary *summary);


/*
 * GL_update_peaklist	Update a peaklist with the specified energy
 *			calibration.  For each peak of type GL_PEAK_CHANNEL,
 *			the energy value is updated.  For each peak of type
 *			GL_PEAK_ENERGY, the channel value is updated.
 *			If ex == NULL, then the updated values are
 *			flagged as invalid instead.
 */

   DLLEXPORT void GL_update_peaklist(const GLEnergyEqn *ex, GLPeakList *peaks);


/*
 * GL_wcalib	calibrates the width equation.  channel, wid, and sigw
 *		are assumed to be arrays of length 'count'.  If
 *		'weighted' is true, then error for each width is
 *		calculated as sigw * wid * 2; otherwise, error for
 *		each width is fixed to be '1'.  The calling routine
 *		must provide space for the answer in 'wx'.
 */

   DLLEXPORT GLRtnCode GL_wcalib(int count, const double *channel,
                                 const double *wid, const double *sigw,
                                 GLWidEqnMode mode, GLboolean weighted,
                                 GLWidthEqn *wx);


#ifdef __cplusplus
}
#endif


#endif  /* GAUSSLIB_H */

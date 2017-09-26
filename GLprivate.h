/*
 * Copyright 2017 Battelle Energy Alliance
 */

/*
 * gaussprv.h - contains prototypes for utilities private to the Gauss
 *		Algorithms library.
 */

#ifndef GAUSSPRV_H
#define GAUSSPRV_H

#define	GP_511PK_THRESHOLD	.5

/*
 * I got these macros from the now obsolete /usr/include/macros.h
 * Do not use anything too fancy for arguments...  such as incrementing.
 *
 * I recommend that arguments a and b be enclosed in parentheses.
 */

#define GPmax(a,b) 	(a<b ? b : a)
#define GPmin(a,b) 	(a>b ? b : a)


/*
 * private subroutine return codes
 */

typedef enum
   {
   GP_CYCLE_DONE = 0,
   GP_CYCLE_DLT = -1,
   GP_CYCLE_ADD = 1,
   GP_CYCLE_CONTINUE = 2
   } GPCycleType;


/*
 * new structure for holding flags previously held in vector "nop".
 * These flags control which items in GLFitInfo structure are allowed
 * to vary in the fitting.  They are also used to indicate which
 * "varied values" are to be copied to/from vectors.
 */

   typedef struct
      {
      GLboolean		height;
      GLboolean		centroid;
      GLboolean		addwidth_511;
      } GPPeakVary;

   typedef struct
      {
      int		listlength;	/* storage length of peakvary */
      GLboolean		intercept;
      GLboolean		slope;
      GLboolean		step_height;
      GLboolean		avg_width;
      int		npeaks;
      GPPeakVary	*peakvary;
      } GPFitVary;


/*
 * dpmpar_	used by MINPACK lmder() and its subroutines
 */

   double dpmpar_(int *choice);


/*
 * GP_calc_functions	used by fcn() to set fvec[].
 */

   void GP_calc_functions(const GLChanRange *chanrange,
                          const GLSpectrum *spectrum, const GLFitInfo *fitinfo,
                          double *resid);


/*
 * GP_calc_jacobian	used by fcn() to set fjac[][].
 *			Currently, only malloc errors detected.
 */

   GLRtnCode GP_calc_jacobian(int chanrange_first, int regn_width,
                              const GLSpectrum *spectrum, GPFitVary *fitvary,
                              const GLFitInfo *fitinfo, double *fjac,
                              int ldfjac);


/*
 * GP_checkfit		determine whether to add or delete a peak & recycle.
 */

   GPCycleType GP_checkfit(const GLSpectrum *spectrum, GLFitRecord *fit,
                           GPFitVary *fitvary, const GLCurve *curve,
                           int *valid_npeaks, GLboolean *arg_peak511,
                           double init_peakwidth);


/*
 * GP_cycle	use lmder_() to fit with current peaklist.  The
 *		number of times lmder_() is called depends on
 *		the split_mode parameter in the fit record.
 */

   GPCycleType GP_cycle(const GLSpectrum *spectrum, GPFitVary *fitvary,
                        double xtol, double ftol, GLFitRecord *fit,
                        GLCurve *curve);


/*
 * GP_info_count	return count of flags set to GL_TRUE.
 */

   int GP_info_count(GPFitVary *fitvary);


/*
 * GP_info_to_vector	using fitvary flags, copy from fitinfo to vector.
 */

   void GP_info_to_vector(GPFitVary *fitvary, GLFitInfo *fitinfo,
                          double *vector);


/*
 * lmder	answer stored in ?
 *		this is the stuff from MINPACK (marquardt algorithm)
 *
 *		fjac is expected to be a 2 dimensional fortran array
 *		fjac[ldfjac][*]
 */

   void lmder_(void (* fcn) (), int *m, int *n, double *x, double *fvec,
               double *fjac, int *ldfjac, double *ftol, double *xtol,
               double *gtol, int *maxfev, double *diag, int *mode,
               double *factor, int *nprint, int *info, int *nfev,
               int *njev, int *ipvt, double *qtf, double *wa1,
               double *wa2, double *wa3, double *wa4);


/*
 * GP_vector_to_info	using fitvary flags, copy vector to fitinfo.
 */

   void GP_vector_to_info(GPFitVary *fitvary, double *vector,
                          GLFitInfo *fitinfo);


#endif  /* GAUSSPRV_H */

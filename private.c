/*
 * Copyright 2017 Battelle Energy Alliance
 */

/*
 *  private.c - contains miscellaneous private utilities for the Gauss
 *		Algorithms library.
 */

#include "GaussLib.h"
#include "GLprivate.h"

/*
 * GP_info_to_vector	global private routine that uses the GPFitVary
 *			flags to copy information from the GLFitInfo
 *			structure to 'vector'.
 */

void GP_info_to_vector(GPFitVary *fitvary, GLFitInfo *info, double *vector)
{
int	i, j;

/* currently assuming enough storage in vector */

i = 0;

if (fitvary->intercept == GL_TRUE)
   vector[i++] = info->intercept;

if (fitvary->slope == GL_TRUE)
   vector[i++] = info->slope;

if (fitvary->step_height == GL_TRUE)
   vector[i++] = info->step_height;

if (fitvary->avg_width == GL_TRUE)
   vector[i++] = info->avg_width;

for (j = 0; j < info->npeaks; j++)
   {
   if (fitvary->peakvary[j].height == GL_TRUE)
      vector[i++] = info->peakinfo[j].height;

   if (fitvary->peakvary[j].centroid == GL_TRUE)
      vector[i++] = info->peakinfo[j].centroid;

   if (fitvary->peakvary[j].addwidth_511 == GL_TRUE)
      vector[i++] = info->peakinfo[j].addwidth_511;
   }
}

/*
 * GP_vector_to_info	global private routine that uses the GPFitVary
 *			flags to copy information from 'vector' to the
 *			GLFitInfo structure.
 */

void GP_vector_to_info(GPFitVary *fitvary, double *vector, GLFitInfo *info)
{
int	i, j;

/* currently assuming enough storage in info */

i = 0;

if (fitvary->intercept == GL_TRUE)
   info->intercept = vector[i++];

if (fitvary->slope == GL_TRUE)
   info->slope = vector[i++];

if (fitvary->step_height == GL_TRUE)
   info->step_height = vector[i++];

if (fitvary->avg_width == GL_TRUE)
   info->avg_width = vector[i++];

info->npeaks = fitvary->npeaks;

for (j = 0; j < info->npeaks; j++)
   {
   if (fitvary->peakvary[j].height == GL_TRUE)
      info->peakinfo[j].height = vector[i++];

   if (fitvary->peakvary[j].centroid == GL_TRUE)
      info->peakinfo[j].centroid = vector[i++];

   if (fitvary->peakvary[j].addwidth_511 == GL_TRUE)
      info->peakinfo[j].addwidth_511 = vector[i++];
   }
}

/*
 * GP_info_count	global private routine that returns the count
 *			of true flags in the GPFitVary structure.
 */

int GP_info_count(GPFitVary *fitvary)
{
int	i, j;

i = 0;

if (fitvary->intercept == GL_TRUE)
   i++;

if (fitvary->slope == GL_TRUE)
   i++;

if (fitvary->step_height == GL_TRUE)
   i++;

if (fitvary->avg_width == GL_TRUE)
   i++;

for (j = 0; j < fitvary->npeaks; j++)
   {
   if (fitvary->peakvary[j].height == GL_TRUE)
      i++;

   if (fitvary->peakvary[j].centroid == GL_TRUE)
      i++;

   if (fitvary->peakvary[j].addwidth_511 == GL_TRUE)
      i++;
   }

return(i);
}

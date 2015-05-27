/* 
 * mkNrrd : This code save in nrrd file data from FTLE or FSLE 3D ascii file
 *          using teem lib (http://teem.sourceforge.net/)
 * Author : P. Monroy
 * Organization: IFISC
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <teem/nrrd.h>
#include <teem/air.h>

char *info = ("generate a nrrd file from a data in ASCII format.");

int main(int argc, const char *argv[]) 
{
  const char *me;
  char *err, *out, *in;
  int size[3], xi, yi, zi;
  hestOpt *hopt;
  hestParm *hparm;
  airArray *mop;
  double min[3], max[3], spacing[3], *data;
  Nrrd *nout;
  FILE *finput;

  me = argv[0];
  mop = airMopNew();
  hparm = hestParmNew();
  hopt = NULL;
  airMopAdd(mop, hparm, (airMopper)hestParmFree, airMopAlways);
  hestOptAdd(&hopt, "i", "filename", airTypeString, 1, 1, &in, NULL,
             "input data file to read data in ASCII format");
  hestOptAdd(&hopt, "s", "sx sy sz", airTypeInt, 3, 3, size, NULL,
             "dimensions of input volume");
  hestOptAdd(&hopt, "min", "x y z", airTypeDouble, 3, 3, min, NULL,
             "lower bounding corner of volume");
  hestOptAdd(&hopt, "max", "x y z", airTypeDouble, 3, 3, max, NULL,
             "upper bounding corner of volume");
  hestOptAdd(&hopt, "o", "filename", airTypeString, 1, 1, &out,"out.nrrd",
             "file to write output nrrd to");
  hestParseOrDie(hopt, argc-1, argv+1, hparm,
                 me, info, AIR_TRUE, AIR_TRUE, AIR_TRUE);
  airMopAdd(mop, hopt, (airMopper)hestOptFree, airMopAlways);
  airMopAdd(mop, hopt, (airMopper)hestParseFree, airMopAlways);

  nout = nrrdNew();
  airMopAdd(mop, nout, (airMopper)nrrdNuke, airMopAlways);
  if (nrrdAlloc_va(nout, nrrdTypeDouble, 3,
                   AIR_CAST(size_t, size[0]),
                   AIR_CAST(size_t, size[1]),
                   AIR_CAST(size_t, size[2]))) {
    airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: problem allocating volume:\n%s\n", me, err);
    airMopError(mop); 
    return 1;
  }

  data = (double *)nout->data;
  if((finput=fopen(in,"r"))==NULL)
    {
      fprintf(stderr, "%s: The file %s does not exist.\n", me, in);
      airMopError(mop); 
      return 1;
    }
  /*
   ******** AIR_AFFINE(i,x,I,o,O)
   **
   ** given intervals [i,I], [o,O] and a value x which may or may not be
   ** inside [i,I], return the value y such that y stands in the same
   ** relationship to [o,O] that x does with [i,I].  Or:
   **
   **    y - o         x - i
   **   -------   =   -------
   **    O - o         I - i
   **
   ** It is the callers responsibility to make sure I-i and O-o are
   ** both non-zero.  Strictly speaking, real problems arise only when
   ** when I-i is zero: division by zero generates either NaN or infinity
   **
   ** NOTE that "x" is evaluated only once (which makes this more useful),
   ** as is "I" and "O" (usually not so important); "i" and "o" are each
   ** evaluated twice
   */
  for (zi=0; zi<size[2]; zi++) 
    {
      for (yi=0; yi<size[1]; yi++) 
	{
	  for (xi=0; xi<size[0]; xi++) 
	    {
	      if(fscanf(finput,"%lf", data)==EOF)
		{
		  fprintf(stderr, "%s: The file %s has not enough data.\n", me, in);
		  airMopError(mop); 
		  return 1; 
		}
	      data += 1;
	    }
	}
    }
  

  nrrdAxisInfoSet_va(nout, nrrdAxisInfoCenter,
		     1, 1, 1);
  nrrdAxisInfoSet_va(nout,nrrdAxisInfoMin, 
		     min[0], min[1], min[2]);
  nrrdAxisInfoSet_va(nout,nrrdAxisInfoMax, 
		     max[0], max[1], max[2]);

  spacing[0] = AIR_DELTA(1.0, 1.0,AIR_CAST(double, size[0]), min[0], max[0]);
  spacing[1] = AIR_DELTA(1.0, 1.0,AIR_CAST(double, size[1]), min[1], max[1]);
  spacing[2] = AIR_DELTA(1.0, 1.0,AIR_CAST(double, size[2]), min[2], max[2]);
  /*
   ******** AIR_DELTA(i,x,I,o,O)
   **
   ** given intervals [i,I] and [o,O], calculates the number y such that
   ** a change of x within [i,I] is proportional to a change of y within
   ** [o,O].  Or:
   **
   **      y             x
   **   -------   =   -------
   **    O - o         I - i
   **
   ** It is the callers responsibility to make sure I-i and O-o are
   ** both non-zero
   **
   ** NOTE that all arguments are evaluated only once
   */
  nrrdAxisInfoSet_va(nout,nrrdAxisInfoSpacing, 
		     spacing[0], spacing[1], spacing[2]);

  if (nrrdSave(out, nout, NULL)) 
    {
      airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
      fprintf(stderr, "%s: problem saving output:\n%s\n", me, err);
      airMopError(mop); return 1;
    }

  airMopOkay(mop);
  exit(0);
}

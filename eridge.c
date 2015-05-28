/* 
 * eridge : This code extract ridges from FTLE or FSLE 3D data (in nrrd file) 
 *          using teem lib (http://teem.sourceforge.net/)
 * Author : P. Monroy
 * Organization: IFISC
 *
 */

#include <teem/seek.h>
#include <teem/air.h>

char *info = ("RIDGE SURFACE EXTRACTION");

int main(int argc, const char *argv[]) 
{
  /* COMMAND LINE ARGUMETS */
  const char *me;  // Executable name
  double strength; // Strength
  double scaling[3]; // Scaling, amount by which to up/down-sample on each spatial axis
  Nrrd *nin;       // Input nrrd file
  char *outvtk;    // File output, vtk format
  hestOpt *hopt=NULL;

  me = argv[0];
  hestOptAdd(&hopt, "s", "strength", airTypeDouble, 1, 1, &strength, NULL,
             "strength"); 
  hestOptAdd(&hopt, "c", "scaling", airTypeDouble, 3, 3, scaling, "1 1 1",
             "amount by which to up/down-sample on each spatial axis");
  hestOptAdd(&hopt, "i", "input nrrd", airTypeOther, 1, 1, &nin, NULL, 
	     "input nrrd file to analyze", 
	     NULL, NULL, nrrdHestNrrd);
  hestOptAdd(&hopt, "o", "output vtk", airTypeString, 1, 1, &outvtk, "out.vtk",
             "output vtk file to save ridges into");
  hestParseOrDie(hopt, argc-1, argv+1, NULL, me, info, AIR_TRUE, AIR_TRUE, AIR_TRUE);

  airArray *mop = airMopNew(); // For memory management
  airMopAdd(mop, nin, (airMopper)nrrdNuke, airMopAlways);
  airMopAdd(mop, hopt, (airMopper)hestOptFree, airMopAlways);
  // Often the examples in teem have this command, but here it causes a annoying memory error: 
  // airMopAdd(mop, hopt, (airMopper)hestParseFree, airMopAlways); 

  /* NRRD: info and compute the maxlength (it is needed to rescaled ridges to original nrrd size) */
  double length[nin->dim];
  double maxlength;
  length[0]=nin->axis[0].spacing * ((double) (nin->axis[0].size-1));
  maxlength=length[0]; 
  int i; // Loop index
  for(i=1; i<nin->dim; i++)
    {
      length[i]=nin->axis[i].spacing * ((double) (nin->axis[i].size-1));
      if(maxlength<length[i])
	maxlength=length[i];
    }

  printf("----------Data Info from input Nrrd----------\n");
  printf("Space info:\n");
  printf(" Dim: %u\n",nin->dim);
  printf("Data info:\n");
  printf(" Type: %d\n",nin->type);
  printf(" Content: %s\n",nin->content);  
  printf("Axis info:\n");
  for(i=0; i<nin->dim; i++)
    printf(" size[%d]=%lu\n",i, nin->axis[i].size);
  for(i=0; i<nin->dim; i++)
    printf(" spacing[%d]: %f\n",i, nin->axis[i].spacing);
  for(i=0; i<nin->dim; i++)
    printf(" min[%d], max[%d]: %f %f\n", i, i, nin->axis[i].min, nin->axis[i].max);
  for(i=0; i<nin->dim; i++)
      printf(" length[%d]=%lf\n", i,length[i]);
  printf(" max. length: %lf\n", maxlength);
  printf("---------------------------------------------\n");

  /* GAGE:
   * it is used to make convolution-based measurements in your input. We will tell it to use a second-order
   * continuous cubic B-spline as the interpolation kernel and to measure all quantities needed for crease 
   * extraction, including Hessian eigenvectors and -values
   */
  gageContext *gctx=NULL;
  gagePerVolume *pvl;
  double kparm[3];
  gctx = gageContextNew();
  airMopAdd(mop, gctx, (airMopper)gageContextNix, airMopAlways);
  ELL_3V_SET(kparm, 1.0, 1.0, 0.0); 
  /* 
   * kparm[0] = 1.0 (Kernel scaling or scale parameter, in units of samples)
   * (kparm[1],kparm[2])=(B,C)=(1,0) -> Kernel parameters 
   */

  if (!(pvl = gagePerVolumeNew(gctx, nin, gageKindScl))
      || gagePerVolumeAttach(gctx, pvl)
      || gageKernelSet(gctx, gageKernel00, nrrdKernelBCCubic, kparm) // Values -> Uniform cubic B-spline (B,C)=(1,0) 
      || gageKernelSet(gctx, gageKernel11, nrrdKernelBCCubicD, kparm) // First Deriv. -> Uniform cubic B-spline 
      || gageKernelSet(gctx, gageKernel22, nrrdKernelBCCubicDD, kparm) // Second Deriv. -> Uniform cubic B-spline
      || gageQueryItemOn(gctx, pvl, gageSclValue) // Measure -> Scalar value 
      || gageQueryItemOn(gctx, pvl, gageSclNormal) // Measure -> Gradient vector, normalized
      || gageQueryItemOn(gctx, pvl, gageSclHessEval) // Measures -> Hessian's eigenvalues
      || gageQueryItemOn(gctx, pvl, gageSclHessEvec) // Measures -> Hessian's eigenvectors
      || gageQueryItemOn(gctx, pvl, gageSclHessEval2) // Measures -> Hessian's 3rd eigenvalue = Ridge strength
      || gageUpdate(gctx)) 
    {
      char *err;
      airMopAdd(mop, err = biffGetDone(GAGE), airFree, airMopAlways);
      fprintf(stderr, "ERROR while setting up Gage:\n%s\n", err);
      airMopError(mop); return 1;
    }

  /* SEEK: Set up the extraction itself */

  seekContext *sctx;
  sctx = seekContextNew();
  airMopAdd(mop, sctx, (airMopper)seekContextNix, airMopAlways);

  limnPolyData *pld;
  pld = limnPolyDataNew(); 
  airMopAdd(mop, pld, (airMopper)limnPolyDataNix, airMopAlways);

  seekVerboseSet(sctx, 10); // set verbose level
  int E=0;
  size_t samples[3];
  ELL_3V_SET(samples,
           2*nin->axis[0].size, 2*nin->axis[1].size, 2*nin->axis[2].size);
  if (!E) E |= seekSamplesSet(sctx, samples);
  if (!E) E |= seekDataSet(sctx, NULL, gctx, 0);
  if (!E) E |= seekItemGradientSet(sctx, gageSclGradVec);
  if (!E) E |= seekItemEigensystemSet(sctx, gageSclHessEval, gageSclHessEvec);
  if (!E) E |= seekItemHessSet(sctx, gageSclHessian);
  if (!E) E |= seekItemNormalSet(sctx, gageSclNormal);
  if (!E) E |= seekStrengthUseSet(sctx, AIR_TRUE);
  if (!E) E |= seekStrengthSet(sctx, -1, strength);
  if (!E) E |= seekItemStrengthSet(sctx, gageSclHessEval2);
  if (!E) E |= seekNormalsFindSet(sctx, AIR_TRUE);
  if (!E) E |= seekTypeSet(sctx, seekTypeRidgeSurfaceT);
  if (!E) E |= seekUpdate(sctx);
  if (!E) E |= seekExtract(sctx, pld);
  if (E) 
    {
      char *err;
      airMopAdd(mop, err = biffGetDone(SEEK), airFree, airMopAlways);
      fprintf(stderr, "ERROR during surface extraction:\n%s\n", err);
      airMopError(mop); return 1;
    }
  fprintf(stderr, " extraction time = %g\n", sctx->time);

  /* Post-filtering */
  Nrrd *nval = nrrdNew();
  airMopAdd(mop, nval, (airMopper) nrrdNuke, airMopAlways);
  if (1==seekVertexStrength(nval, sctx, pld)) 
    {
      char *err;
      airMopAdd(mop, err = biffGetDone(SEEK), airFree, airMopAlways);
      fprintf(stderr, "ERROR during surface probing:\n%s\n", err);
      airMopError(mop); 
      return 1;
    }
  limnPolyDataClip(pld, nval, sctx->strength);

  /* Classify the output by primitives */
  if (limnPolyDataCCFind(pld))
    {
      char *err;
      err = biffGetDone(LIMN);
      fprintf(stderr, "%s: trouble sorting:\n%s", me, err);
      free(err);
    }

  /*  WRITE VTK FILE FROM PLD: Triangle Soup with connectivity in cell atribute*/

  double alpha, beta[3]; // alpha and beta is needed in order to rescale pld to original nrrd grid
  alpha = 0.5 * maxlength; 
  beta[0] = 0.5 * length[0];
  beta[1] = 0.5 * length[1];
  beta[2] = 0.5 * length[2];

  FILE *filevtk;
  filevtk = fopen(outvtk,"w");
  fprintf(filevtk,"# vtk DataFile Version 3.0\n");
  fprintf(filevtk,"vtk output\n");
  fprintf(filevtk,"ASCII\n");
  fprintf(filevtk,"DATASET POLYDATA\n");
  fprintf(filevtk,"POINTS %d float\n", pld->xyzwNum);

  double *xyz;
  xyz = AIR_MALLOC(3*pld->xyzwNum, double);// TO DO: include the free command of xyz in airmop 

  int j; // Loop index
  for(i=0, j=0; i < 4*pld->xyzwNum; i+=4, j+=3)
    {
      xyz[j] = alpha * pld->xyzw[i] + beta[0] + nin->axis[0].min;
      xyz[j+1] = alpha * pld->xyzw[i+1] + beta[1] + nin->axis[1].min;
      xyz[j+2] = alpha * pld->xyzw[i+2] + beta[2] +  nin->axis[2].min;

      fprintf(filevtk,"%10.7f ", xyz[j]);//print x
      fprintf(filevtk,"%10.7f ", xyz[j+1]);//print y
      fprintf(filevtk,"%10.7f\n", xyz[j+2]);//print z
    }
  fprintf(filevtk,"POLYGONS %d %d\n", pld->indxNum/3, (pld->indxNum/3)*4);
  for(i=0; i < pld->indxNum; i+=3)
    {
      fprintf(filevtk,"3 %u", pld->indx[i]);//print indx
      fprintf(filevtk, " %u", pld->indx[i+1]);//print indx
      fprintf(filevtk, " %u\n", pld->indx[i+2]);//print indx
    }
  fprintf(filevtk,"CELL_DATA %d\n", pld->indxNum/3);
  fprintf(filevtk,"SCALARS Connectivity unsigned_int\n");
  fprintf(filevtk,"LOOKUP_TABLE default\n");
  for(i=0; i < pld->primNum; i++)
    {
      for(j=0; j<pld->icnt[i]/3; j++)
	fprintf(filevtk,"%u\n",pld->icnt[i]/3);
    }
  
  fprintf(filevtk,"POINT_DATA %d\n", pld->xyzwNum);
  fprintf(filevtk,"SCALARS Scalar_Value float\n");
  fprintf(filevtk,"LOOKUP_TABLE default\n");

  const double *val;
  double xi,yi,zi;

  val = gageAnswerPointer(gctx, pvl, gageSclValue);
  for(i=0; i < 3*pld->xyzwNum; i+=3)
    {
      xi = (xyz[i] - nin->axis[0].min)/nin->axis[0].spacing;
      if (xi < 0.0)
	xi = (double) 0;
      else if (xi > (double) (nin->axis[0].size-1))
	xi = (double) (nin->axis[0].size-1);

      yi = (xyz[i+1] - nin->axis[1].min)/nin->axis[1].spacing;
      if (yi < 0.0)
	yi = (double) 0;
      else if (yi > (double) (nin->axis[1].size-1))
	yi = (double) (nin->axis[1].size-1);

      zi = (xyz[i+2] - nin->axis[2].min)/nin->axis[2].spacing;        
      if (zi < 0.0)
	zi = (double) 0;
      else if (zi > (double) (nin->axis[2].size-1))
	zi = (double) (nin->axis[2].size-1);

      if (gageProbe(gctx, xi, yi, zi)) 
	{
	  fprintf(stderr, "%s: trouble:\n(%d) %s\n", me, gctx->errNum, gctx->errStr);
	  airMopError(mop); 
	  return 1;
	}
      fprintf(filevtk,"%g\n",*val);
    }
  fclose(filevtk);
  
  /* LOG: Prints some values of PLD (TO DO: give it format and print more things in the command line)*/
  printf("----------Data Info from output *pld----------\n");
  printf("Number of points = %u\n", pld->xyzwNum); 
  printf("Number of triangles = %u\n", pld->indxNum/3);
  printf("Number of primitives = %u\n", pld->primNum); 
  printf("----------------------------------------------\n");

  /* Freeing memory  */
  airFree(xyz); // TO DO: include the free command of xyz in airmop 
  airMopOkay(mop);
  exit(0);
}

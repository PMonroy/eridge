/*
  TO COMPILE:
  gcc -O2 eridge.c -o eridge -lteem -lm

*/

#include <teem/seek.h>
#include <teem/air.h>

char *info = ("RIDGE SURFACE EXTRACTION");

int main(int argc, const char *argv[]) 
{
  int i,j; // Loop index
  limnPolyData *pld;
  gageContext *gctx=NULL;
  gagePerVolume *pvl;
  seekContext *sctx;

  double kparm[3];




  FILE *filevtk;


  double alpha, beta[3];

  /* COMMAND LINE ARGUMETS */
  const char *me;  // Executable name
  double strength; // Strength
  Nrrd *nin;       // Input nrrd file
  char *outvtk;    // File output, vtk format
  hestOpt *hopt=NULL;

  me = argv[0];
  hestOptAdd(&hopt, "s", "strength", airTypeDouble, 1, 1, &strength, NULL,
             "strength");  
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
  gctx = gageContextNew();
  airMopAdd(mop, gctx, (airMopper)gageContextNix, airMopAlways);
  ELL_3V_SET(kparm, 1.0, 1.0, 0.0);
  if (!(pvl = gagePerVolumeNew(gctx, nin, gageKindScl))
      || gagePerVolumeAttach(gctx, pvl)
      || gageKernelSet(gctx, gageKernel00, nrrdKernelBCCubic, kparm)
      || gageKernelSet(gctx, gageKernel11, nrrdKernelBCCubicD, kparm)
      || gageKernelSet(gctx, gageKernel22, nrrdKernelBCCubicDD, kparm)
      || gageQueryItemOn(gctx, pvl, gageSclNormal) // Gradient vector, normalized
      || gageQueryItemOn(gctx, pvl, gageSclHessEval) // Hessian's eigenvalues
      || gageQueryItemOn(gctx, pvl, gageSclHessEvec) // Hessian's eigenvectors
      || gageQueryItemOn(gctx, pvl, gageSclHessEval2) // Hessian's 3rd eigenvalue=ridge strength
      || gageUpdate(gctx)) 
    {
      char *err;
      airMopAdd(mop, err = biffGetDone(GAGE), airFree, airMopAlways);
      fprintf(stderr, "ERROR while setting up Gage:\n%s\n", err);
      airMopError(mop); return 1;
    }

  /* SEEK: Set up the extraction itself */
  sctx = seekContextNew();
  airMopAdd(mop, sctx, (airMopper)seekContextNix, airMopAlways);
  pld = limnPolyDataNew(); 
  airMopAdd(mop, pld, (airMopper)limnPolyDataNix, airMopAlways);

  seekVerboseSet(sctx, 10); // set verbose level
  int E=0;
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
  alpha = 0.5 * maxlength; // alpha and beta is needed in order to rescale pld to original nrrd grid
  beta[0] = 0.5 * length[0];
  beta[1] = 0.5 * length[1];
  beta[2] = 0.5 * length[2];

  filevtk = fopen(outvtk,"w");
  fprintf(filevtk,"# vtk DataFile Version 3.0\n");
  fprintf(filevtk,"vtk output\n");
  fprintf(filevtk,"ASCII\n");
  fprintf(filevtk,"DATASET POLYDATA\n");
  fprintf(filevtk,"POINTS %d float\n", pld->xyzwNum);
  for(i=0; i < 4*pld->xyzwNum; i+=4)
    {
      fprintf(filevtk,"%10.7f ", alpha * pld->xyzw[i] + beta[0]);//print x
      fprintf(filevtk,"%10.7f ", alpha * pld->xyzw[i+1] + beta[1]);//print y
      fprintf(filevtk,"%10.7f\n",alpha * pld->xyzw[i+2] + beta[2]);//print z
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
  fclose(filevtk);
  
  /* LOG: Prints some values of PLD (TO DO: give it format and print more things in the command line)*/
  printf("number of primitives = %u\n", pld->primNum); 
  printf("number of index = %u\n", pld->indxNum);
  printf("number of coordenates = %u\n", 4*pld->xyzwNum); 

  airMopOkay(mop);
  exit(0);
}

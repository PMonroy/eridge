/*
  TO COMPILE:
  gcc -O2 eridge.c -o eridge -lteem -lm

*/

#include <teem/seek.h>
#include <teem/air.h>

char *info = ("RIDGE SURFACE EXTRACTION");

int main(int argc, const char *argv[]) 
{
  const char *me;
  hestOpt *hopt=NULL;
  airArray *mop = airMopNew(); // For memory management


  Nrrd *nin; // Input nrrd file
  double strength; // Strength
  char *outvtk; // File output, vtk format

  limnPolyData *pld;
  gageContext *gctx=NULL;
  gagePerVolume *pvl;
  seekContext *sctx;

  double kparm[3];
  char *err;

  int E;

  FILE *filevtk;


  double alpha, beta[3];

  /* COMMAND LINE ARGUMETS */
  me = argv[0];
  hestOptAdd(&hopt, "s", "strength", airTypeDouble, 1, 1, &strength, NULL,
             "strength");  
  hestOptAdd(&hopt, "i", "input nrrd", airTypeOther, 1, 1, &nin, NULL, 
	     "input nrrd file to analyze", 
	     NULL, NULL, nrrdHestNrrd);
  hestOptAdd(&hopt, "o", "output vtk", airTypeString, 1, 1, &outvtk, "out.vtk",
             "output vtk file to save ridges into");
  hestParseOrDie(hopt, argc-1, argv+1, NULL, me, info, AIR_TRUE, AIR_TRUE, AIR_TRUE);


  /* TEST NRRD*/
  double length[nin->dim];
  double maxlength;
  printf("Nrrd: nin");
  printf("type=%d\n", nin->type);
  printf("dim=%u\n", nin->dim);
  printf("content=%s\n", nin->content);
  int i;
  for(i=0; i<nin->dim; i++)
    {     
      printf("axis[%d].size=%lu\n",i, nin->axis[i].size);
      printf("axis[%d].spacing=%lf\n",i, nin->axis[i].spacing);
      printf("axis[%d].thickness=%lf\n",i, nin->axis[i].thickness);
      printf("axis[%d].min=%lf\n",i, nin->axis[i].min);
      printf("axis[%d].max=%lf\n",i, nin->axis[i].max);
      printf("axis[%d].center=%d\n",i, nin->axis[i].center);
      printf("axis[%d].kind=%d\n",i, nin->axis[i].kind);
      printf("axis[%d].label=%s\n",i, nin->axis[i].label);
      printf("axis[%d].units=%s\n",i, nin->axis[i].units);
      length[i]=nin->axis[i].spacing * ((double) (nin->axis[i].size-1));
      printf("length[%d]=%lf\n", i,length[i]);
    }

  maxlength= length[0];

  for(i=1; i<nin->dim; i++)
    {
      if(maxlength<length[i])
	maxlength=length[i];
    }

  printf("maxlength=%lf\n", maxlength);
  airMopDebug(mop);
  
  /* INITIALIZATION */
  pld = limnPolyDataNew();// is not added a comand to the airMop stack to pld 

  sctx = seekContextNew();
  airMopAdd(mop, sctx, (airMopper)seekContextNix, airMopAlways);

  gctx = gageContextNew();
  airMopAdd(mop, gctx, (airMopper)gageContextNix, airMopAlways);

  airMopAdd(mop, nin, (airMopper)nrrdNuke, airMopAlways);
  airMopDebug(mop);

  /* GAGE */
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
      airMopAdd(mop, err = biffGetDone(GAGE), airFree, airMopAlways);
      fprintf(stderr, "ERROR while setting up Gage:\n%s\n", err);
      airMopError(mop); return 1;
    }

  /* EXTRACTION: Set up the extraction itself */

  seekVerboseSet(sctx, 10);
  
  E = 0;
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
      airMopAdd(mop, err = biffGetDone(SEEK), airFree, airMopAlways);
      fprintf(stderr, "ERROR during surface probing:\n%s\n", err);
      airMopError(mop); return 1;
    }
  limnPolyDataClip(pld, nval, sctx->strength);

  airMopOkay(mop);// memory is freed except pld

  /*  WRITE VTK FILE FROM PLD*/
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
  fclose(filevtk);


  


  exit(0);
}

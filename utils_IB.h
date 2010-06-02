// Information Bottleneck code 
// incl. annealing
// Susanne Still, 2004.
// modified to use hash tables (histogram_t) 
// for large dimensional data. Lisa Miller 2010

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "utils_ps.h"

void ini_pygc(int NCLUSTnow, int YLENGTH, int NCLUST, histogram_t *pygc);
void ini_pcgx(int NCLUST, int XLENGTH, histogram_t *pcgx);
void marg_px(int XLENGTH, int YLENGTH, histogram_t *pxy, double px[XLENGTH]);
void marg_py(int XLENGTH, int YLENGTH, histogram_t *pxy, double py[YLENGTH]);
void conditional(int XLENGTH, int YLENGTH, histogram_t *pxy, double px[XLENGTH], histogram_t *pygx);
double I_XY(int XLENGTH, int YLENGTH, histogram_t *pxy, double px[XLENGTH], double py[YLENGTH]);
double I_YgX(int XLENGTH, int YLENGTH, double **pygc, double px[XLENGTH], double py[YLENGTH]);
double H_X(int XLENGTH, double px[XLENGTH]);
//bool Bayes(int XLENGTH, int YLENGTH, double **pygx, double px[XLENGTH],double **pxgy);
bool Bayes(int XLENGTH, int YLENGTH, histogram_t *pygx, double px[XLENGTH], histogram_t *pxgy);
bool P_c(int XLENGTH, int NCLUST, double **pcx, double px[XLENGTH], double pc[NCLUST]);
//==================================================================================================
// Input: XLENGTH, YLENGTH, NCLUST, p(c|x), p(x), p(y|x)
// Output: p(y|c)
// Mstep(XLENGTH, YLENGTH, NCLUST, pygx, px, pcgx, pygc);
//bool Mstep(int XLENGTH, int YLENGTH, int NCLUST, histogram_t *pygx, double px[XLENGTH], double **pcgx, double **pygc);
bool Mstep(int XLENGTH, int YLENGTH, int NCLUST, histogram_t *pygx, double px[XLENGTH], histogram_t *pcgx, histogram_t *pygc);
//==================================================================================================
// rescaled DKL
// Input: XLENGTH, YLENGTH, NCLUST, p(y|x), p(y|c)
// Output: DKL[ p(y|x) || p(y|c)]
//oid DKL_rescaled(int XLENGTH, int YLENGTH, int NCLUST,  histogram_t *pygx, histogram_t *pygc, histogram_t *DKL);
//===========================================================================================================
void DKL_Prod(int XLENGTH, int YLENGTH, int NCLUST,  histogram_t *pygx, histogram_t *pygc, histogram_t *DKL, double beta);
//==================================================================================================
// DKL (of two vectors containing probability densities)
// Input: XLENGTH, YLENGTH, NCLUST, p(y|x), p(y|c)
// Output: DKL[ p(y|x) || p(y|c)]
double D_KL(int XLENGTH, double px1[XLENGTH], double px2[XLENGTH]);
//==================================================================================================
// Input: XLENGTH, NCLUST, DKL[ p(y|x) || p(y|c)], p(c), beta
// Output: p(c|x), p(c)
bool Estep(int XLENGTH, int NCLUST, double beta, histogram_t *DKL, histogram_t *pc, histogram_t *pcgx);
//==================================================================================================
// perturb cluster centers
void perturb_pygc(double pert, int YLENGTH, int NCLUST, double **pygc);
//==================================================================================================
// compute self-consistent IB equations 
// NEEDS NORMALIZED HISTOGRAM AS INPUT:  p(x,y)
// Input: XLENGTH, YLENGTH,  NCLUST, p(x,y), p(c|x) -- initialization
// Output: p(c|x), p(y|c)
// F = EM(QPN, 2, NCLUSTnow, NCLUSTmax, beta, P_history_future, P_future_given_history, P_nextstate, clustercenter, pc);
double EM(int XLENGTH, int YLENGTH, int NCLUST, int NCLUSTmax, double beta,
	  histogram_t *pxy,
	  histogram_t *histo,
	  histogram_t *FULLpcgx,
	  histogram_t *FULLpygc, 
	  histogram_t *FULLpc,
	  double infocurve[4],
	  bool printme);

//=================================================================================================
// split clusters
bool split(int NCLUSTnow[], int NCLUSTmax, int YLENGTH, double **pygc);
//=================================================================================================
// merge clusters
int merge(double beta, int NCLUSTnow, int NCLUSTmax, int YLENGTH, double **pygc);
//=================================================================================================
//check if cluster assignments are (almost) deterministic
bool det(int NCLUSTnow, int NCLUSTmax, int XLENGTH, double **pcgx);
//=================================================================================================
//check if clusters get used
bool used(double usemin, int NCLUST, double pc[]);
//==================================================================================================
//==================================================================================================
double anneal(int XLENGTH, int YLENGTH, int NCLUST, double beta,
		      histogram_t *pxy,
		      histogram_t *histo,
		      histogram_t *pcgx,
		      histogram_t *pygc, 
		      histogram_t *pc,
		      double infocurve[4], 
		      double annealingrate, bool plotme);
//==================================================================================================


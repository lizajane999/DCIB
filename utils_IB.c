// Information Bottleneck code 
// incl. annealing
// Susanne Still, 2004.
// Modified to use histo2 hash tables for matrices
// Lisa Miller 2010

#include "utils.h"
#include "utils_IB.h"

//==================================================================================================
// make (random) input cluster centers p(y|c)
void ini_pygc(int NCLUSTnow, int YLENGTH, int NCLUST, histogram_t *pygc)
{ 
  int i, k;
  double v,sum, checksum;
  bool debug = true;
  
  printf("\n==================================\nini_pygc:INITIALIZING %d CLUSTER CENTERS\n==================================\n", NCLUSTnow);
 
 
  v = 1/(double)YLENGTH;
  for(k=0; k< NCLUST; k++){
    for(i = 0; i< YLENGTH; i++){
      histogram_set(pygc,i,k,v); 
    }
  }
    
  for(k=0; k< NCLUSTnow; k++){
    sum = 0;
    for(i = 0; i< YLENGTH; i++){
      v = drand_uniform(0.1, 0.2) + EPS;
      histogram_set(pygc,i,k,v);
      sum += v;
    }
    checksum = 0;
    for(i = 0; i< YLENGTH; i++){
      histogram_divide(pygc,i,k,sum);
      v = histogram_get(pygc,i,k);
      if(debug){ printf("pygc[%d][%d] = %f ", i, k, v);}
      checksum += v;
     }
    if(debug){printf("\n");}
    if((double)checksum > (double)1 + (double)EPS || (double)checksum < (double)1 - (double)EPS){
      printf("ini_pygc: error normalizing p(y|c). checksum = %f\n", checksum);
    }
   }
  if(debug){
    printf("ini_pygc: output: NCLUST = %d NCLUSTnow = %d INITIAL centers:\n", NCLUSTnow, NCLUST);  
    for(k=0; k< NCLUST; k++){
      for(i = 0; i< YLENGTH; i++){
	printf("%f\t",histogram_get(pygc,i,k));
      }
      printf("\n");	
    }
  }
}

//==================================================================================================
// pic input cluster centers p(y|c) from input distribution p(y|x) at random. and then perturb a bit.
void pic_pygc(int XLENGTH, int YLENGTH, int NCLUST, double pygc[YLENGTH][NCLUST], double pygx[YLENGTH][XLENGTH])
{ 
  int i, j, k, l;
  int list[100];
  bool redo = false;
  double vec1[YLENGTH];
  double vec2[YLENGTH];
  double D = 100;
  double Dmin = 100000;

  // compute min DKL between any two example distributions
  for(k=0; k<NCLUST; k++){
    for(j = 0; j< YLENGTH; j++){
      vec1[j] = pygx[j][k];
    }
    for(l=0; l<k; l++){
      for(j = 0; j< YLENGTH; j++){
	vec2[j] = pygc[j][l];
      }
      D =  D_KL(YLENGTH, vec1, vec2);
      if(Dmin > D){
	Dmin = D;
      }
    }
  }

  k = 0;
  while(k< NCLUST){
    i = irand_uniform(0, XLENGTH - 1);
    redo = false;
    if(k>0){ 
      for(l=0; l<k; l++){
	if(list[l]==i){ // did this cluster get picked before?
	  redo = true;
	}	    
      }
      if(redo == false){ // is the center similar to one that got picked before?
	for(j = 0; j< YLENGTH; j++){
	  vec1[j] = pygx[j][i];
	}
	for(l=0; l<k; l++){
	  for(j = 0; j< YLENGTH; j++){
	    vec2[j] = pygc[j][l];
	  }
	  D =  D_KL(YLENGTH, vec1, vec2);
	  if(D <= Dmin+Dmin/5){ // D_kl should not smaller than smallest D_kl +20% (arbitrary threshold value)  
	    redo = true;
	  }
	}
      }
    }
    if(redo == false){
      list[k] = i;
      for(j = 0; j< YLENGTH; j++){
	pygc[j][k] = pygx[j][i];
      }
      k++;
    }
  }
 
  // perturb_pygc(0.001, YLENGTH, NCLUST, pygc);
}

//==================================================================================================
// make (random) input assignment p(c|x) 
void ini_pcgx(int NCLUST, int XLENGTH, histogram_t *pcgx)
{ 
  int i, k;
  double sum, v;

   printf("ini_pcgx: INITIAL assignment rule:\n");
   pcgx = histogram_create(NCLUST*XLENGTH);
  // random matrix
  for(i=0; i< XLENGTH; i++){ 
    sum = 0;    
    for(k=0; k< NCLUST; k++){
      v = drand_uniform(0, 1) + EPS; 
      histogram_set(pcgx,k,i,v);
      sum += v;
    }
    for(k=0;k<NCLUST;k++){
      histogram_divide(pcgx,k,i,sum);
      printf("%f\t", histogram_get(pcgx,k,i));
    }
    printf("\n");	
  }
}

//==================================================================================================
// marginalize
// input: p(x,y), XLENGTH, YLENGTH,
// output: p(x)
// marg_px(QPN, 2, P_history_future, P_history);
void marg_px(int XLENGTH, int YLENGTH, histogram_t *pxy, double px[XLENGTH])
{
  int i, j;
  double p,checksum = 0;

    printf("\n HEY this is marg_px:\n"); 
/*    printf("\n Printing P(x,y):\n"); */
/*    print_mat(XLENGTH, YLENGTH, pxy); */
  for(i = 0; i < XLENGTH; i++){
    px[i] = 0;
    for(j = 0; j < YLENGTH; j++){
      if((p = histogram_get(pxy,i,j))){
	px[i] += p; 
	//printf("pxy[%d][%d] = %f\n",i,j,p);
      }
    }
   // printf("px[%d] = %f\n",i,px[i]);
    checksum += px[i]; 
  }
  if(checksum > 1 + 0.00001 || checksum < 1 - 0.00001){
    printf("marg_px: error normalizing p(x). sum = %f\n", checksum);
  } 
    //       else{
      // 	printf("%f ",z);
      //       }
  }

//==================================================================================================
// marginalize
// input: p(x,y), XLENGTH, YLENGTH,
// output: p(y)
void marg_py(int XLENGTH, int YLENGTH,  histogram_t *pxy, double py[YLENGTH])
{
  int i, j;
  double p, checksum = 0.0;
  
  for(j = 0; j < YLENGTH; j++){
    py[j] = 0;
    for(i = 0; i < XLENGTH; i++){
      if((p = histogram_get(pxy,i,j))){
	py[j] += p;
      }
    }
    checksum += py[j]; 
  }

  if(checksum > 1 + 0.00001 || checksum < 1 - 0.00001){
    printf("marg_py: error normalizing p(y). sum = %f\n", checksum);
  }
}

//==================================================================================================
// conditional probability
// Input p(x,y); p(x)
// Output: p(y|x)
void conditional(int XLENGTH, int YLENGTH, histogram_t *pxy, double px[XLENGTH], histogram_t *pygx)
{
  int i, j;
  double p,checksum = 0;
 
  for(j = 0; j < YLENGTH; j++){
    for(i = 0; i < XLENGTH; i++){
      if((p = histogram_get(pxy,i,j))){
	histogram_set(pygx, i, j, (p/px[i]));
      }
    }
  }

  for(i = 0; i < XLENGTH; i++){
    checksum = 0;
    for(j = 0; j < YLENGTH; j++){
      if((p = histogram_get(pygx,i,j))){
	checksum += p; 
      }
    }
    if(checksum > (1 + 0.00001) || checksum < (1 - 0.00001)){
      printf("conditional: error normalizing p(y|x). sum = %f\n", checksum);
    } 
  }

}

//==================================================================================================
// compute total mutual information I(X;Y)
// Input = XLENGTH, YLENGTH, pxy, px, py
// Output = I(X;Y) in bits

double I_XY(int XLENGTH, int YLENGTH, histogram_t *pxy, double px[XLENGTH], double py[YLENGTH])
{ 
  double Ixy = 0;
  double p;
  int i, j;
  for(i=0;i<XLENGTH;i++)
    {
      for(j=0;j<YLENGTH;j++)
	{
	  if((p = histogram_get(pxy,i,j)) && px[i]>0 && py[j]>0)
	    Ixy += p * log( p / ( px[i]*py[j] ) );
	}
    }  
  Ixy /= log(2.0);
  return Ixy;
}

//==================================================================================================
// compute total mutual information I(X;Y)
// Input = XLENGTH, YLENGTH, py_given_x, px, py
// Output = I(X;Y) in bits

double I_YgX(int XLENGTH, int YLENGTH, double **pygx, double px[XLENGTH], double py[YLENGTH])
{ 
  double Ixy = 0;
  int i, j;
  for(i=0;i<XLENGTH;i++)
  {
    for(j=0;j<YLENGTH;j++)
    {
      if(pygx[j][i]>0 && px[i]>0 && py[j]>0)
	Ixy += pygx[j][i] * px[i] * log(pygx[j][i] / py[j]);
    }
  }  
  Ixy /= log(2.0);
  return Ixy;
}

//==================================================================================================
// compute total entropy H(X)
// Input = XLENGTH, px
// Output = H(X) in bits
double H_X(int XLENGTH, double px[XLENGTH])
{  
  double Hx = 0;
  int i;
  for(i=0;i<XLENGTH;i++){
    if(px[i] > 0)
      Hx -= px[i]*log(px[i]);
  }
  Hx /= log(2.0);
  return Hx;
}
//*************TODO*********************************
//==================================================================================================
// BAYES RULE
// Input XLENGTH, YLENGTH, p(y|x), p(x)
// Output p(x|y)
bool Bayes(int XLENGTH, int YLENGTH, double **pygx, double px[XLENGTH],double **pxgy)
{
  int i, j;
  double sum = 0;
  double checksum = 0;
  bool t = false;
  
  printf("Bayes:\n");
  // p(x|y) = p(y|x)*p(x) / (sum_x p(y|x)*p(x) )
  for(j=0;j<YLENGTH;j++){ 
    sum = 0;
 //   printf("C = %d\t", j);
    for(i=0;i<XLENGTH;i++){
      //printf("X = %d : px[i] = %f : pcgx = %f \t", i, px[i], pygx[i][j]);
      pxgy[i][j] = pygx[j][i] * px[i];
      sum += pygx[j][i] * px[i];
    }
    //printf("\n");
    if(sum < EPS){
      printf("Bayes: Error in input distribution: sum_x p(y|x)p(x) = %f at y = %d", sum, j);
      t = true;
    }
    checksum = 0;
    for(i=0;i<XLENGTH;i++){
      pxgy[i][j] /= sum;
      checksum += pxgy[i][j];
    }
    if(checksum > 1 + 0.000000001 || checksum < 1 - 0.000000001){
      printf("Bayes: P(x|y) not normalized properly. checksum =  %f\n",checksum);
    }
  }
  return t;
}

//*************TODO******************************
//==================================================================================================
// compute p(c) from p(c|x) and p(x)
bool P_c(int XLENGTH, int NCLUST, double **pcx, double px[XLENGTH], double pc[NCLUST])
{
  int i,k;
  double checksum = 0.0;  
  bool t = false;

  for(k=0; k<NCLUST; k++){ 
    pc[k] = 0;
    for(i=0; i<XLENGTH; i++){
      pc[k] += pcx[k][i] * px[i];
    }
    if(pc[k] == 0){
      t = true;
   //   printf("P_c: pc[%d] = %f\n", k, pc[k]);
    }
    checksum += pc[k]; 
  }

  if(checksum > 1.0001 || checksum < 1 - 0.0001){
    printf("P_c: error normalizing p(c). sum = %f\n", checksum);
    t = true;
  }

  return t;
}

//*****************TODO********************************************
//==================================================================================================
// M-STEP in IB
// Input: XLENGTH, YLENGTH, NCLUST, p(c|x), p(x), p(y|x)
// Output: p(y|c)
// Mstep(XLENGTH, YLENGTH, NCLUST, pygx, px, pcgx, pygc);
bool Mstep(int XLENGTH, int YLENGTH, int NCLUST, histogram_t *pygx, double px[XLENGTH], double **pcgx, double **pygc)
{ 
 // double pxgc[XLENGTH][NCLUST];
  int i,j,k;
  double p,sum = 0;
  bool alarm = false;
  bool t = false;
  bool debug = false;
  // if the probability of one cluster is absolutely zero, this will result in a NaN (+inf) in p(x|c). 
  // to avoind that, have to ad small noise to p(c|x) in Estep or EM !!

  double** pxgc;
  //allocate pointer memory for first dimension
  pxgc = (double**)malloc(XLENGTH*sizeof(double));
  if(NULL == pxgc){free(pxgc); printf("Memory allocation failed while allocating for pxgc[].\n"); exit(-1);}
  /*allocate memory for second dimension */
  for(i = 0; i < XLENGTH;i++)
  {
    pxgc[i] = (double *) malloc( NCLUST * sizeof(double) );
    if(NULL == pxgc[i]){free(pxgc[i]); printf("Memory allocation failed while allocating for matrix[x][].\n"); exit(-1);}
  }
  // compute p(x|c)
  alarm = Bayes(XLENGTH, NCLUST, pcgx, px, pxgc);
  if (alarm == true){
    printf("Mstep: Bayes returns alarm. p(c|x) or p(x) must have zeros.\n");
  }
  printf("Mstep: \n"); 
  if(debug){
    printf("\n Mstep: p(x|c)\n");     
    for(i=0; i<XLENGTH; i++){
      for(k=0; k<NCLUST; k++){
	printf("%f \t", pxgc[i][k]); 
      }
      printf("\n");
    }
  }
  // compute p(y|c) = \sum_x p(y|x) p(x|c)
  for(j=0; j<YLENGTH; j++){
    for(k=0; k<NCLUST; k++){ 
      pygc[j][k] = 0;
    }
  } 
  t = false;
  // printf("\np(y|c); p(y|x); p(x|c):\n");    
  for(k=0; k<NCLUST; k++){ 
    sum = 0;    
    for(j=0; j<YLENGTH; j++){
      for(i=0; i<XLENGTH; i++){
	if((p = histogram_get(pygx,i,j))){
	  pygc[j][k] += (p * pxgc[i][k]);
	}
	if(debug) printf("p(x|c)[%d][%d] = %f\n",i,k, pxgc[i][k]);
      }
      if(debug) printf("\n p(y|c)[%d][%d] = %f \n\n", j,k,pygc[j][k]);
      sum += pygc[j][k];
      if(pygc[j][k] < EPS){ // check if p(y|c) = 0
	t = true;
	//printf("Mstep: alarm. pygc[%d][%d] = %f\n", j, k, pygc[j][k]); 
      }
    }
    // printf("\n");
    if(sum > 1 + 0.0001 || sum < 1 - 0.0001){
      printf("Mstep: P(y|c) not normalized properly. sum =  %f\n",sum);
    }
  }
  if (t == true){
    alarm = add_EPS_pxgy(YLENGTH, NCLUST, pygc);
  }
  if (alarm == true){
    printf("Mstep: add_EPS_pxgy returns alarm. p(y|c) has zeros.\n");
  }
  for(i = 0; i < XLENGTH; i++)
    free(pxgc[i]);
  free(pxgc);
  
  return alarm;
}

//replace dkl rescaled with susanna's below !!!!!!
//==================================================================================================
// rescaled DKL
// Input: XLENGTH, YLENGTH, NCLUST, p(y|x), p(y|c)
// Output: DKL[ p(y|x) || p(y|c)]
void DKL_rescaled(int XLENGTH, int YLENGTH, int NCLUST,  histogram_t *pygx, histogram_t *pygc, histogram_t *DKL)
{
 int i,j,k;
 double pyx, pyc, v;
 double mini[XLENGTH];
 bool debug = true;
 // reasoning: define m(x) = min_c DKL[p(y|x)||p(y|c)]
 // p(c|x) = ( p(c)/Z(x,beta) ) * exp{-beta*DKL[p(y|x)||p(y|c)]} 
 //        = ( p(c)/Z(x,beta) ) * exp{-beta*(DKL[p(y|x)||p(y|c)]-m(x))} * exp{beta*m(x)} 
 //        = ( p(c)/Z'(x,beta) ) * exp{-beta*(DKL[p(y|x)||p(y|c)]-m(x))}
 // with Z'(x,beta) = Z(x,beta) * exp{-beta*m(x)}
 // and Z'(x,beta) = sum_c p(c) * exp{-beta*(DKL[p(y|x)||p(y|c)]-m(x))} which is what
 // we use later in "Estep" to calculate p(c|x). 

 for(i=0;i<XLENGTH;i++){
   for(k=0;k<NCLUST;k++){
     //DKL[i][k] = 0;
     for(j=0;j<YLENGTH;j++){
       //if both of these are not zeros
       if((pyx = histogram_get(pygx,i,j)) && (pyc = histogram_get(pygc,j,k))){
	// printf("p= %f, pygc[%d][%d] = %f\n", p, j, k, pygc[j][k]);
	v = pyx * log( pyx / pyc );
	histogram_add(DKL,i,k,v);
	// DKL[i][k] += p * log( p / pygc[j][k] );
       }
     }
   }
 }
 if(debug)
 {
    printf("\n DKL_rescaled: DKL = \n");
    histogram_print(DKL);
 }
    
 // rescale
 for(i=0;i<XLENGTH;i++){
   mini[i] = 100000000;
   for(k=0;k<NCLUST;k++){
     v = histogram_get(DKL,i,k);//if DKL is not set or zero, v = 0
     if(mini[i] > v){
       mini[i] = v;
     }
   }
 }
 if(debug)
 {
    printf(" DKL_rescaled: minimum = \n");
    print_vec(XLENGTH, mini);
 } 
 for(i=0;i<XLENGTH;i++){
   for(k=0;k<NCLUST;k++){
     if((v = histogram_get(DKL,i,k))){
       histogram_set(DKL,i,k,(v-mini[i]));
     //DKL[i][k] -= mini[i];
     }
   }
 }
 if(debug)
 {
  printf("\n DKL_rescaled: rescaled DKL = \n");
  histogram_print(DKL);
 }
 
}
//**************** TODO ************** send beta to this function *****************!!!!!!!!!!!!!
//======================== NEW FUNCTION
//==================================================================================================
// rescaled DKL
// Input: XLENGTH, YLENGTH, NCLUST, p(y|x), p(y|c)
// Output: DKL[ p(y|x) || p(y|c)]
void DKL_Prod(int XLENGTH, int YLENGTH, int NCLUST,  histogram_t *pygx, histogram_t *pygc, histogram_t *DKL, double beta)
{
   int i,j,k;
   double p,q,v;
   bool debug = true;
 
  // reasoning: p(c|x) ~ p_(c) * e ^ {-1/T \sum_{y} (p(y|x)* log p(y|x)) - \sum_{y}(p(y|x) * log p(y|c))}
  // 			  \sum_{y} (p(y|x)* log p(y|x))  is constant --> Z (norm)
  //        	p(c|x) ~ p(c) * e ^ (1/T \sum_{y} (p(y|x) * log p(y|c)))
  //        		~ p(c) \prod_{y} e^(1/T * p(y|x)* log p(y|c))
  // 			~ p(c) \prod_{y} (p(y|c))^{p(y|x)/T}
  // tf: if p(y|x) = 0, stuff inside prod = 1 --> no need to calc 
  // we use later in "Estep" to calculate p(c|x). 
  
  for(i=0;i<XLENGTH;i++){
    for(k=0;k<NCLUST;k++){
      for(j=0;j<YLENGTH;j++){
	if((p = histogram_get(pygx,i,j)) && (q = histogram_get(pygc,j,k))){//only make non zero ones
	  // printf("p= %f, pygc[%d][%d] = %f\n", p, j, k, pygc[j][k]);
	  p = p*beta;//divide p(y|c) by temp
	  if((v = histogram_get(DKL,i,k)))
	  {//have to check for zero or not set yet
	    histogram_set(DKL,i,k,(v * pow(q,p)));
	  }else{//don't multiply by zero
	    histogram_set(DKL,i,k, (pow(q,p)));
	  }
	}
      }
    }
  }
  if(debug)
  {
    printf("\n DKL_Prod: DKL = \n");
    print_mat(XLENGTH, NCLUST, DKL);
  }
  
}

//==================================================================================================
// DKL (of two vectors containing probability densities)
// Input: XLENGTH, YLENGTH, NCLUST, p(y|x), p(y|c)
// Output: DKL[ p(y|x) || p(y|c)]

double D_KL(int XLENGTH, double px1[XLENGTH], double px2[XLENGTH])
{
 int i;
 double Dkl = 0;
 bool debug = false;

 if(debug) printf("\nD_KL:\n");
 for(i=0;i<XLENGTH;i++){
   if(px1[i] > 0 && px2[i] > 0 ){
     Dkl += px1[i] * log( px1[i] / px2[i] );
     if(debug)printf("px1 = %f; px2 = %f; Dkl = %f\n",px1[i],px2[i],Dkl);
   }
 }
 return Dkl;
}

//==================================================================================================
// E-STEP in IB
// Input: XLENGTH, NCLUST, DKL[ p(y|x) || p(y|c)], p(c), beta
// Output: p(c|x)
bool Estep(int XLENGTH, int NCLUST, double beta, histogram_t *DKL, histogram_t *pc, histogram_t *pcgx)
{
  int i,k;
  double p,q,v;
  double sumup[XLENGTH];
  bool t = false;
  bool alarm = false;
  bool debug = true;
  bool useDKLReScale = false;

   printf("\nEstep: computing pcx... exp(-beta*DKL):\n");
   //*******************************TODO********************************************************
//   // compute p(c|x)*Z = p(c)*e^(-beta*DKL(c))
  if(!useDKLReScale){
    for(i=0;i<XLENGTH;i++){
      sumup[i] = 0.0;
      if(debug) {printf("x = %d\n",i);}
	for (k=0; k<NCLUST; k++){
	  /************ Using Product ******************************/ 
	  //check for empties
	  // this is p(c) * (prod_{y} (p(y|c))^p(y|x))^beta
	  if((p = histogram_get(pc,k,0)) && (q = histogram_get(DKL,i,k)))
	  {
	    v = p * q;//underflowing!!
	    if(v){
	      histogram_set(pcgx,k,i,v);
	      sumup[i] += v;
	    }
	/// old: pcgx[k][i] = pc[k] * exp( - beta * DKL[i][k]); // + EPS;
	    if(debug){ 
	      printf("p(c = %d) = %e; DKL(x = %d,c = %d) = %e; p(c|x) = %e\n", k, p, i,k,q,v);
	    }
//**************** ALLOWING ZEROS HERE! ************************************
	    // 	    if(v < EPS){
// 	      t = true;
// 	      printf("setting t to true!\n");
	  } 
      }//for k loop
      if(debug) printf("\n");
    }//i-loop
  }
  else{ //using DKL-rescaled
    // printf("\nEstep: computing pcx... exp(-beta*DKL):\n");
    // compute p(c|x)*Z = p(c)*e^(-beta*DKL(c)) 
    for(i=0;i<XLENGTH;i++){
      sumup[i] = 0.0;
       printf("x = %d\n",i);
      for (k=0; k<NCLUST; k++){
	if((p = histogram_get(pc,k,0)) && (q = histogram_get(DKL,i,k)))
	{
	  // OLD WAY pcgx[k][i] = pc[k] * exp( - beta * DKL[i][k]); // + EPS;
	  v = exp( - beta * q);
	  histogram_set(pcgx,k,i,(p * v));
	  printf("p(c = %d) = %e; exponent = %e; p(c|x) = %e\n", k, p, v, (p*v));
	  if(v < EPS){
	    t = true;
	    printf("setting t to true!\n");
	  }
	  sumup[i] += (p*v);
	}
      }
      // printf("\n");
    }//i-loop
  }
  
   if(debug)
   {
    printf("\nEstep: BEFORE NORMALIZATION:\n");
   // print_hist_inv(NCLUST,XLENGTH,pcgx);
   }
//   if(t == true){
//     alarm = add_EPS_pxgy(NCLUST, XLENGTH, pcgx); // redo and add eps
//     if(alarm == true){
//       printf("Estep: add_EPS_pxgy returns alarm. p(c|x) has zeros.\n");
//     }
//   }
//   if (t ==false){ 
// normalize
    for(i=0;i<XLENGTH;i++){ 
      for (k=0; k<NCLUST; k++){
	if((v = histogram_get(pcgx,k,i))){
	  histogram_set(pcgx,k,i,(v /sumup[i]));
       }
     }   
   }
  return alarm;
}

//==================================================================================================
// perturb cluster centers
void perturb_pygc(double pert, int YLENGTH, int NCLUST, double **pygc){
  int j,k;
  double sum = 0;
  double checksum = 0;
  
  printf("perturb_pygc\n");

  for(k=0;k<NCLUST;k++){
    sum = 0;
    for(j=0;j<YLENGTH;j++){
      pygc[j][k] += drand_uniform(0, 1)*pert;  // drand_uniform(0, 0.05); 
      sum += pygc[j][k];
    }
    if(sum < EPS){
      printf("perturb_pygc: ERROR. sum = %f \n", sum);
    }
    else{
      for(j=0;j<YLENGTH;j++){
	pygc[j][k] /= sum; 
      }
    }
  }
  
  for(k=0;k<NCLUST;k++){
    checksum = 0;
    for(j=0;j<YLENGTH;j++){
      checksum += pygc[j][k];
    }
    if(checksum < 1 - EPS || checksum > 1 + EPS){
      printf("perturb_pygc: NORMALIZATION ERROR: checksum = %f \n", sum);
    }
  }
}

//=================================================================================================
// split clusters
bool split(int NCLUSTnow[], int NCLUSTmax, int YLENGTH, double **pygc)
{
  int j, k;
  double checksum;
  bool Full = false;
  bool debug = true;
  int NCLUSTold = NCLUSTnow[0];
  int NCLUSTnew = NCLUSTnow[0] * 2;
  double pert = 0.005;

  printf("SPLIT:\n");
  if(debug) printf("split: NCLUSTnew = %d \t copy %d clusters: \n",NCLUSTnew, NCLUSTnow[0]);
  // make two copies centers: copy first half of pygc onto second half.
  if(NCLUSTnew <=  NCLUSTmax){
    for(k=0; k<NCLUSTold; k++){
      for(j=0;j<YLENGTH;j++){
	pygc[j][k] = pygc[j][k];
	pygc[j][k+NCLUSTold] = pygc[j][k];
      }
    }
    for(k=NCLUSTnew; k<NCLUSTmax; k++){
      for(j=0;j<YLENGTH;j++){
	pygc[j][k] = 1 / ((double) YLENGTH ); 
      }
    }
   // if(debug) print_mat_inv(YLENGTH, NCLUSTmax, pygc);
    
    perturb_pygc(pert, YLENGTH, NCLUSTmax, pygc);
    if(debug){
      printf("split: calling  perturb_pygc:\n");
      //print_mat_inv(YLENGTH, NCLUSTmax, pygc);
    }
    // check normaliztion
    for(k=0;k<NCLUSTnew;k++){
      checksum = 0;
      for(j=0;j<YLENGTH;j++){
	checksum += pygc[j][k];
      }
      if(checksum > 1+EPS || checksum < 1-EPS){
	printf("split: ERROR. normalization wrong at cluster %d. sum = %f\n", k, checksum);
      }
    }
    
    NCLUSTnow[0] = NCLUSTnew;
  } // end if
  else{
    Full = true;
  }
  if(debug){
    printf("split: Number of NEW custers: %d center distributions:\n",NCLUSTnew);
    //print_mat_inv(YLENGTH, NCLUSTmax, pygc);
  }
 
  return Full;
}

//==================================================================================================
// compute self-consistent IB equations 
// NEEDS NORMALIZED HISTOGRAM AS INPUT:  p(x,y) 
// computes p(y|x) from that input.
// Input: XLENGTH, YLENGTH,  NCLUST, p(x,y), {initial p(c|x) or p(y|c)} -- initialization
// Output: p(c|x), p(y|c), p(y|x)
// give arrays of size 2*NCLUST as input/output
	     
double EM(int XLENGTH, int YLENGTH, int NCLUST, int NCLUSTmax, double beta,
	  histogram_t *pxy,
	  histogram_t *histo,
	  histogram_t *FULLpcgx,
	  histogram_t *FULLpygc, 
	  histogram_t *FULLpc,
	  double infocurve[4],
	  bool printme)
{
  int i, j, k;
  int count;
  bool alarm = false;
  bool debug = true; //for verbose debug printing
  bool useDKLReScale = false;
  double v,conv;
  double sum = 0;

  double px[XLENGTH];
  double py[YLENGTH];
  
  double F = 0;
  double Icy = 0;
  double Icx = 0;
  double Hx = 0;
  double Ixy = 0;
  
  // copy input to smaller arrays of size of current number of clusters
  histogram_t *pcgx;
  pcgx = histogram_create(NCLUST*XLENGTH);
  printf("EM: made pcgx\n");
  histogram_t *pygc;
  pygc = histogram_create(YLENGTH*NCLUST);
  printf("EM: made pygc\n");
  histogram_t *pc;
  histogram_t *pc_old;
  pc = histogram_create(NCLUST);
  pc_old = histogram_create(NCLUST);
  
  histogram_t *DKL;
  DKL = histogram_create(XLENGTH*NCLUST);
  printf("EM: made DKL\n");
  
  for(k=0; k<NCLUST; k++){
    //init pc to flat distr
    histogram_set(pc,k,0,(1/(double)NCLUST));
    for(i=0; i<XLENGTH; i++){
      //put current valid pcgxs into pcgx
      if((v = histogram_get(FULLpcgx,k,i)))
      {
	histogram_set(pcgx,k,i,v);
      }
    }
    for(j=0;j<YLENGTH; j++){
      if((v = histogram_get(FULLpygc,j,k)))
      {
	histogram_set(pygc,j,k,v);
      }
    }
  }
  if(debug)
  {
    printf("p(c):\n");
    histogram_print(pc);
    printf("p(c|x):\n");
    histogram_print(pcgx);
    printf("p(y|c):\n");
    histogram_print(pygc);
  }

  printf("EM: NCLUST = %d\n", NCLUST);

  // get p(y) from p(x,y) assuming that p(x,y) is normalized properly
  marg_py(XLENGTH, YLENGTH, pxy, py);
  
  alarm = checknorm_vec(YLENGTH, py);
  if(alarm == true){
    printf("EM: after marg_py: checknorm_vec returns alarm.\n   p(y) is not normalized properly.\n");
    alarm = false;
  }
  printf("\nEM: Made p(y)\n");
 // print_vec(YLENGTH, py);
  
  // get p(x) from p(x,y)

  if(debug){
    printf("\nEM input P(x,y) = \n");
    //print_mat(XLENGTH, YLENGTH, pxy);
  }  
  
  
  marg_px(XLENGTH, YLENGTH, pxy, px);
  alarm = checknorm_vec(XLENGTH, px);
  if(alarm == true){
    printf("EM: after marg_px: checknorm_vec returns alarm.\n   p(x) is not normalized properly.\n");
  }
  printf("\n EM: made p(x)\n");
  // print_vec(XLENGTH, px);
  
  // compute conditional probability p(y|x) = p(x,y)/p(x)
//  conditional(XLENGTH, YLENGTH, pxy, px, pygx);
  
  if(printme == true){
    printf("\n EM: input to while: p(y|x)\n");
   // print_hist_inv(XLENGTH, YLENGTH, histo);
  }
  conv = 1;
  count = 0;
 
  if(debug){
    printf("EM: initial cluster centers: \n");
   // print_hist_inv(YLENGTH, NCLUST, pygc);
  }

  // core loop.
 //********* while(conv >= 0.1/beta){ // 0.000001){ // scale by beta.
    
    // compute DKL. input: p(y|x) p(y|c)
  if(useDKLReScale){
    DKL_rescaled(XLENGTH, YLENGTH, NCLUST, histo, pygc, DKL);  
  }else{
    DKL_Prod(XLENGTH, YLENGTH, NCLUST, histo, pygc, DKL, beta);  
  }
     if(printme == true){
      printf("\n EM: DKL at count = %d\n", count);
    //  print_mat(XLENGTH, NCLUST, DKL);
    
      printf("\n EM: p(c) at count = %d\n", count);
    //  print_vec(NCLUST, pc);
    } 
    // compute p(c|x) with Estep. input: DKL, p(c)
    alarm = Estep(XLENGTH, NCLUST, beta, DKL, pc, pcgx);
    if(alarm == true){
      printf("EM: Estep returns alarm. At iteration %d\n", count);
      alarm = false;
    }
/*
    if(printme == true){
      printf("\n EM: p(c|x) at count = %d\n", count);
     // print_mat2(NCLUST, XLENGTH, pcgx);
    }
    
    // compute p(y|c). input: p(y|x), p(c|x), p(x).
    alarm = Mstep(XLENGTH, YLENGTH, NCLUST, histo, px, pcgx, pygc);
    if(alarm == true){
      printf("EM: Mstep returns alarm. At iteration %d\n", count);
      alarm = false;
    }

    if(printme == true){
      printf("\n EM: p(y|c) at count = %d\n", count);
     // print_mat2(YLENGTH, NCLUST, pygc);
    }

    //save p(c) in p(c_old)
    copy_vec(NCLUST, pc, pc_old);
    
    //compute new p(c). input: p(c|x), p(x).
    alarm = P_c(XLENGTH, NCLUST, pcgx, px, pc);
    if(alarm == true){
      printf("EM: computing p(c). iteration %d. P_c returns alarm. p(c) either not normalized or has zeros.\n", count);
      alarm = false;
    }

    // convergence criterium: sum_c DKL[p(c)||p_old(c)]
    conv = D_KL(NCLUST, pc, pc_old);

    count=count+1;
    printf("EM: iteration %d; conv = %f; beta = %f\n\n", count, conv,beta);
  } //conv loop
  
  // copy back 
  for(k=0; k<NCLUST; k++){
    FULLpc[k] = pc[k];
    for(i=0; i<XLENGTH; i++){
      FULLpcgx[k][i] = pcgx[k][i];
    }
    for(j=0;j<YLENGTH; j++){
      FULLpygc[j][k] = pygc[j][k]; 
    }
  }
  for(k=NCLUST; k<NCLUSTmax; k++){ //fill rest with zeros
    FULLpc[k] = EPS;
    for(i=0; i<XLENGTH; i++){
      FULLpcgx[k][i] = EPS;
    }
    for(j=0;j<YLENGTH; j++){
      FULLpygc[j][k] = EPS; 
    }
  }
      
  // compute H(X)
  Hx = H_X(XLENGTH, px);
  infocurve[0] = Hx;

  // compute I(X;Y)
  Ixy = I_XY(XLENGTH, YLENGTH, pxy, px, py);
  infocurve[1] = Ixy;

  // compute I(C;X)
  Icx = I_YgX(XLENGTH, NCLUST, pcgx, px, pc);
  infocurve[2] = Icx;
 
  // compute I(C;Y)
  Icy = I_YgX(NCLUST, YLENGTH, pygc, pc, py);
  infocurve[3] = Icy;
  
  //F = Icy - Icx/beta; // F unscaled
  F = Icy / Ixy - Icx / (Hx * beta); // F scaled
  //  printf("EM result: F = %f \n", F);
*/
//clear hashes
  histogram_destroy(DKL);
  histogram_destroy(pygc);
  histogram_destroy(pcgx);
  histogram_destroy(pc); 
  histogram_destroy(pc_old);   
  return F;
}

//=================================================================================================
// merge clusters
int merge(double beta, int NCLUSTnow, int NCLUSTmax, int YLENGTH, double **pygc)
{
  int i,j, k, l, n;
  int NCLUSTnew = 0;
  double tol = 0.0001/beta; //0.01; // have to scale that by the temperature 
  // maybe the tolerance should be a parameter we hand to the function.
  bool different = true;
 // double pygc_tmp[YLENGTH][NCLUSTmax]; 
  /*************** allocate memory for large cluster probability array **************/
  double** pygc_tmp;
  //allocate pointer memory for first dimension
  pygc_tmp = (double**)malloc(YLENGTH*sizeof(double));
  if(NULL == pygc_tmp){free(pygc_tmp); printf("Memory allocation failed while allocating for pygc[].\n"); exit(-1);}
  
  /*allocate memory for second dimension */
  for(i = 0; i < YLENGTH;i++)
  {
    pygc_tmp[i] = (double *) malloc( NCLUSTmax * sizeof(double) );
    if(NULL == pygc_tmp[i]){free(pygc_tmp[i]); printf("Memory allocation failed while allocating for matrix[x][].\n"); exit(-1);}
  }
  
  //cut off tol
  if(tol < 0.000001)
    tol = 0.000001;
  
  // copy first cluster center
  for(j=0;j<YLENGTH;j++){
    pygc_tmp[j][NCLUSTnew] =  pygc[j][0];
  } 
  NCLUSTnew++;
  
  // compare cluster center to all previous ones, if different, copy.
  for(k=1; k<NCLUSTnow;k++){
    different = true;
    for(l=0;l<k;l++){
      n=0;
      for(j=0;j<YLENGTH;j++){
	if((pygc[j][k] < pygc[j][l] + tol) && (pygc[j][k] > pygc[j][l] - tol)){
	  n++;
	}
      }
      if (n == YLENGTH){ // center same
	different = false;
      }  
    }
    if(different == true){ // if different, copy
      for(j=0;j<YLENGTH;j++){
	pygc_tmp[j][NCLUSTnew] =  pygc[j][k];
      } 
      NCLUSTnew++;
    }
  }
  
  for(k=0; k<NCLUSTnew; k++){
    for(j=0;j<YLENGTH;j++){
      pygc[j][k] = pygc_tmp[j][k];
    }
  } 
  // set unused cluster centers to zero 
  for(k=NCLUSTnew; k<NCLUSTmax; k++){
    for(j=0;j<YLENGTH;j++){
      pygc[j][k] = 0; //1 / ( (double) YLENGTH );
    }
  }
  for(i = 0; i < YLENGTH; i++)
    free(pygc_tmp[i]);
  free(pygc_tmp);
  return NCLUSTnew;
}

//=================================================================================================
//check if cluster assignments are (almost) deterministic
bool det(int NCLUSTnew, int NCLUSTmax, int XLENGTH, double **pcgx){
  bool De = false;
  int n = 0;
  int i, k;
  double deviation =  0.000001;
  for(k=0; k<NCLUSTnew; k++){
    for(i=0;i<XLENGTH;i++){
      if(pcgx[k][i] > 1 - deviation || pcgx[k][i] < deviation ){
	n++;
      }
    }
  }
  if(n == NCLUSTnew * XLENGTH)
    De = true;
  return De;
}

//=================================================================================================
//check if cluster assignments are (almost) deterministic
bool used(double usemin, int NCLUST, double pc[]){
  bool use = true;
  int k;
  
  for(k=0;k<NCLUST;k++){
    if(pc[k] < usemin){
      use = false;
    }
  }
  return use;
}

//=================================================================================================

//=================================================================================================
// annealing loop. uses function EM() as core loop. the last value indicates if output is printed while computing.

double anneal(int XLENGTH, int YLENGTH, int NCLUST, double beta_start,
	       histogram_t *pxy,
	       histogram_t *histo,
	       histogram_t *pcgx,
	       histogram_t *pygc, 
	       histogram_t *pc,
	       double infocurve[4],
	       double annealingrate, bool plotme)
{
  printf("ANNEAL:\n");
  // starting inverse temperature: b = 1/T
  double b = beta_start;
  // value of objective function
  double F = 0;
  // variables to indicate: if true then
  bool Full = false; // clusters have split up to a pre-set number of K_max = NCLUST
  bool De = false; // solution is deterministic
  bool Use = false; // all clusters are being used in solution (no cluster is dropped)
  bool debug = false; //verbose debug output
  // counting the number of clusters
  int NCLUSTnow[1];
  int NCLUSTnew = 1;
  // counters. use i for x-dim j for y-dim and k for clusters
  int i, j, k;
  
  if(debug)
  {
    printf("made all ints\n");
    printf("NCLUST = %d\n",NCLUST);
    printf("XLENGTH= %d\n",XLENGTH);
  }
  
  // hashes to host solution, as it is splitting and merging (have to be of size 2*K_max)

  histogram_t *FULLpcgx;  
  histogram_t *FULLpygc;
  histogram_t *FULLpc;
  
  FULLpcgx = histogram_create(2*NCLUST*XLENGTH);
  FULLpygc = histogram_create(YLENGTH*2*NCLUST);
  FULLpc = histogram_create(2*NCLUST);
  
  printf("made FULL matrices\n");
  // the variable "usemmin" serves to determine if all clusters are used. 
  double usemin = 0.00001;

  const char* assignName = "assignments.txt";
  const char* centersName = "centers.txt";
  FILE * aFile;
  FILE * cFile;

  if(plotme == true){ 
    printf("START ANNEALING\n==================\n");
    if(debug)
    {
      printf("Input distribution (histograms):\n");
      print_mat(XLENGTH,YLENGTH,histo);
    }
  }
  // start with 1 cluster. initialize FULLpygc. only the first column is meaningful, but i initialize the whole thing to something random
  if(plotme) {printf("init FULLpygc\n");}
  ini_pygc(2*NCLUST,YLENGTH,2*NCLUST,FULLpygc); 
  
  // run the OCI/IB/EM algorithm. It uses only one cluster here! NCLUSTnew = 1.
  F = EM(XLENGTH, YLENGTH, NCLUSTnew, 2*NCLUST, b, pxy, histo, FULLpcgx, FULLpygc, FULLpc, infocurve, plotme);
  /*
  if(plotme == true){  
    printf("Done with first IB step in annealing.\nINITIAL CENTER (unused centers in this array are set to zero):\n"); 
   // print_mat_inv(YLENGTH,2*NCLUST,FULLpygc);
    printf("==================================\n");
  }
   // let clusters split, anneal until all NCLUST clusters are being used.
  int loopc = 0;
  while(Full == false && NCLUSTnew <= NCLUST ){
    loopc += 1;
    if(NCLUSTnew == NCLUST)
      break;
    NCLUSTnow[0] = NCLUSTnew;
    Full = split(NCLUSTnow, 2*NCLUST, YLENGTH, FULLpygc);

    NCLUSTnew = NCLUSTnow[0];
    b *= annealingrate;	// change temperature
    printf("Anneal:finished splitting, calling EM\n");
    F = EM(XLENGTH, YLENGTH, NCLUSTnew, 2*NCLUST, b, pxy, histo, FULLpcgx, FULLpygc, FULLpc, infocurve, modus, false);
    if(debug){
      printf("output EM. ASSIGNMENTS before merging. Pcgx. x=cluster: \n");
      //print_mat2(2*NCLUST, XLENGTH, FULLpcgx);
    }
    if(plotme == true){   
	printf("Loop %d: b = %f\n", loopc, b);
 	printf(" CENTERS before merging:\n");
 	//print_mat_inv(YLENGTH,2*NCLUST,FULLpygc);
    }
    NCLUSTnew = merge(b,NCLUSTnew, 2*NCLUST, YLENGTH, FULLpygc);
    if(plotme == true){
      printf(" CENTERS after merging:\n");
     // print_mat_inv(YLENGTH,2*NCLUST,FULLpygc);
      printf("Number of clusters = %d. DONE with loop %d\n",NCLUSTnew, loopc);
    }
//Full = true;
  }//end first while
  
  if(plotme == true) printf("ANNEAL: done splitting and merging\n");
  // one more iteration after final merging, but using the same temperature. 
  F = EM(XLENGTH, YLENGTH, NCLUSTnew, 2*NCLUST, b, pxy, histo, FULLpcgx, FULLpygc, FULLpc, infocurve, modus, false);
  // copy solution back
  for(k=0;k<NCLUST;k++){
    pc[k] = FULLpc[k];
    for(j=0; j<YLENGTH; j++){
      pygc[j][k] = FULLpygc[j][k];
    }
    for(i=0; i<XLENGTH; i++){
      pcgx[k][i] = FULLpcgx[k][i];
    }
  }
  if(plotme == true){ 
    printf("\nDONE FILLING UP CLUSTERS!!!\nNumber of clusters: NCLUSTnew = %d\n F = %f\nCenters are:\n", NCLUSTnew, F);
   // print_mat_inv(YLENGTH,NCLUST,pygc);
  }
  De = false;
  Use = false;
  // check if solution is deterministic
  De = det(NCLUSTnew, NCLUST, XLENGTH, pcgx);
  
  if(plotme == true){
    if(De == false)
      printf("Solution is NOT deterministic.\n");
    if(De == true)
      printf("Solution is deterministic.\n");
    printf("ASSIGNMENTS (solution): \n");
   // print_mat_inv(NCLUST, XLENGTH, pcgx);
    printf("P(c) = \n");
  //  print_vec(NCLUST,pc);
    printf("=====================================================\nContinue annealing until solution is deterministic.\n=====================================================\n");
  }
  // continue annealing until either --- assignments are deterministic
  //                          or --- p(c) converges
  double oldPc[NCLUST]; 
  double conv = 100.0;
  while(De == false && conv > 0.001){
    //printf("beta = %f\n", b);
    if(b > 100.0)
      die("Anneal: Not converging to hard assignments. EXITING");
    // allow splitting and merging to avoid degeneracies
    NCLUSTnow[0] = NCLUSTnew;
    Full = split(NCLUSTnow, 2*NCLUST, YLENGTH, FULLpygc);
    NCLUSTnew = NCLUSTnow[0];
    // save old p(c)
    copy_vec(NCLUST,pc,oldPc);
    
    b *= annealingrate;	// change temperature    
    F = EM(XLENGTH, YLENGTH, NCLUSTnew, 2*NCLUST, b, pxy, histo, FULLpcgx, FULLpygc, FULLpc, infocurve, modus, false);
    
    if(plotme == true){   
      printf("Loop %d: b = %f\n", loopc, b);
      printf(" CENTERS before merging:\n");
     // print_mat_inv(YLENGTH,2*NCLUST,FULLpygc);
    }
    
    NCLUSTnew = merge(b,NCLUSTnew, 2*NCLUST, YLENGTH, FULLpygc);
    
    if(plotme == true){
      printf(" CENTERS after merging:\n");
      //print_mat_inv(YLENGTH,2*NCLUST,FULLpygc);
    }
    
    // one more iteration after final merging, but using the same temperature. 
    F = EM(XLENGTH, YLENGTH, NCLUSTnew, 2*NCLUST, b, pxy, histo, FULLpcgx, FULLpygc, FULLpc, infocurve, modus, false);
    // copy solution back
    for(k=0;k<NCLUST;k++){
      pc[k] = FULLpc[k];
      for(j=0; j<YLENGTH; j++){
	pygc[j][k] = FULLpygc[j][k];
      }
      for(i=0; i<XLENGTH; i++){
	pcgx[k][i] = FULLpcgx[k][i];
      }
    }
    
    // check if assignments are deterministic:
    De = det(NCLUSTnew, NCLUST, XLENGTH, pcgx);
    //De = det(NCLUST, NCLUST, XLENGTH, pcgx);
    if(plotme == true){
      printf("b = %f \t F = %f\n ASSIGNMENTS: \n", b, F);
     // print_mat_inv(NCLUST, XLENGTH, pcgx);
      printf("\n CENTERS:\n");
     // print_mat_inv(YLENGTH,NCLUST,pygc);
     // printf("\n");
    }
    // check if p(c) converged
    // Calculate error between P(s) and P(s')
    double min_error = 0.0; 
    double temp_error = 999999.9;
    double partial_error=0.0;
    for(i=0; i<NCLUST; i++)
    {
      temp_error = 999999.9;
      for(j=0; j<NCLUST; j++)
      {
	partial_error = fabs(pc[j] - oldPc[i]);
	if(partial_error < temp_error)
	  temp_error = partial_error;
      }
      min_error += temp_error;
    }    
    conv = min_error;
    printf("Conv = %d\n", conv);
  }//end while loop 2

// write results to files
  aFile = fopen(assignName, "w");
  cFile = fopen(centersName, "w");
  if(aFile){
 // print_mat_inv(NCLUST, XLENGTH, pcgx);
    for(j=0;j<XLENGTH;j++){
      fprintf(aFile,"%d: ",j);
      for(i=0;i<NCLUST;i++){
	fprintf(aFile,"%f ",pcgx[i][j]);
      }
      fprintf(aFile,"\n");
    }
  }  
  else{
    printf("ERROR OPENING assignments.txt\n");
  }
  if(cFile){
  //printf("\n CENTERS:\n");
    for(j=0;j<NCLUST;j++){
      fprintf(cFile,"%d: ",j);
      for(i=0;i<YLENGTH;i++){
	fprintf(cFile,"%f ",pygc[i][j]);
      }
      fprintf(cFile,"\n");
    }
  }
  else{
    printf("ERROR OPENING centers.txt\n");
  }
  fclose(cFile);
  fclose(aFile);
  
  for(i = 0; i < YLENGTH; i++)
    free(FULLpygc[i]);
  free(FULLpygc);
  for(i = 0; i < 2*NCLUST; i++)
    free(FULLpcgx[i]);
  free(FULLpcgx);
  
  */
    if(plotme == true){
      printf("DONE ANNEALING\n==================\n");
    }
    return F;
   
}

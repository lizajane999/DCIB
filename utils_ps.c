// Some useful Functions
// For 1-dim. temporal ising-like model
// finding predictive states
// Susanne Still, 2004.

#include "utils.h"
#include "utils_ps.h"

//===========================================================
void print_vec(int length, double * vec)
{
  int j;
  for(j=0;j<length;j++)
    printf("%f\n",vec[j]);
}

//===========================================================
void print_vec_INT(int length, int * vec)
{
  int j;
  for(j=0;j<length;j++)
    printf("%d\n",vec[j]);
}

//===========================================================
void print_mat(int xlength, int ylength, histogram_t *histo)
{
  int i, j;
  double p;
  for(i=0;i<xlength;i++){
    printf("x%d: ",i);
    for(j=0;j<ylength;j++){
      if((p = histogram_get(histo,i,j))){
	printf("%d:%f ",j,p);
      }
//       else{
// 	printf("%f ",z);
//       }
    }
    printf("\n");
  }
}

//===========================================================
void print_mat2(int xlength, int ylength, double mat[xlength][ylength])
{
  int i, j;
  for(i=0;i<xlength;i++){
    printf("x%d: ",i);
    for(j=0;j<ylength;j++){
      printf("%f ",mat[i][j]);
    }
    printf("\n");
  }
}
//===========================================================
void print_mat_INT(int xlength, int ylength, int mat[xlength][ylength])
{
  int i, j;
  for(i=0;i<xlength;i++){
    for(j=0;j<ylength;j++){
      printf("%d ",mat[i][j]);
    }
    printf("\n");
  }
}

//===========================================================
void print_hist_inv(int xlength, int ylength, histogram_t *mat)
{
  int i, j;
  double p;
  for(j=0;j<ylength;j++){
    printf("y%d: ",j);
    for(i=0;i<xlength;i++){
      if((p = histogram_get(mat,i,j))){
	printf("x%d:%f ", i,p);
      }
    }
    printf("\n");
  }
}
//===========================================================
void print_mat_inv(int xlength, int ylength, double **mat)
{
  int i, j;
  for(j=0;j<ylength;j++){
    printf("y%d: ",j);
    for(i=0;i<xlength;i++){
      printf("%f ",mat[i][j]);
    }
    printf("\n");
  }
}

//===========================================================
void write_mat(int xlength, int ylength, double mat[xlength][ylength], char * filename)
{
  int i, j;
  FILE * fid;
  fid = fopen(filename,"a");
  for(j=0;j<ylength;j++){
    for(i=0;i<xlength;i++){
      fprintf(fid,"%f\t",mat[i][j]);
    }
    fprintf(fid,"\n");
  }
  fclose(fid);
}

//===========================================================
void copy_mat(int xlength, int ylength, double mat[xlength][ylength], double newmat[xlength][ylength])
{
  int i, j;
  for(j=0;j<ylength;j++){
    for(i=0;i<xlength;i++){
      newmat[i][j] = mat[i][j];
    }
  }
}

//===========================================================
void invert_mat(int xlength, int ylength, double mat[xlength][ylength], double newmat[ylength][xlength])
{
  int i, j;
  for(j=0;j<ylength;j++){
    for(i=0;i<xlength;i++){
      newmat[j][i] = mat[i][j];
    }
  }
}

//===========================================================
void copy_vec(int length, double * vec, double * newvec)
{
  int j;
  for(j=0;j<length;j++)
    newvec[j] = vec[j];
}

//===========================================================
double sum_mat(int xlength, int ylength, double mat[xlength][ylength])
{
  double sum = 0;
  int i, j;
  for(j=0;j<ylength;j++){
    for(i=0;i<xlength;i++){
      sum += mat[i][j];
    }
  }
  return sum;
}

//===========================================================
double sum_vec(int length, double * vec)
{
  double sum = 0;
  int j;
  for(j=0;j<length;j++)
    sum += vec[j];
  return sum;
}

//===========================================================
double add_EPS_pxy(int xlength, int ylength, histogram_t *pxy)
{
  double p, checksum = 0;
  double sum = 0;
  int i, j;
  for(j=0;j<ylength;j++){
    for(i=0;i<xlength;i++){
      if((p = histogram_get(pxy,i,j))){
	histogram_add(pxy, i, j, EPS);
	sum += (p + EPS);
      }
    }
  }
  for(j=0;j<ylength;j++){
    for(i=0;i<xlength;i++){
      if((p = histogram_get(pxy,i,j))){
	histogram_set(pxy, i, j, (p/sum));
	checksum += p/sum;
      }
    }
  }
  if(checksum > 1+EPS || checksum < 1-EPS){
    printf("add_EPS_pxy: ERROR renormalizing: sum = %f\n", checksum);
  }
  return checksum;
}

//===========================================================
bool add_EPS_pxgy(int xlength, int ylength, double **pxgy)
{
  double checksum[ylength];
  double sum[ylength];
  int i, j;
  bool t = false;

  for(j=0;j<ylength;j++){
    sum[j] = 0;
    checksum[j] = 0;
    for(i=0;i<xlength;i++){
      pxgy[i][j] += EPS;
      sum[j] += pxgy[i][j];
    }
  }
  for(j=0;j<ylength;j++){
    for(i=0;i<xlength;i++){
      pxgy[i][j] /= sum[j];
      checksum[j] += pxgy[i][j];
    }
    if(((double)checksum[j] > ((double)1+EPS)) || ((double)checksum[j] < ((double)1-EPS))){
      printf("add_EPS_pxy: ERROR renormalizing: sum[%d] = %f\n", j, checksum[j]);
      t = true;
    }
  }
  return t;
}

//=========================================================
bool checknorm_vec(int xlength, double px[xlength])
{
  double checksum = 0.0;
  int i;
  bool alarm = false;
  
  for(i=0;i<xlength;i++){
    checksum += px[i];
  }
  if((double)checksum > (double)1 + (double)EPS || (double)checksum < (double)1 - (double)EPS){
    alarm = true;
    printf("checknorm_vec: not normalized to 1. sum = %f\n", checksum);
  }
  return alarm;
}

//=======================================================



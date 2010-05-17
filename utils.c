// Some useful Functions
// Susanne Still, 2004.

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
// #include "fpu_control.h"
#include "utils.h"

double 
drand_uniform(double low, double high)
{
  return low + (double)random() * (high - low) / RAND_MAX;
}

int 
irand_uniform(int low, int high)
{
  return (int) floor(drand_uniform(low, high+1));
}

double 
drand_gauss(double mean, double sdev)
{
  double r, theta;
  // Make a positive number distributed as  r exp(-r*r/2)
  // by applying the inverse of the cumulative 1.0 - exp(-r*r/2) 
  // to a uniform number in range [0,1).
  r = sqrt( -2.0 * log (1.0 - drand_uniform(0,1)));
  // Make a number between -PI and PI
  theta = drand_uniform(-M_PI, M_PI);
  // Take the x coordinate and rescale
  return mean + sdev * r * sin(theta);
}

int 
randspin_p(double p)
{
  double r;
  int spin;
  r = (double)random() / RAND_MAX;
  // printf("randspin_p: r = %f\n",r);
  if (r < p)
    spin = 1;
  else 
    spin = -1;
  return spin;
}


int 
max_int(int x, int y)
{
  int k = x;
  if(y > x) 
    k = y;
  return k;
}

int 
min_int(int x, int y)
{
  int k = x;
  if(y < x) 
    k = y;
  return k;
}

void 
enable_sigfpe(bool enable)
{
//  fpu_control_t cw;
//  _FPU_GETCW(cw);
//  cw |= (_FPU_MASK_IM|_FPU_MASK_ZM|_FPU_MASK_OM);
//  if (enable)
//    cw &= ~(_FPU_MASK_IM|_FPU_MASK_ZM|_FPU_MASK_OM);
//  _FPU_SETCW(cw);
}

void die(char *s)
{
  fprintf(stderr,"%s\n",s);
  exit(10);
}

void integer_to_bitvector(int x, int length, int vec[]) // first entry in vec == least significant bit
{
  int i;
  for (i=0; i<length; i++)
    {
      vec[i] = x % 2;
      x = x / 2;
    }
}

int bitvector_to_integer(int length, int vec[])
{
  int x;
  int i;
  double xd = 0;
  for (i=0; i<length; i++)
    {
      xd += vec[i] * pow(2, (double) i);
    }
  x = (int) xd;
  return x; 	
}

void sort_label(int n, double ar[], int label[]) 
{
  int i;
  double tmp;
  int tmplabel;

  for(i=0; i<n; i++){
    label[i] = i; 
  }

  i = 0;
  while (i < n) {
    if (i == 0 || ar[i-1] <= ar[i]) 
      i++;
    else 
      {
	tmp = ar[i];   tmplabel = label[i];
	ar[i] = ar[i-1];   label[i] = label[i-1]; label[i-1] = tmplabel; 
	ar[--i] = tmp; 
      }
  }
}

// histogram returns totoal # of data.
int hist_double(const char * infile, int lbins, double edges[], double h[])
{ 
  FILE * fin = fopen(infile,"r");
  double x;
  int j, k; 
  int N = 0; 
  // ASSERT(fin);
  
  for( j = 0; j < lbins; j++)
    h[j] = 0;
  
  while( fread( &x, sizeof(double), 1, fin) != EOF){
    N += 1;
    k = 0;
    if(x < edges[0])
      k = 0;
    else{
      for( j = 0; j < lbins-1; j++){
	if(x > edges[j] && x < edges[j+1])
	  k = j+1;
  
      }
      if(x > edges[lbins-1]){
	printf("hist_double: ERROR: upper edge of bins must be max of vector to bin.\n");
      }
    }
    h[k] = h[k]+1;
  }
  fclose(fin);
  return N;
}

float distance(float a, float b){
  
  float dist = (a - b);	
  dist *=dist;
  dist = sqrt(dist);
  return dist;

}

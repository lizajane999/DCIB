// Document Classification using IB with annealing
// Lisa Miller
// Susanne Still
// University of Hawaii at Manoa, 2010.
// Thanks to Ben Karsin and to Leon Bottou for contributing routines


#include "utils.h"
#include "utils_IB.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

// Names of files containing preprocessor data
const char* histoName = "data/histograms.txt";
const char* probName = "data/probs.txt";


//=================================
void myusage(void)
{
  fprintf(stderr,"Usage: dcib\n");
  fprintf(stderr,"	XLENGTH\n");
  fprintf(stderr,"	YLENGTH\n");
  fprintf(stderr,"	NCLUST\n");
  exit(10);
}
//===========================================================


int main(int argc, char * argv[])
{
  //============================================
  // Set input parameters
   if (argc < 4|| argc > 4){
     printf("Wrong number of inputs. You entered %d. You should enter 4\n.", argc);
     myusage();
   }
  int XLENGTH = atoi(argv[1]);
  int YLENGTH = atoi(argv[2]);
  int NCLUST = atoi(argv[3]);
  int HLENGTH;
  double px[XLENGTH];
  double py[YLENGTH];
  double beta = 30.0;
  double F;
  double infocurve[4] = {1,1,1,1};
  //int modus = 0;
  // modus = 0 : use p(y|c) cluster centers as initialization 
  // modus = 1 : use p(c|x) assignment ruls as initialization 
  double annealingrate = 1.1;
  bool plotme = true; //some debug printing
  bool debug = false;//verbose debug printing
  

  FILE * histFile;
  FILE * probFile;
  histogram_t *histo;
  histogram_t *pxy;
  histogram_t *pygc;
  histogram_t *pcgx;
  histogram_t *pc;
  
  if(plotme)
  {
    printf("XLENGTH= %d\n", XLENGTH);
    printf("YLENGTH= %d\n", YLENGTH);
    printf("NCLUST = %d\n", NCLUST);
  }
  
  histFile = fopen(histoName, "r");
  if(histFile == NULL)
  {
    printf("Error: can't access %s\n", histoName);
    exit(1);
  }
  float v;
  int j,i = 0;
  int y;
  int x;
  char s [YLENGTH*20];
  char * tok;
  
//   //allocate memory for large cluster probability arrays
//   double** pygc;
//   //allocate pointer memory for first dimension
//   pygc = (double**)malloc(YLENGTH*sizeof(double));
//   if(NULL == pygc){free(pygc); printf("Memory allocation failed while allocating for pygc[].\n"); exit(-1);}
//   
//   /*allocate memory for second dimension */
//   for(i = 0; i < YLENGTH;i++)
//   {
//     pygc[i] = (double *) malloc( NCLUST * sizeof(double) );
//     if(NULL == pygc[i]){free(pygc[i]); printf("Memory allocation failed while allocating for matrix[x][].\n"); exit(-1);}
//   }
// 
//   double** pcgx;
//   //allocate pointer memory for first dimension
//   pcgx = (double**)malloc(NCLUST*sizeof(double));
//   if(NULL == pcgx){free(pcgx); printf("Memory allocation failed while allocating for pcgx[].\n"); exit(-1);}
//   
//   /*allocate memory for second dimension */
//   for(i = 0; i < NCLUST;i++)
//   {
//     pcgx[i] = (double *) malloc( XLENGTH * sizeof(double) );
//     if(NULL == pcgx[i]){free(pcgx[i]); printf("Memory allocation failed while allocating for matrix[x][].\n"); exit(-1);}
//   }
  
  //count how many entries are actually in the histogram file -- these are only the non-zero p(y|x)
  while(fscanf(histFile, "%f", &v)!= EOF)
  {
     i++;
  }
   rewind(histFile);
   if(debug){printf("Counted %d entries in histogram file\n", i);}
   //make histogram size 1/2 number entries counted because 1/2 are the indices (y values)
   HLENGTH = i/2;
   histo = histogram_create(HLENGTH);
   x = 0;
  //make p(y|x) histogram
  //YLENGTH*20 is assuming most words are less than 20 chars long
  while(fgets(s, (YLENGTH*20), histFile)!= NULL)
  {
    tok = strtok(s,"\t");
    while(tok != NULL)
    {
      if(debug){printf("Found word: %d\n", atoi(tok));}
      y = atoi(tok);
      tok = strtok(NULL, "\t");
      if(tok != NULL) //hits null in the word 
      {
	if(debug){printf("Found value: %f\n", atof(tok));}
	v = atof(tok);
	tok = strtok(NULL, "\t");
	histogram_set(histo,x,y,v);
      }       
    }
    x++;
  }
  if(debug)
  {
    printf("Checking histogram\n");
    for(x = 0; x < XLENGTH; x++)
    {
      for(y = 0; y < YLENGTH; y++)
      {
	  if((v = histogram_get(histo, x, y)))
	  {
	    printf("Found doc %d - word %d in histo: %f\n",x,y,v); 
	  }
      }
    }
  }
  fclose(histFile);
  // now read in the other initial probs
  probFile = fopen(probName, "r");
  if(probFile == NULL)
  {
    printf("Error: can't access %s\n", histoName);
    exit(1);
  }
//*********** p(x) is first line *******************
  printf("read in p(x)\n");
  fgets(s, (XLENGTH*20), probFile);
  //printf("s = %s\n",s);
 // printf("length of s = %d\n", strlen(s));
  tok = strtok(s,"\t");
  x = 0;
  while(tok != NULL)
  {
    //printf("tok = %s\n", tok);
    if(strcmp(tok,"\n") != 0)
    {
      px[x] = atof(tok);
      //printf("px(%d) = %f\n", x, px[x]);
      x++;
    }
    tok = strtok(NULL,"\t");
  }
  //second line is blank
  fgets(s, (YLENGTH*20), probFile);
  //************** p(y) is third line ************************
  fgets(s, (YLENGTH*20), probFile);
  printf("read in p(y)\n");
 // printf("s = %s\n",s);
 // printf("length of s = %d\n", strlen(s));
  tok = strtok(s,"\t");
  x = 0;
  while(tok != NULL)
  {
    if(strcmp(tok,"\n") != 0)
    {
      py[x] = atof(tok);
     // printf("py(%d) = %f\n", x, py[x]);
      x++;
    }
    tok = strtok(NULL,"\t");
  }
  //fourth line is blank
  fgets(s,(YLENGTH*20), probFile);
  //******************** fifth line on is p(xy)
  printf("Making p(x;y)\n");
  pxy = histogram_create(HLENGTH);
  x = 0;
  while(fgets(s, (YLENGTH*20), probFile)!= NULL)
  {
    //get line for doc
    //printf("Found %s\n", s);
    //split y and p
     tok = strtok(s,"\t");
     while(tok != NULL)
     {
      //printf("Found word: %d\n", atoi(tok));
      y = atoi(tok);
      tok = strtok(NULL, "\t");
      if(tok != NULL) //hits null in the word 
      {
	//printf("Found value: %f\n", atof(tok));
	v = atof(tok);
	tok = strtok(NULL, "\t");
	histogram_set(pxy,x,y,v);
      }       
    }
    x++;
  }
  fclose(probFile);
  if(debug){
   printf("Checking pxy\n");
    for(x = 0; x < XLENGTH; x++)
    {
      for(y = 0; y < YLENGTH; y++)
      {
	if((v = histogram_get(pxy, x, y)))
	{
	  printf("Found doc %d - word %d in histo: %f\n",x,y,v); 
	}
      }
    }
  }
//****************** initialize cluster probs *****************************
printf("Making p(c|x)\n");
ini_pcgx(NCLUST, XLENGTH, pcgx);
printf("Now going to anneal\n");
//******************* NOW GO TO ANNEAL **************************
F = anneal(XLENGTH, YLENGTH, NCLUST, beta, pxy, histo, pcgx, pygc, pc,
		   infocurve, annealingrate, plotme);
  
printf("Returned from anneal\n");  
  //cleanup
  histogram_destroy(pxy);
  histogram_destroy(histo);
  histogram_destroy(pygc);
  histogram_destroy(pcgx);
  histogram_destroy(pc);  
//   for(i = 0; i < YLENGTH; i++)
//     free(pygc[i]);
//   free(pygc);
//   for(i = 0; i < NCLUST; i++)
//     free(pcgx[i]);
//   free(pcgx);
  
  exit(1);
}

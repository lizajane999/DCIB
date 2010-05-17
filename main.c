// Document Classification using IB with annealing
// Lisa Miller
// Susanne Still
// University of Hawaii at Manoa, 2010.
// Thanks to Ben Karsin and to Leon Bottou for contributing routines

//#include "histo2.h"
#include "utils.h"
#include "utils_ps.h"
#include "utils_IB.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

// Names of files containing preprocessing data
const char* histoName = "data/histograms.txt";
//const char* pastName = "Online-OCI_operating/past.txt";
//const char* StateName = "Online-OCI_operating/states.txt";

//===========================================================
void myusage(void)
{
  fprintf(stderr,"Usage: NEW-OCI-FL\n");
  fprintf(stderr,"	<inputfile> (including the path)\n");
  fprintf(stderr,"	<alphabet size>\n");
  fprintf(stderr,"	<future length> \n");
  fprintf(stderr,"	<history length> \n");
  fprintf(stderr,"	<annealing rate>\n");
  fprintf(stderr,"	<---------------------- optional:> \n");
  fprintf(stderr,"	<max number of clusters to allow, K (optional)> \n");
  fprintf(stderr,"	<maximum inverse temperature> \n");
  exit(10);
}
//===========================================================


int main(int argc, char * argv[])
{
  //============================================
  // Set input parameters
  if (argc < 6 || argc > 8){
    printf("Wrong number of inputs. You entered %d. You should enter either 6 or 7 or 8\n.", argc);
    myusage();
  }
  return 1;
}

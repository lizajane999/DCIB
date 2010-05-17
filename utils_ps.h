// Some useful Functions
// For 1-dim. temporal ising-like model
// finding predictive states
// Susanne Still, 2004.

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "histo2.h"

void print_vec(int length, double * vec);
//void print_mat(int xlength, int ylength, double mat[xlength][ylength]);
void print_mat(int xlength, int ylength, histogram_t *histo);
void print_mat2(int xlength, int ylength, double mat[xlength][ylength]);
void print_mat_INT(int xlength, int ylength, int mat[xlength][ylength]);
void print_hist_inv(int xlength, int ylength, histogram_t *mat);
void print_mat_inv(int xlength, int ylength, double **mat);
void write_mat(int xlength, int ylength, double mat[xlength][ylength], char * filename);
void copy_mat(int xlength, int ylength, double mat[xlength][ylength], double newmat[xlength][ylength]);
void invert_mat(int xlength, int ylength, double mat[xlength][ylength], double newmat[ylength][xlength]);
void copy_vec(int length, double * vec, double * newvec);
double sum_mat(int xlength, int ylength, double mat[xlength][ylength]);
double sum_vec(int length, double * vec);
double add_EPS_pxy(int xlength, int ylength, histogram_t *pxy);
bool add_EPS_pxgy(int xlength, int ylength, double **pxgy);
bool checknorm_vec(int xlength, double px[xlength]);

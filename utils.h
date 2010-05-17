// Some useful Functions
// Susanne Still, 2004.

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define EPS 0.000001


//typedef enum { false, true } bool;

// few functions for generating random data

// --- uniform integer in [low,high]
int irand_uniform(int low, int high);

// --- uniform double in [low,high)
double drand_uniform(double low, double high);

// --- gaussian number with specified mean and sdev
double drand_gauss(double mean, double sdev);

// --- spin (+1 or -1) == +1 with probability p
int randspin_p(double p);

// --- maximum of two integers
int max_int(int x, int y);

// --- maximum of two integers
int min_int(int x, int y);

// make sure to get floating point exception SIGFPE instead of Nan
void enable_sigfpe(bool enable);

// terminate and print a message
void die(char *s);

// --- turn integer into bit vector
void integer_to_bitvector(int x, int length, int vec[]);

// --- turn bit vector into integer
int bitvector_to_integer(int length, int vec[]);

// macro to check an error condition

#define ASSERT(condition) \
   do { if (!(condition)) \
     die("Assertion failed"); \
   } while (0)

// simple insertion sort. modified from gnome search
void sort_label(int len, double ar[], int label[]);

// histogram of doubles in a file using lbins bins specified by the upper edges. return histogram in h.
int hist_double(const char * infile, int lbins, double edges[], double h[]);

// euclid distance between 2 numbers
float distance(float a, float b);

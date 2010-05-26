#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "histo2.h"
#include "utils.h"

typedef struct record_s record_t;

struct record_s {
  record_t *next;
  int x,y;
  double w;
};

struct histogram_s {
  int hsize;
  record_t **buckets;
};

static inline int
hash_function(int x, int y)
{
  return ((x<<10) ^ y) & 0x7fffffff;
}

histogram_t *histogram_create(int hsize)
{
  int i;
  histogram_t *self = malloc(sizeof(histogram_t));
  if(NULL == self){free(self); printf("Memory allocation failed while allocating for histogram.\n"); exit(-1);}
  if(hsize<1 || hsize>1000000000){free(self); printf("Histogram too small or too big, failed to create.\n"); exit(-1);}
  self->hsize = hsize;
  self->buckets = malloc(hsize * sizeof(record_t*));
  if(NULL == self->buckets){free(self); printf("Memory allocation failed while allocating for histogram buckets.\n"); exit(-1);}
  for (i=0; i<hsize; i++)
    self->buckets[i] = 0;
  return self;
}

void histogram_clear(histogram_t *self)
{
  int i;
  for (i=0; i<self->hsize; i++) 
    {
      record_t *d = self->buckets[i];
      while (d) 
        {
          record_t *next = d->next;
          free(d);
          d = next;
        }
      self->buckets[i] = 0;
    }
}

void histogram_destroy(histogram_t *self)
{
  if(self != NULL){
    histogram_clear(self);
    free(self->buckets);
    free(self);
  }
}

static record_t *lookup(histogram_t *self,
			int x, int y, int create)
{
  record_t *d;
  int h = hash_function(x,y) % self->hsize;
  for (d = self->buckets[h]; d; d=d->next)
    if (d->x == x && d->y == y)
      return d;
  if (! create)
    return 0;
  d = malloc(sizeof(record_t));
  if(NULL == d){free(d); printf("Memory allocation failed while allocating for record in histogram lookup.\n"); exit(-1);}
  d->x = x;
  d->y = y;
  d->w = 0.0;
  d->next = self->buckets[h];
  self->buckets[h] = d;
  return d;
}

void histogram_add(histogram_t *self, int x, int y, double w)
{
  record_t *d = lookup(self, x, y, 1);
  d->w += w;
}

double histogram_get(histogram_t *self, int x, int y)
{
  record_t *d = lookup(self, x, y, 0);
  if (d)  
    return d->w;
  return 0.0;
}

void histogram_set(histogram_t *self, int x, int y, double w)
{
  record_t *d = lookup(self, x, y, 1);
  d->w = w;
}

void histogram_divide(histogram_t *self, int x, int y, double w)
{
  record_t *d = lookup(self, x, y, 1);
  d->w /= w;
}

void histogram_mult(histogram_t *self, int x, int y, double w)
{
  record_t *d = lookup(self, x, y, 1);
  d->w *= w;
}

double histogram_entropy(histogram_t *self)
{
  int h;
  record_t *p;
  double total = 0;
  double clogc = 0.0;
  /* loop on the buckets */
  for (h=0; h<self->hsize; h++) {
    /* loop on records in list */
    for (p = self->buckets[h]; p; p=p->next) {
    //  ASSERT(p->w >= 0);
      total += p->w;
      clogc += p->w * log( p->w );
    }
  }
  /* compute entropy */
  return log(total) - clogc / total;
}

void histogram_print(histogram_t *self)
{
  int h;
  record_t *p;
  /* loop on the buckets */
  for (h=0; h<self->hsize; h++) {
    /* loop on records in list */
    for (p = self->buckets[h]; p; p=p->next) {
      //  ASSERT(p->w >= 0);
      if(p){
	printf("x:%d,y:%d value:%f\t",p->x, p->y,p->w);
      }
    }
    printf("\n");
  }
}








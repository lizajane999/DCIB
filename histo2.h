typedef struct histogram_s histogram_t;

histogram_t *histogram_create(int hashtablesize);
void histogram_destroy(histogram_t *self);
void histogram_clear(histogram_t *self);
void histogram_add(histogram_t *self, int x, int y, double w);
double histogram_get(histogram_t *self, int x, int y);
void histogram_set(histogram_t *self, int x, int y, double w);
void histogram_divide(histogram_t *self, int x, int y, double w);
void histogram_mult(histogram_t *self, int x, int y, double w);
double histogram_entropy(histogram_t *self);
void histogram_print(histogram_t *self);


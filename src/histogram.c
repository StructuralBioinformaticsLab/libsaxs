#include "histogram.h"

histogram *create_histogram(int size, histogram_type ht)
{
  histogram *hist = (histogram *) malloc(sizeof(histogram));

  switch (ht)
  {
  case INT_HIST:
    hist->data.int_bin = (int *) malloc(size * sizeof(int));
    memset(hist->data.int_bin, 0, size * sizeof(int));
    break;
  case FLOAT_HIST:
    hist->data.float_bin = (double *) malloc(size * sizeof(double));
    memset(hist->data.float_bin, 0, size * sizeof(double));
    break;
  default:
    free(hist);
    return NULL;
  }
  hist->type = ht;
  hist->num_bins = size;

  return hist;
}

void destroy_histogram(histogram *hist)
{
  if (hist != NULL)
  {
    if (hist->data.int_bin != NULL && hist->type == INT_HIST)
    {
      free(hist->data.int_bin);
    }
    else if (hist->data.float_bin != NULL && hist->type == FLOAT_HIST)
    {
      free(hist->data.float_bin);
    }

    free(hist);
  }
}

void print_histogram(histogram *hist)
{
  fprint_histogram(stdout, hist);
}

void fprint_histogram(FILE *fp, histogram *hist)
{
  for (int i = 0; i < hist->num_bins; ++i)
  {
    if (hist->type == INT_HIST)
      fprintf(fp, "%d\t%d", i, hist->data.int_bin[i]);
    else if (hist->type == FLOAT_HIST)
      fprintf(fp, "%d\t%f", i, hist->data.float_bin[i]);
    fprintf(fp, "\n");
  }
}

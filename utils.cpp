#include "utils.h"

void swap (double *a, const int i, const int j)
{
  double b = a[i];
  a[i] = a[j];
  a[j] = b;
}

void swap_arrays (double *a, double *b, const int n)
{
  double buf[8];
  int i;
  int k = n & 7;
  for (i = 0; i < n - 8; i+= 8)
    {
      buf[0] = a[i]; 
      buf[1] = a[i + 1]; 
      buf[2] = a[i + 2]; 
      buf[3] = a[i + 3]; 
      buf[4] = a[i + 4]; 
      buf[5] = a[i + 5]; 
      buf[6] = a[i + 6]; 
      buf[7] = a[i + 7]; 
      a[i] = b[i];
      a[i + 1] = b[i + 1];
      a[i + 2] = b[i + 2];
      a[i + 3] = b[i + 3];
      a[i + 4] = b[i + 4];
      a[i + 5] = b[i + 5];
      a[i + 6] = b[i + 6];
      a[i + 7] = b[i + 7];
      b[i] = buf[0];
      b[i + 1] = buf[1];
      b[i + 2] = buf[2];
      b[i + 3] = buf[3];
      b[i + 4] = buf[4];
      b[i + 5] = buf[5];
      b[i + 6] = buf[6];
      b[i + 7] = buf[7];
    }
  for (i = n - k; i < n; i++)
    {
      buf[0] = a[i];
      a[i] = b[i];
      b[i] = buf[0];
    }

}

void log_error (const int rank, const char *msg, FILE* fout)
{
  fprintf (fout, "[PROC %d:ERROR]: %s\n", rank, msg);
}

void log_info (const int rank, const char *msg, FILE* fout)
{
  fprintf (fout, "[PROC %d:INFO]: %s\n", rank, msg);
}

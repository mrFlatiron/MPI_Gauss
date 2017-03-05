#include <math.h>
#include "init_functions.h"

double gilbert (const int i, const int j, const int n)
{
  (void)n;
  return 1 / (double)(1 + i + j);
}

double abs_i_minus_j (const int i, const int j, const int n)
{
  (void)n;
  return fabs (i - j);
}

double consecutive (const int i, const int j, const int n)
{
  return i * n + j;
}

double all_null (const int i, const int j, const int n)
{
  (void)i;
  (void)j;
  (void)n;

  return 0;
}

double all_const (const int i, const int j, const int n)
{
  (void)i;
  (void)j;
  (void)n;
  return CONST;
}

double I (const int i, const int j, const int n)
{
  (void)n;
  return (i == j) ? 1 : 0;
}

double m1 (const int i, const int j, const int n)
{
  double ar[16] = {5, 7, 6, 5,
                   7, 10, 8, 7,
                   6, 8, 10, 9,
                   5, 7, 9, 10};
  return ar[j + i * n];
}

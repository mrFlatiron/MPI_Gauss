#include "matrix.h"
#include "const.h"
#include "utils.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


double matrix_sq_norm (double *a, const int n)
{
  double max = 0;
  double sum[4];
  int r = n  & 3;
  int j;
  for (int i = 0; i < n; i++)
    {
      sum[0] = 0;
      sum[1] = 0;
      sum[2] = 0;
      sum[3] = 0;
      for (j = 0; j + 3 < n; j += 4)
        {
          sum[0] += fabs (a[i*n + j]);
          sum[1] += fabs (a[i*n + j + 1]);
          sum[2] += fabs (a[i*n + j + 2]);
          sum[3] += fabs (a[i*n + j + 3]);
        }
      sum[0] += sum[1] + sum[2] + sum[3];
      for (int kr = 0; kr < r; kr++)
        sum[0] += fabs (a[i*n + j  +kr]);
      if (max < sum[0])
        max = sum[0];
    }
  return max;
}

double matrix_norm (double *a, const int rows, const int cols)
{
  double max = 0;
  for (int i = 0; i < rows; i++)
    {
      double sum = 0;
      for (int j = 0; j < cols; j++)
        sum += fabs (a[i*cols + j]);
      if (max < sum)
        max = sum;
    }
  return max;

}

void matrix_sq_init_e (double *a, const int n)
{
  int r;
  for (int i = 0; i < n; i++)
    a[i*n + i] = 1;
  for (int i = 0; i < n; i++)
    {
      int j;
      for (j = i + 1; j + 1 < n; j += 2)
        {
          a[i*n + j] = 0;
          a[j*n + i] = 0;
          a[(j+1)*n + i] = 0;
          a[i*n + j + 1] = 0;
        }
      r = (n - i - 1) & 1;
      for (int kr = 0; kr < r; kr++)
        {
          a[i*n + j + kr] = 0;
          a[(j+kr)*n + i] = 0;
        }
    }
}

void matrix_subtr_rw (double *a, double *b, const int rows, const int cols)
{
  int i, r = rows * cols & 3;
  for (i = 0; i + 3 < rows * cols; i += 4)
    {

      a[i] -= b[i];
      a[i + 1] -= b[i  +1];
      a[i + 2] -= b[i + 2];
      a[i + 3] -= b[i + 3];
    }
  for (int kr = 0; kr < r; kr++)
    a[i + kr] -= b[i + kr];
}

int matrix_reverse (double *a, double *rev, int *sp ,const int n)
{
  matrix_sq_init_e (rev, n);
  for (int i = 0; i < n; i++)
    {
      int mi = i;         /* max index */
      double mv = fabs (a[i*n + i]); /* max value */
      for (int j = i; j < n; j++)
        if (mv < fabs (a[i*n + j]))
          {
            mi = j;
            mv = fabs (a[i*n + j]);
          }
      double buf_swap[4];
      int _i;
      int n2 = n * n;
      int inc = 4 * n;
      for (_i = 0; _i + 3 * n < n2; _i += inc)
        {

          buf_swap[0] = a[_i + i];
          buf_swap[1] = a[_i + n + i];
          buf_swap[2] = a[_i + 2 * n + i];
          buf_swap[3] = a[_i + 3 * n + i];
          a[_i  + i] = a[_i + mi];
          a[_i + n + i] = a[_i + n + mi];
          a[_i + 2 * n + i] = a[_i + 2 * n + mi];
          a[_i + 3 * n + i] = a[_i + 3 * n + mi];
          a[_i  + mi] = buf_swap[0];
          a[_i + n + mi] = buf_swap[1];
          a[_i + 2 * n + mi] = buf_swap[2];
          a[_i + 3 * n + mi] = buf_swap[3];
        }
      int r = n & 3;
      for (int kr = 0; kr < r; kr++)
        {
          buf_swap[0] = a[_i + kr * n + i];
          a[_i + kr * n + i] = a[_i + kr * n + mi];
          a[_i + kr * n + mi] = buf_swap[0];
        }
      if (fabs (mv) < EPS)
        return 1;

      sp[i] = mi;
      mv = 1 / a[i*n + i];
      for (int j = 0; j < n; j++)
        rev[i*n + j] = mv * rev[i*n +j];
      for (int j = i; j < n; j++)
        a[i*n + j] = mv * a[i*n + j];
      for (int k = i + 1; k < n; k++)
        {
          double coef = a[k*n+i];
          int j;
          for (j = 0; j + 3 < n; j += 4)
            {
              rev[k*n + j]     -= coef * rev[i*n + j];
              rev[k*n + j + 1] -= coef * rev[i*n + j + 1];
              rev[k*n + j + 2] -= coef * rev[i*n + j + 2];
              rev[k*n + j + 3] -= coef * rev[i*n + j + 3];
            }
          int r = n & 3;
          for (int kr = 0; kr < r; kr++)
            rev[k*n + j + kr] -= coef * rev[i*n + j + kr];
          for (j = i + 1; j + 3 < n; j += 4)
            {
              a[k*n + j]     -= coef * a[i*n + j];
              a[k*n + j + 1] -= coef * a[i*n + j + 1];
              a[k*n + j + 2] -= coef * a[i*n + j + 2];
              a[k*n + j + 3] -= coef * a[i*n + j + 3];
            }
          r = (n  - i - 1) & 3;
          for (int kr = 0; kr < r; kr++)
            a[k*n + j + kr] -= coef * a[i*n + j + kr];
        }
    }

  for (int i = n - 1; i >= 0; i--)
    {
      for (int _i = i - 1; _i >= 0; _i--)
        {
          double coef = a[_i * n + i];
          int j;
          for (j = 0; j + 3 < n; j += 4)
            {
              rev[_i*n + j] -= coef * rev[i*n + j];
              rev[_i*n + j + 1] -= coef * rev[i*n + j + 1];
              rev[_i*n + j + 2] -= coef * rev[i*n + j + 2];
              rev[_i*n + j + 3] -= coef * rev[i*n + j + 3];
            }
          int r = n & 3;
          for (int kr = 0; kr < r; kr++)
            rev[_i * n + j + kr] -= coef * rev[i*n + j + kr];
        }
    }
  for (int i = n - 1; i >= 0; i--)
    if (sp[i] != i)
      for (int j = 0; j < n; j++)
        swap (rev, sp[i]*n + j, i*n + j);
  return 0;
}


void matrix_product (double *a, double *b, double *c, const int a_rows, const int a_cols, const int b_cols)
{
  double unroll[8];
  double sum[9];
  int i = 0,
      j = 0;
  for (; i < a_rows - 2; i += 3)
    {
      for (j = 0; j < b_cols - 2; j += 3)
        {
          for (int iter = 0; iter < 9; iter++)
            sum[iter] = 0;
          for (int k = 0; k < a_cols; k++)
            {
              sum[0] += a[i*a_cols + k]       * b[k*b_cols +  j];
              sum[1] += a[(i + 1)*a_cols + k] * b[k*b_cols +  j];
              sum[2] += a[(i + 2)*a_cols + k] * b[k*b_cols +  j];
              sum[3] += a[i*a_cols + k]       * b[k*b_cols +  j + 1];
              sum[4] += a[(i + 1)*a_cols + k] * b[k*b_cols +  j + 1];
              sum[5] += a[(i + 2)*a_cols + k] * b[k*b_cols +  j + 1];
              sum[6] += a[i*a_cols + k]       * b[k*b_cols +  j + 2];
              sum[7] += a[(i + 1)*a_cols + k] * b[k*b_cols +  j + 2];
              sum[8] += a[(i + 2)*a_cols + k] * b[k*b_cols +  j + 2];
            }
          for (int j_bl = 0; j_bl < 3; j_bl++)
            for (int i_bl = 0; i_bl < 3; i_bl++)
              c[(i + i_bl) * b_cols + j + j_bl] = sum[j_bl * 3 + i_bl];

        }
    }
  int i_stop = i,
      j_stop = j;
  for (; i < a_rows; i++)
    for (j = 0; j < b_cols; j++)
      {
        for (int iter = 0; iter < 8; iter++)
          unroll[iter] = 0;
        int k;
        for (k = 0; k < a_cols - 7; k += 8)
          {
            unroll[0] += a[i*a_cols + k] * b[k*b_cols +  j];
            unroll[1] += a[i*a_cols + k + 1] * b[(k + 1)*b_cols +  j];
            unroll[2] += a[i*a_cols + k + 2] * b[(k + 2)*b_cols +  j];
            unroll[3] += a[i*a_cols + k + 3] * b[(k + 3)*b_cols +  j];
            unroll[4] += a[i*a_cols + k + 4] * b[(k + 4)*b_cols +  j];
            unroll[5] += a[i*a_cols + k + 5] * b[(k + 5)*b_cols +  j];
            unroll[6] += a[i*a_cols + k + 6] * b[(k + 6)*b_cols +  j];
            unroll[7] += a[i*a_cols + k + 7] * b[(k + 7)*b_cols +  j];
          }
        for (; k < a_cols; k++)
          unroll[0] += a[i*a_cols + k] * b[k*b_cols + j];
        c[i * b_cols + j] = unroll[0] + unroll[1] + unroll[2] + unroll[3] 
         + unroll[4] + unroll[5] + unroll[6] + unroll[7]; 
      }
  for (i = 0; i < i_stop; i++)
    for (j = j_stop; j < b_cols; j++)
      {
        for (int iter = 0; iter < 8; iter++)
          unroll[iter] = 0;

        int k;
        for (k = 0; k < a_cols - 7; k += 8)
          {
            unroll[0] += a[i*a_cols + k] * b[k*b_cols +  j];
            unroll[1] += a[i*a_cols + k + 1] * b[(k + 1)*b_cols +  j];
            unroll[2] += a[i*a_cols + k + 2] * b[(k + 2)*b_cols +  j];
            unroll[3] += a[i*a_cols + k + 3] * b[(k + 3)*b_cols +  j];
            unroll[4] += a[i*a_cols + k + 4] * b[(k + 4)*b_cols +  j];
            unroll[5] += a[i*a_cols + k + 5] * b[(k + 5)*b_cols +  j];
            unroll[6] += a[i*a_cols + k + 6] * b[(k + 6)*b_cols +  j];
            unroll[7] += a[i*a_cols + k + 7] * b[(k + 7)*b_cols +  j];
          }
        for (; k < a_cols; k++)
          unroll[0] += a[i*a_cols + k] * b[k*b_cols + j];
        c[i * b_cols + j] = unroll[0] + unroll[1] + unroll[2] + unroll[3] 
         + unroll[4] + unroll[5] + unroll[6] + unroll[7]; 
      }

}

void matrix_product_rw (double *a, double *b, double *c, const int a_rows, const int a_cols, const int b_cols)
{
  double sum[4];
  int r = a_cols & 3;
  for (int j = 0; j < b_cols; j++)
    {
      int l;
      for (l = 0; l + 3 < a_cols; l += 4)
        {
          c[l] = b[l*b_cols + j];
          c[l+1] = b[(l+1)*b_cols + j];
          c[l+2] = b[(l+2)*b_cols + j];
          c[l+3] = b[(l+3)*b_cols + j];
        }
      for (int kr = 0; kr < r; kr++)
        c[l + kr] = b[(l + kr)*b_cols + j];
      for (int i = 0; i < a_rows; i++)
        {
          sum[0] = 0;
          sum[1] = 0;
          sum[2] = 0;
          sum[3] = 0;
          int k;
          for (k = 0; k + 3 < a_cols; k += 4)
            {
              sum[0] += a[i*a_cols + k] * c[k];
              sum[1] += a[i*a_cols + k + 1] * c[k  +1];
              sum[2] += a[i*a_cols + k + 2] * c[k + 2];
              sum[3] += a[i*a_cols + k + 3] * c[k + 3];
            }
          for (int kr = 0; kr < r; kr++)
            sum[0] += a[i*a_cols +k + kr] * c[k + kr];
          b[i * b_cols + j] = sum[0] + sum[1] + sum[2] + sum[3];
        }
    }
}

void matrix_copy (double *source, double *dest, const int rows, const int cols)
{
  int j;
  for (j = 0; j + 3 < cols * rows; j += 4)
    {
      dest[j] = source[j];
      dest[j + 1] = source[j + 1];
      dest[j + 2] = source[j + 2];
      dest[j + 3] = source[j + 3];
    }
  int r = cols * rows & 3;
  for (int kr = 0; kr < r; kr++)
    dest[j + kr] = source[j + kr];
}

int matrix_read_file(double *a, int n, int m, const char *filename)
{
  FILE *fin = fopen(filename, "r");
  if (!fin)
    return FILE_NOT_FOUND;
  if (n < m)
    {
      fprintf(stderr, "ERROR: m must be less then n in matrix_read_file\n");
      return ERR_LOGGED;
    }
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
        int index = _get_ij_index(n, m, i, j);
        if (1 != fscanf (fin, "%lf", a + index))
          {
            if (ferror (fin))
              {
                fclose (fin);
                return WRONG_FORMAT;
              }
            else
              {
                fclose(fin);
                return UNEXP_EOF;
              }
          }
      }
  return 0;
}

int _block_ij(int n, int m, int i, int j, int *rows, int *cols)   // size of the block is returned through rows and cols
{
  int s = n % m;
  int l = (n - s) / m;
  int ret = 0;
  int _rows = 0;
  int _cols = 0;

  if (i == l && j == l)
    {
      _rows = s;
      _cols = s;
      ret = m * m * l * l + s * m * 2 * l;
    }
  if (i < l && j == l)
    {
      _rows = m;
      _cols = s;
      ret = m * m * l * l + s * m * (l + i);
    }
  if (i == l && j < l)
    {
      _rows = s;
      _cols = m;
      ret = m * m * (l * j + i) + s * m * j;
    }
  if (i < l && j < l)
    {
      _rows = m;
      _cols = m;
      ret = m * m * (l * j + i) + s * m * j;
    }
  if (rows)
    *rows = _rows;
  if (cols)
    *cols = _cols;
  return ret; 
}

int _get_ij_index(int n, int m, int i, int j)
{
  int block_i = i / m,
      block_j = j / m,
      rows, cols,
      block = _block_ij(n, m, block_i, block_j, &rows, &cols);
  int
   block_row = i % m,
             block_col = j % m;
  return block + cols * block_row + block_col;
}


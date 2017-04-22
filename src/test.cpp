#include "../include/utils.h"
#include "../include/matrix.h"
#include "../include/sbc_matrix.h"
#include "../include/init_functions.h"

#include <stdio.h>
#include <stdlib.h>

int main ()
{
  int n = 5;
  int m = 2;
  int p = 1;
  int rank = 0;
  SBC_storage st;
  SBC_init (&st, n, m, p, rank, consecutive);
  SBC_print (st.loc_matrix, n, m, stdout);
  double *coef_matrix_array = (double *) malloc ((n - m) * m * sizeof (double));
  double *subtrahend_ = (double *) malloc (m * m * sizeof (double));
  for (int i = 0; i < (n - m) * m; i++)
    coef_matrix_array[i] = st.loc_matrix[i + m * m];
  SBC_gauss_subtr_all (&st, coef_matrix_array, subtrahend_, 0);
  SBC_print (st.loc_matrix, n, m, stdout);
  SBC_destroy (&st);
}

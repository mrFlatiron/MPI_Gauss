#pragma once

void   matrix_sq_init_e(double *a, const int n);

void   matrix_product(double *a, double *b, double *c, const int a_rows, const int a_cols, const int b_cols);

void   matrix_product_rw(double *a, double *b, double *c, const int a_rows, const int a_cols, const int b_cols);


int    matrix_reverse(double *a, double *rev, int *sp, const int n); /*sp (int[n]) - buffer for internal algorithm permutations*/

void   matrix_copy(double *source, double *dest, const int rows, const int cols);

void   matrix_subtr_rw(double *a, double *b, const int rows, const int cols);

double matrix_sq_norm(double *a, const int n);
double matrix_norm(double *a, const int rows, const int cols);

int    matrix_read_file(double *a, int n, int m, const char *filename);

int _block_ij(int n, int m, int i, int j, int *rows, int *cols);
int _get_ij_index(int n, int m, int i, int j);



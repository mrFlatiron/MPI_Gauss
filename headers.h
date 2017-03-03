#pragma once
#include <math.h>
#include "const.h"
#include <stdlib.h>
#include <ctime>
#include <sstream>
#include <string>
#include <string.h>
#include <cstdio>

using namespace std;

void   matrix_init(double *a, int n, int m);
void   matrix_init_e(double *a, const int n);
void   matrix_product(double *a, double *b, double *c, const int a_rows, const int a_cols, const int b_cols);
int    matrix_read_file(double *a, int n, int m, const char *filename);
void   matrix_print(double *a, int n, int m, const int corner);
int    matrix_reverse(double *a, double *rev, int *sp, const int n);
void   matrix_copy(double *source, double *dest, const int rows, const int cols);
void   matrix_subtr_rw(double *a, double *b, const int rows, const int cols);
void   matrix_vector_product(double *a, double *x, double *b, const int n, const int m);
double matrix_norm(double *a, const int n);
double matrix_norm(double *a, const int rows, const int cols);

void swap(double *a, const int i, const int j);
int _block_ij(int n, int m, int i, int j, int *rows, int *cols);
int _get_ij_index(int n, int m, int i, int j);

int read_doubles(double *a, const int size, const char *filename);
void print_doubles(double *a, const int size);

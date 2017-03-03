#pragma once
#include <stdio.h>

void swap (double *a, const int i, const int j);
void swap_arrays (double *a, double *b, const int n);

void log_error (const int rank, const char *msg, FILE* = stderr);
void log_info (const int rank, const char *msg, FILE* = stdout);

#pragma once

#include "const.h"
#include <mpi.h>
#include <stdio.h>


/*
   SBC - Stored By Column

   VERSION
   1.0
   DESCRIPTION This is a C++ library for manipulating a square matrix in MPI programs with a storage, organized by block columns.  
   CRUCIAL CONVENTION
   It is assumed that a RANK process stores all block columns with 
   numbers RANK + K * P, where K is an integer number.

   VARIABLES NAMES
   p - size of communicator
   rank - MPI rank of a process
   n - size of global matrix
   m - size of block
   l = n / m
   s = n - m * l

   COMMENT
   В начале хотел на си делать, поэтому вместо класса - префикс и т. д.

*/

#define SBC_STORAGE_CRIPPLED    -1
#define SBC_STORAGE_DECLARED    0
#define SBC_STORAGE_ALLOCATED   1 
#define SBC_STORAGE_INITIALIZED 2


struct pivot_c_struct
{
  int glob_pivot_index;
  double norm_of_rev;
};



struct SBC_storage
{
  int state;
  int number_of_m_cols;
  int has_s_col;
  int total_size;
  int n;
  int m;
  int l;
  int s;
  int p;
  int rank;
  double *loc_matrix;
  double *b;
};

void SBC_pivot_reduce (
                       pivot_c_struct *in, 
                       pivot_c_struct *inout,
                       int *len,
                       MPI_Datatype *SBC_pivot_candidate
                      );

void SBC_pivot_op_create (
                          MPI_Op *SBC_pivot_op
                         );

void SBC_pivot_candidate_create (
                                 MPI_Datatype *new_type
                                );

void SBC_set_meta (
                   SBC_storage *loc_storage,
                   const int n,
                   const int m,
                   const int p, 
                   const int rank
                  );

void SBC_allocate (
                    SBC_storage *loc_storage,
                    const int n,
                    const int m,
                    const int p,
                    const int rank
                   );


void SBC_init (
               SBC_storage *loc_storage, 
               const int n, 
               const int m, 
               const int comm_size,
               const int rank,
               double (*f) (const int glob_i, const int glob_j, const int n)
              ); 

void SBC_fill_in (
                 SBC_storage *loc_storage,
                 double (*f) (const int glob_i, const int glob_j, const int n)
                );

void SBC_MPI_init (
                   SBC_storage *loc_storage,
                   const int n, 
                   const int m, 
                   const int p,
                   const int rank,
                   MPI_Comm comm,
                   double *glob_matrix,
                   const char *file
                  );

void SBC_MPI_fill_in (
                   SBC_storage *loc_storage,
                   MPI_Comm comm,
                   double *glob_matrix,
                   const char *file
                   );


void SBC_simple_print (
                double *glob_matrix,
                const int n,
                const int m,
                const int corner,
                FILE *fout
               );  /* Ineffective. Use only for occasional logging */

void SBC_MPI_print (
                    SBC_storage *loc_storage,
                    double *glob_matrix,
                    double *buf_matrix,
                    MPI_Comm comm,
                    const int corner = DEF_PRINT_SIZE,
                    FILE *fout = stdout
                   ); /*Ineffective */

int SBC_index (
               const int i, 
               const int j, 
               const int n, 
               const int m
              ); /* Ineffective. Use only for occasional logging */

int SBC_get_loc_size (
                      const int p, 
                      const int rank, 
                      const int n, 
                      const int m
                     );

void SBC_place_in_glob (
                        SBC_storage *loc_storage,
                        double *glob_matrix
                       );


void SBC_gauss_find_candidate_pivot ( 
                                     SBC_storage *loc_storage,
                                     const int i_main,
                                     int *glob_candidate_pivot,
                                     double *norm_of_rev,
                                     double *buf, //2 * m * m size
                                     int *permutation //m size
                                    );

void SBC_gauss_multi_row (
                          SBC_storage *loc_storage,
                          double *coef_matrix,
                          double *buf_matrix, //m * m size
                          const int i_main
                         );
void SBC_gauss_subtr_all (
                          SBC_storage *loc_storage,
                          double *coef_matrix_array,
                          double *subtrahend_, //m * m size
                          const int i_main
                         );

void SBC_gauss_subtr_backward (
                               SBC_storage *loc_storage,
                               const int i_main,
                               double *coef_col,
                               double *buf //m * m size
                              );


void SBC_gauss_MPI_swap (
                         SBC_storage *loc_storage,
                         const int pivot,
                         const int i_main,
                         MPI_Comm comm
                        );

void SBC_gauss_MPI_find_pivot (
                               SBC_storage *loc_storage,
                               const int i_main,
                               int *glob_pivot,
                               MPI_Comm comm,
                               MPI_Datatype SBC_pivot_candidate,
                               MPI_Op SBC_pivot_op,
                               double *buf, //2 * m * m size
                               int *permutation
                              );

void SBC_gauss_MPI_multi_row (
                              SBC_storage *loc_storage,
                              const int i_main,
                              MPI_Comm comm,
                              double *buf_, //2 * m * m size
                              int *permutation
                             );


void SBC_gauss_MPI_subtr_all (
                              SBC_storage *loc_storage,
                              const int i_main,
                              MPI_Comm comm,
                              double *buf_ // n * m + m * m size
                             );

void SBC_gauss_MPI_backward (
                             SBC_storage *loc_storage,
                             const int i_main,
                             MPI_Comm comm,
                             double *buf_ //n * m + m * m size
                            );

void SBC_MPI_multi_vector (
                           SBC_storage *loc_storage,
                           double *x,
                           double *b_save, //out
                           MPI_Comm comm,
                           double *buf //m * m  + n size
                          );
void SBC_destroy (
                  SBC_storage *loc_storage
                 );

#include "mpi.h"
#include "../include/utils.h"
#include "../include/matrix.h"
#include "../include/sbc_matrix.h"
#include "../include/init_functions.h"
#include "../include/const.h"

#include <stdlib.h>

#include <stdio.h>


int main (int argc, char *argv[])
{

double (*init_function) (const int, const int, const int)
 = I;


  int n;
  int m; 
  int l;
  int p;
  int rank;
  int glob_pivot;


  SBC_storage  loc_storage; MPI_Datatype SBC_pivot_candidate;
  MPI_Op       SBC_pivot_op;



  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &p);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  SBC_pivot_candidate_create (&SBC_pivot_candidate);
  SBC_pivot_op_create (&SBC_pivot_op);

  if (argc < 3)
    {
      if (rank == 0) 
        fprintf (stderr, "Usage: %s <size> <size of block> [matrix file]\n", argv[0]);

      MPI_Finalize ();
      return 0;
    }

  n  = atoi (argv[1]);
  m = atoi (argv[2]);
  l = n / m;

  if (!n || !m)
    {
      if (rank == 0)
        fprintf (stderr, "size and size of block must be integer > 0\n");

      MPI_Finalize ();
      return 0;
    }

  if (m > n)
    {
      if (rank == 0)
        fprintf (stderr, "m must be less or equal  than n\n");

      MPI_Finalize ();
      return 0;
    }


  double *buf_              = new double[n * m + m * m];
  double *x                 = new double[n];
  double *b_save            = new double[n];
  double *b_untouchable     = new double[n];
  int *int_buf              = new int[m];
  int *permutation          = new int[n];

  double *b;
  double *glob_matrix;
  double *buf_print;

  if (rank == 0)
    {
      glob_matrix = new double[n * n];
      buf_print = new double [(l / p + 1) * n * m];
    }
  else
    {
      glob_matrix = NULL;
      buf_print = NULL;
    }

  if (!buf_ ||
      !x ||
      !int_buf ||
      !permutation ||
      !b_save ||
      !b_untouchable ||
      (rank == 0 && !glob_matrix))
    {
      MPI_Abort (MPI_COMM_WORLD, BAD_ALLOC);
      return 0;
    }

  if (argc >= 4)
    SBC_MPI_init (&loc_storage, 
                  n, m, p, rank, 
                  MPI_COMM_WORLD,
                  glob_matrix, 
                  argv[3]);
  else
    SBC_init (&loc_storage, n, m, p, rank, init_function); 

  b = loc_storage.b;

  if (rank == 0)  
    matrix_copy (b, b_untouchable, n, 1);

  if (loc_storage.state == SBC_STORAGE_CRIPPLED)
    {
      fprintf (stderr, "[RANK]: %d -> [ERROR]: Allocation failure\n", rank);
      MPI_Abort (MPI_COMM_WORLD, -1);
      return 0;
    }

  if (loc_storage.state == SBC_STORAGE_DECLARED)
    {
      fprintf (stderr, "[RANK]: %d -> [ERROR]: No job for the process. Adjust the number of processes\n", rank);
      MPI_Abort (MPI_COMM_WORLD, -1);
      return 0;
    }

  if (loc_storage.state < SBC_STORAGE_INITIALIZED)
    {
      if (rank == 0)
        fprintf (stderr, "Matrix is not initialized. Nothing to solve\n");
      MPI_Abort (MPI_COMM_WORLD, -1);
      return 0;
    }

  fflush (stdout);

  if (rank == 0)
    printf ("=================================BEGIN RESOLVING===================================\n");


  SBC_MPI_print (&loc_storage, glob_matrix, buf_print, MPI_COMM_WORLD);

  MPI_Barrier (MPI_COMM_WORLD);
  double start = MPI_Wtime ();

  for (int i_main = 0; i_main < l; i_main++)
    {

      if (rank == 0)
        printf ("iteration = %d\n", i_main);

      SBC_gauss_MPI_find_pivot (
                                &loc_storage,
                                i_main,
                                &glob_pivot,
                                MPI_COMM_WORLD,
                                SBC_pivot_candidate,
                                SBC_pivot_op,
                                buf_,
                                int_buf
                               ); // all_reduce 
      if (glob_pivot == -1)
        {
          delete[] buf_;
          delete[] x;
          delete[] b_save;
          delete[] b_untouchable;
          delete[] int_buf;
          delete[] permutation;
          if (glob_matrix)
            delete[] glob_matrix;
          if (buf_print)
            delete[] buf_print;
          SBC_destroy (&loc_storage);

          if (rank == 0)
            fprintf (stderr, "[ERROR]: NO SOLUTION\n");

          MPI_Finalize ();
          return NO_SOLUTION;
        }
      else 
        permutation[i_main] = glob_pivot;

      SBC_gauss_MPI_swap (
                          &loc_storage,
                          glob_pivot,
                          i_main,
                          MPI_COMM_WORLD
                         );  //sendrecv_replace 
      SBC_gauss_MPI_multi_row (
                               &loc_storage,
                               i_main,
                               MPI_COMM_WORLD,
                               buf_,
                               int_buf
                              ); //bcast

      SBC_gauss_MPI_subtr_all (
                               &loc_storage,
                               i_main,
                               MPI_COMM_WORLD,
                               buf_
                              );  //bcast

    }


  SBC_gauss_MPI_multi_row (
                           &loc_storage,
                           l,
                           MPI_COMM_WORLD,
                           buf_,
                           int_buf
                          ); //send

  for (int i_main = l; i_main > 0; i_main --)
    SBC_gauss_MPI_backward (
                            &loc_storage,
                            i_main,
                            MPI_COMM_WORLD,
                            buf_
                           ); //send


  if (rank == 0)
    {
      for (int i = l - 1; i >= 0; i--)
        if (permutation[i] != i)
          swap_arrays(b + permutation[i] * m, b + i * m, m);

      matrix_copy (b, x, n, 1);
    }

  if (rank == 0)
    printf ("=================================FINISHED RESOLVING================================\n");


  if (argc < 4)
    SBC_fill_in (&loc_storage, init_function);
  else 
    SBC_MPI_fill_in (&loc_storage, MPI_COMM_WORLD, glob_matrix, argv[3]);

  SBC_MPI_multi_vector (&loc_storage, x, b_save, MPI_COMM_WORLD, buf_);

  if (rank == 0)
    {
      matrix_subtr_rw (x, b_untouchable, n, 1);
      double norm = 
       matrix_norm (x, n, 1 );
      printf ("discrepancy = %e\n", norm);
    }
  double loc_end = MPI_Wtime ();
  MPI_Barrier (MPI_COMM_WORLD);
  double glob_end = MPI_Wtime ();

  if (rank == 0)
    printf ("GLOB WTIME = %f sec\n", glob_end - start);

  printf ("RANK %d WTIME = %f sec\n", rank, loc_end - start);

  delete[] buf_;
  delete[] x;
  delete[] b_save;
  delete[] b_untouchable;
  delete[] int_buf;
  delete[] permutation;
  if (glob_matrix)
    delete[] glob_matrix;
  if (buf_print)
    delete[] buf_print;

  SBC_destroy (&loc_storage);
  MPI_Finalize ();
  return 0;
}

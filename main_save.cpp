#include "utils.h"
#include "matrix.h"
#include "sbc_matrix.h"
#include "init_functions.h"

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main (int argc, char *argv[])
{

  MPI_Init(&argc, &argv);

  int n;
  int m;
  int p;
  int rank;

  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc < 3)
    {
      if (rank == 0) 
        fprintf (stderr, "Usage: %s <size> <size of block> [matrix file]\n", argv[0]);
      MPI_Finalize();
      return 0;
    }

  n  = atoi (argv[1]);
  m = atoi (argv[2]);

  if (!n || !m)
    {
      if (rank == 0)
        fprintf (stderr, "size and size of block must be integer > 0\n");
      MPI_Finalize();
      return 0;
    }


  SBC_storage loc_storage;
  SBC_init (&loc_storage, n, m, p, rank, consecutive); 

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


  int loc_size = loc_storage.total_size;
  
  printf ("[RANK]: %d -> [INFO]: Good init. Size = %d!\n", rank, loc_size);

/*
********************************************
*/
/*
  double *buf_multi_matrix = (double *) malloc (m * m * sizeof (double));
  
  double *buf_save_matrix = (double *) malloc (m * m * sizeof (double));

  int *permutation = (int *) malloc (m * sizeof (int));

  int glob_candidate_pivot;
  double loc_norm_of_rev;
  
  printf ("[%d]->[INFO]: pos = %d | val = %lf\n", rank, glob_candidate_pivot, loc_norm_of_rev);

  double *coef_matrix_array = (double *) malloc ((n - m) * m * sizeof (double));

  for (int i = m * m; i < n * m; i++)
    coef_matrix_array[i - m * m] = loc_storage.loc_matrix[i];
*/

  SBC_gauss_MPI_swap (&loc_storage, 0, 1, MPI_COMM_WORLD);

  if (rank == 0)
    {
      double *glob_matrix = (double *) malloc (n * n * sizeof (double));
      double *buf_matrix = (double *) malloc ((loc_size + n * m) * sizeof (double));

      SBC_place_in_glob (glob_matrix, loc_storage.loc_matrix, p, 0, n, m);
      MPI_Status st; 
     

      for (int i = 1; i < p; i++)
        {
          int out_size = SBC_get_loc_size (p, i, n , m); 
          MPI_Recv (buf_matrix, out_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &st); 
          SBC_place_in_glob (glob_matrix, buf_matrix, p, i, n, m);
        }
     

      SBC_print (glob_matrix, n, m, stdout); 


      free (glob_matrix);
      free (buf_matrix);
    }
  else
    MPI_Send (loc_storage.loc_matrix, loc_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);


  SBC_destroy (&loc_storage);
  MPI_Finalize ();
  return 0;
}

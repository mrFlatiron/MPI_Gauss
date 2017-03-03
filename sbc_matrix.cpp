#include "sbc_matrix.h"
#include "matrix.h"
#include "utils.h"
#include "const.h"
#include <mpi.h>
#include <stdlib.h>

void SBC_pivot_candidate_create (MPI_Datatype *SBC_pivot_candidate)
{
  MPI_Datatype type[2] = {MPI_INT, MPI_DOUBLE};
  pivot_c_struct pivot_candidate;
  int          blocklen[2] = {1, 1};
  MPI_Aint     disp[2];
  MPI_Aint base;


  MPI_Get_address (&pivot_candidate, &base);
  MPI_Get_address (&pivot_candidate.glob_pivot_index, disp); 
  MPI_Get_address (&pivot_candidate.norm_of_rev, disp + 1); 
  disp[0] = disp[0] - base;
  disp[1] = disp[1] - base - disp[0];

  MPI_Type_struct (2, blocklen, disp, type, SBC_pivot_candidate);

  MPI_Type_commit (SBC_pivot_candidate);
}

void SBC_pivot_reduce (
                       pivot_c_struct *in, 
                       pivot_c_struct *inout,
                       int *len,
                       MPI_Datatype *SBC_pivot_candidate
                      )
{
  if (in->glob_pivot_index >= 0)
    {
      if (inout->glob_pivot_index < 0 || inout->norm_of_rev > in->norm_of_rev)
        {
          inout->glob_pivot_index = in->glob_pivot_index;
          inout->norm_of_rev = in->norm_of_rev;
        }
    }
}



void SBC_pivot_op_create (
                          MPI_Op *SBC_pivot_op
                         )
{
  MPI_Op_create ((MPI_User_function *)SBC_pivot_reduce, 1, SBC_pivot_op);
}

void SBC_set_meta (
                   SBC_storage *loc_storage,
                   const int n,
                   const int m,
                   const int p,
                   const int rank
                  )
{
  int l = n / m;
  int s = n - m * l;


  int loc_size = SBC_get_loc_size (p, rank, n , m);

  loc_storage->l = l;
  loc_storage->s = s;

  loc_storage->n = n;
  loc_storage->m = m;

  loc_storage->p = p;
  loc_storage->rank = rank;

  loc_storage->total_size = loc_size;

  loc_storage->number_of_m_cols = (rank < l % p) ? l/p +1 : l/p;
  loc_storage->has_s_col = (l % p == rank && s != 0) ? 1 : 0;

  loc_storage->state = SBC_STORAGE_DECLARED;
}

void SBC_allocate (
                    SBC_storage *loc_storage,
                    const int n,
                    const int m,
                    const int p,
                    const int rank
                   )
{
  SBC_set_meta (loc_storage, n, m, p, rank);
  
  int loc_size = loc_storage->total_size;
  double *loc_matrix = loc_storage->loc_matrix;
  double *b = loc_storage->b;
  
    if (!loc_size)
    return;

  if (loc_size != 0)
    {
      if (rank == 0)
        b = new double[n];
      else 
        b = NULL;

      loc_matrix = new double[loc_size];
    }

  loc_storage->loc_matrix = loc_matrix;
  loc_storage->b = b;


  if (!loc_matrix || (rank == 0 && !b))
    {
      loc_storage->state = SBC_STORAGE_CRIPPLED;  
      if (loc_matrix)
        delete[] loc_matrix;
      if (b) 
        delete[] b;
      b = NULL;
      return;
    }

  loc_storage->state = SBC_STORAGE_ALLOCATED;
}

void SBC_init (
               SBC_storage *loc_storage, 
               const int n, 
               const int m, 
               const int p,
               const int rank,
               double (*f) (const int glob_i, const int glob_j, const int n)
              ) 
{
  loc_storage->state = SBC_STORAGE_DECLARED;

  SBC_allocate (loc_storage, n, m, p, rank);

  if (loc_storage->state < SBC_STORAGE_ALLOCATED)
    {
      loc_storage->state = SBC_STORAGE_CRIPPLED;
      return ;
    }

  SBC_fill_in (loc_storage, f);
  
  loc_storage->state = SBC_STORAGE_INITIALIZED;
} 

void SBC_MPI_init (
                   SBC_storage *loc_storage,
                   const int n, 
                   const int m, 
                   const int p,
                   const int rank,
                   MPI_Comm comm,
                   double *glob_matrix,
                   const char *file
                  )
{
  loc_storage->state = SBC_STORAGE_DECLARED;

  SBC_allocate (loc_storage, n, m, p, rank);

  if (loc_storage->state < SBC_STORAGE_ALLOCATED)
    {
      loc_storage->state = SBC_STORAGE_CRIPPLED;
      return ;
    }

  int l = n / m;
  int s = n - m * l; 
  double *loc_matrix = loc_storage->loc_matrix;
  int number_of_m_cols = loc_storage->number_of_m_cols;
  int has_s_col = loc_storage->has_s_col;

  int nm = n * m;
  int ns = s * n;

  int dest;

  double *to_send = glob_matrix;
  double *to_recv = loc_matrix;


  if (rank == 0)
    {
      int error = matrix_read_file (glob_matrix, n, m, file);
      if (error)
        {
          loc_storage->state = SBC_STORAGE_ALLOCATED;
          fprintf (stderr, "Error occured during reading the file. Code : %d\n", error);
          return ;
        }
      for (int i = 0; i < l; i++)
        {
          dest = i % p;
          if (dest == 0)
            {
              matrix_copy (to_send, to_recv, n, m); 
              to_recv += nm;
            }
          else
            MPI_Send (to_send, nm, MPI_DOUBLE, dest, 0, comm);
          to_send += nm;
        }
      dest = l % p;
      if (dest == 0)
          matrix_copy (to_send, to_recv, n, s); 
      else
        MPI_Send (to_send, ns, MPI_DOUBLE, dest, 0, comm);
      to_send += nm;
    }
  else 
    {
      MPI_Status st;
      for (int i = 0; i < number_of_m_cols; i++)
        {
          MPI_Recv (to_recv, nm, MPI_DOUBLE, 0, 0, comm, &st);
          to_recv += nm;
        }
      if (has_s_col)
        MPI_Recv (to_recv, ns, MPI_DOUBLE, 0, 0, comm, &st);
    }

  loc_storage->state = SBC_STORAGE_INITIALIZED;
}


void SBC_fill_in (
                  SBC_storage *loc_storage,
                  double (*f) (const int glob_i, const int glob_j, const int n)
                 )
{
  if (loc_storage->state < SBC_STORAGE_ALLOCATED)
    {
      printf ("[ERROR}: wrong use of SBC_fill_in\n");
      return;
    }
  int rank = loc_storage->rank;
  int n = loc_storage->n;
  int m = loc_storage->m;
  int l = loc_storage->l;
  int p = loc_storage->p;
  int s = loc_storage->s;
  double *loc_matrix = loc_storage->loc_matrix;
  double *b = loc_storage->b;

  loc_storage->state = SBC_STORAGE_INITIALIZED; /* Not yet, but it will be by the end of the function */

  int glob_block_col = rank;
  int loc_block_col = 0;
  int glob_i;
  int glob_j;
  int index = 0;

  for (; glob_block_col < l; glob_block_col += p, loc_block_col ++)
    {
      glob_i = 0;
      for (int block = 0; block < l; block++)
        {
          for (int loc_inblock_i = 0; loc_inblock_i < m; loc_inblock_i ++, glob_i ++)
            {
              glob_j = glob_block_col * m;
              for (int loc_inblock_j = 0; loc_inblock_j < m; loc_inblock_j ++, glob_j ++)
                {
                  loc_matrix[index] = f (glob_i, glob_j, n); 
                  index++;
                }
            }
        }
      for (int loc_inblock_i = 0; loc_inblock_i < s; loc_inblock_i ++, glob_i ++)
        {
          glob_j = glob_block_col * m;
          for (int loc_inblock_j = 0; loc_inblock_j < m; loc_inblock_j ++, glob_j ++)
            {
              loc_matrix[index] = f (glob_i, glob_j, n); 
              index++;
            }
        }
    } 

  if (index < loc_storage->total_size)
    {

      glob_i = 0;

      for (int block = 0; block < l; block ++)
        {
          for (int loc_inblock_i = 0; loc_inblock_i < m; loc_inblock_i ++, glob_i ++)
            {
              glob_j = l * m;
              for (int loc_inblock_j = 0; loc_inblock_j < s; loc_inblock_j ++, glob_j ++)
                {
                  loc_matrix[index] = f (glob_i, glob_j, n);
                  index++;
                }
            }
        }

      glob_i = l * m;

      for (int loc_inblock_i = 0; loc_inblock_i < s; loc_inblock_i ++, glob_i ++)
        {
          glob_j = l * m;
          for (int loc_inblock_j = 0; loc_inblock_j < s; loc_inblock_j ++, glob_j ++)
            {
              loc_matrix[index] = f (glob_i, glob_j, n);
              index++;
            }
        }
    }
  if (b)
    {
      for (int glob_i = 0; glob_i < n; glob_i ++)
        {
          b[glob_i] = 0;
          for (int glob_j = 0; glob_j < n; glob_j ++)
            b[glob_i] += (glob_j & 1) ? f (glob_i, glob_j, n) : 0;
        }
    }
}


int SBC_index (const int i, const int j, const int n, const int m)
{
  int l = n / m;
  int s = n - m * l;

  int block_i = i / m;  
  int block_j = j / m;
  int inblock_i = i % m;
  int inblock_j = j % m;
  int block = block_j * (n * m);
  block += (block_j != l) ? block_i * (m * m) : block_i * (s * m);
  int cols = (block_j != l) ? m : s;
  return block + inblock_i * cols + inblock_j;
}

void SBC_simple_print (double *glob_matrix,
                       const int n,
                       const int m,
                       const int corner,
                       FILE *fout)
{
  int i_size = (n < corner) ? n : corner;
  int j_size = (n < corner) ? n : corner;
  for (int i = 0; i < i_size; i++)
    {
      for (int j = 0; j < j_size; j++)
        {
          fprintf (fout, "%.3f\t", glob_matrix[SBC_index (i, j , n , m)]);
        }
      fprintf (fout, "\n");
    }

}

void SBC_MPI_print (
                    SBC_storage *loc_storage,
                    double *glob_matrix,
                    double *buf_matrix,
                    MPI_Comm comm,
                    const int corner,
                    FILE *fout
                   )
{
  int rank = loc_storage->rank;
  int n = loc_storage->n;  
  int m = loc_storage->m;
  int p = loc_storage->p;

  if (rank == 0)
    {
      MPI_Status st; 
      SBC_storage temp;


      SBC_place_in_glob (loc_storage, glob_matrix);

      for (int i = 1; i < p; i++)
        {
          SBC_set_meta (&temp, n, m, p, i);
          int out_size = temp.total_size; 
          MPI_Recv (buf_matrix, out_size, MPI_DOUBLE, i, 0, comm, &st); 
          temp.loc_matrix = buf_matrix;
          SBC_place_in_glob (&temp, glob_matrix);
        }
    }
  else
    MPI_Send (loc_storage->loc_matrix, loc_storage->total_size, MPI_DOUBLE, 0, 0, comm);

  if (rank == 0)
    SBC_simple_print (glob_matrix, n, m, corner, fout);

}


int SBC_get_loc_size (const int p, const int rank, const int n, const int m)
{
  int l = n / m;
  int s = n - m * l;
  int loc_size = l / p * m * n;
  loc_size += (l % p > rank) ? n * m : 0;
  loc_size += (l % p == rank) ? n * s : 0;
  return loc_size;
}

void SBC_place_in_glob (
                        SBC_storage *loc_storage,
                        double *glob_matrix
                       )
{
  int rank = loc_storage->rank;
  int n = loc_storage->n;
  int m = loc_storage->m;
  int p = loc_storage->p;
  int l = loc_storage->l;
  int s = loc_storage->s;
  double *loc_matrix = loc_storage->loc_matrix;

  int glob_index = rank * n * m;
  int loc_index = 0;

  int glob_col = rank;
  int loc_col = 0;

  for (; glob_col < l; glob_col += p, loc_col ++)
    {
      glob_index = glob_col * n * m;
      for (int iter = 0; iter < n * m; iter ++, loc_index ++, glob_index ++)
        {
          glob_matrix[glob_index] = loc_matrix[loc_index];
        }
    }

  if (glob_col == l)
    {
      glob_index = l * n * m;
      for (int iter = 0; iter < n * s; iter ++, loc_index ++, glob_index ++)
        {
          glob_matrix[glob_index] = loc_matrix[loc_index];
        }
    }

}

void SBC_gauss_find_candidate_pivot ( 
                                     SBC_storage *loc_storage,
                                     const int i_main,
                                     int *glob_candidate_pivot,
                                     double *norm_of_rev,
                                     double *buf,
                                     int *permutation
                                    )


{
  int rank = loc_storage->rank;
  int n = loc_storage->n;
  int m = loc_storage->m;
  int l = loc_storage->l;
  int p = loc_storage->p;
  int mm = m * m;
  int nm = n * m;

  double *loc_matrix = loc_storage->loc_matrix;
  double *buf_rev_matrix = buf;
  double *buf_save_matrix = buf + mm;


  int loc_col = i_main / p;
  if (rank < i_main % p)
    loc_col ++;


  double *block = loc_matrix + loc_col * nm + i_main * mm;
  double min_norm = -1;
  int glob_min_col_index = -1;
  double temp_norm;



  for (int glob_col = rank + loc_col * p; glob_col < l; glob_col += p, loc_col ++, block += nm)
    { 
      matrix_copy (block, buf_save_matrix, m, m);
      if (!matrix_reverse (buf_save_matrix, buf_rev_matrix, permutation, m))
        {
          temp_norm = matrix_sq_norm (buf_rev_matrix, m);
          if (temp_norm < min_norm || glob_min_col_index == -1)
            {
              min_norm = temp_norm;
              glob_min_col_index = glob_col;
            }
        }
    }

  *norm_of_rev = min_norm;
  *glob_candidate_pivot = glob_min_col_index;
}

void SBC_gauss_multi_row (
                          SBC_storage *loc_storage,
                          double *coef_matrix,
                          double *buf_matrix,
                          const int i_main
                         )
{
  int rank = loc_storage->rank;
  int n = loc_storage->n;
  int m = loc_storage->m;
  int s = loc_storage->s;
  int p = loc_storage->p;
  int l = loc_storage->l;
  int mm = m * m;
  int nm = n * m;
  double *loc_matrix = loc_storage->loc_matrix;
  double *b = loc_storage->b;
  int has_s_col = loc_storage->has_s_col;
  int number_of_m_cols = loc_storage->number_of_m_cols;

  if (i_main == l && s != 0)
    {
      if (rank == 0)
        {
          matrix_product (coef_matrix, b + l * m, buf_matrix, s, s, 1);
          matrix_copy (buf_matrix, b + l * m, s, 1);
        }
      return ;
    }

  int loc_col = i_main / p;
  if (rank <= i_main % p)
    loc_col ++;

  if (b)
    {
      matrix_product (coef_matrix, b + i_main * m, buf_matrix, m, m, 1); 
      matrix_copy (buf_matrix, b + i_main * m, m, 1);
    }

  double *block = loc_matrix + loc_col * nm; 
  block += i_main * mm; 

  for (; loc_col < number_of_m_cols; loc_col ++, block += nm)
    {
      matrix_product (coef_matrix, block, buf_matrix, m, m, m);
      matrix_copy (buf_matrix, block, m, m);
    }
  if (has_s_col)
    block = loc_matrix +  number_of_m_cols * nm + i_main * s * m; 
  else 
    return;

  matrix_product (coef_matrix, block, buf_matrix, m, m, s); 
  matrix_copy (buf_matrix, block, m, s);
}


void SBC_gauss_subtr_all (
                          SBC_storage *loc_storage,
                          double *coef_matrix_array,
                          double *subtrahend_,
                          const int i_main
                         )
{
  int rank = loc_storage->rank;
  int n = loc_storage->n;
  int m = loc_storage->m;
  int s = loc_storage->s;
  int p = loc_storage->p;
  int l = loc_storage->l;
  double *loc_matrix = loc_storage->loc_matrix;
  double *b = loc_storage->b;
  int has_s_col = loc_storage->has_s_col;
  int number_of_m_cols = loc_storage->number_of_m_cols;

  int loc_col = i_main / p;
  if (rank <= i_main % p)
    loc_col ++;


  int mm = m * m;
  int nm = n * m;
  int sm = s * m;

  double *subtrahend = loc_matrix + loc_col * nm + i_main * mm;
  double *minuend;
  double *coef_matrix;

  for (; loc_col < number_of_m_cols; loc_col ++)
    {
      minuend = subtrahend + mm;
      coef_matrix = coef_matrix_array;
      for (int i = i_main + 1; i < l; i ++)
        {
          matrix_product (coef_matrix, subtrahend, subtrahend_, m, m, m);
          matrix_subtr_rw (minuend, subtrahend_, m, m);
          minuend += mm;
          coef_matrix += mm;
        }
      matrix_product (coef_matrix, subtrahend, subtrahend_, s, m, m);
      matrix_subtr_rw (minuend, subtrahend_, s, m);
      subtrahend += nm;
    }

  if (b)
    {
      coef_matrix = coef_matrix_array;
      subtrahend = b + i_main * m;
      minuend = b + i_main * m + m;

      for (int i = i_main + 1; i < l; i++)
        {
          matrix_product (coef_matrix, subtrahend, subtrahend_, m, m, 1);
          matrix_subtr_rw (minuend, subtrahend_, m, 1);
          minuend += m;
          coef_matrix += mm;
        }


      matrix_product (coef_matrix, subtrahend, subtrahend_, s, m , 1);
      matrix_subtr_rw (minuend, subtrahend_, s, 1);

    }


  if (!has_s_col)
    return;

  subtrahend = loc_matrix + number_of_m_cols * nm + i_main * sm;
  minuend = subtrahend + sm;
  coef_matrix = coef_matrix_array;

  for (int i = i_main + 1; i < l; i++)
    {
      matrix_product (coef_matrix, subtrahend, subtrahend_, m, m, s);

      matrix_subtr_rw (minuend, subtrahend_, m, s);
      minuend += sm;
      coef_matrix += mm;
    }



  matrix_product (coef_matrix, subtrahend, subtrahend_, s, m, s);
  matrix_subtr_rw (minuend, subtrahend_, s, s);

}

void SBC_gauss_MPI_swap (
                         SBC_storage *loc_storage,
                         const int pivot,
                         const int i_main,
                         MPI_Comm comm
                        )
{
  int rank = loc_storage->rank;
  int n = loc_storage->n;
  int m = loc_storage->m;
  int p = loc_storage->p;
  double *loc_matrix = loc_storage->loc_matrix;
  int nm = n * m;

  double *loc_matrix_pivot;
  double *loc_matrix_i_main;
  double *loc_matrix_swap;

  int rank_pivot = pivot % p;
  int rank_i_main = i_main % p;

  int loc_pivot;
  int loc_i_main;
  int loc_swap;

  int dest;

  if (rank != rank_pivot && rank != rank_i_main)
    return;

  if (rank_pivot == rank_i_main)
    {
      loc_pivot = pivot / p;
      loc_i_main = i_main / p; 

      loc_matrix_pivot = loc_matrix + nm * loc_pivot;
      loc_matrix_i_main = loc_matrix + nm * loc_i_main;

      swap_arrays (loc_matrix_pivot, loc_matrix_i_main, nm); 
      return ;
    }
  MPI_Status status;

  dest = rank_pivot + rank_i_main - rank;

  loc_swap = (dest == rank_pivot) ? i_main : pivot;

  loc_swap /= p;

  loc_matrix_swap = loc_matrix + nm * loc_swap;


  MPI_Sendrecv_replace (loc_matrix_swap, nm, MPI_DOUBLE,
                        dest, 0, dest, 0,
                        comm,  &status);
  return ;
}

void SBC_gauss_MPI_find_pivot (
                               SBC_storage *loc_storage,
                               const int i_main,
                               int *glob_pivot,
                               MPI_Comm comm,
                               MPI_Datatype SBC_pivot_candidate,
                               MPI_Op SBC_pivot_op,
                               double *buf,
                               int *permutation
                              )
{
  pivot_c_struct pivot_candidate;
  pivot_c_struct pivot;

  SBC_gauss_find_candidate_pivot (
                                  loc_storage,
                                  i_main,
                                  &(pivot_candidate.glob_pivot_index),
                                  &(pivot_candidate.norm_of_rev),
                                  buf,
                                  permutation
                                 );

  MPI_Allreduce (&pivot_candidate, &pivot, 1, SBC_pivot_candidate, SBC_pivot_op, MPI_COMM_WORLD);  
  *glob_pivot = pivot.glob_pivot_index;
}

void SBC_gauss_MPI_subtr_all (
                              SBC_storage *loc_storage,
                              const int i_main,
                              MPI_Comm comm,
                              double *buf_,
                              int *permutation
                             )
{
  int rank = loc_storage->rank;
  int n = loc_storage->n;
  int m = loc_storage->m;
  int p = loc_storage->p;
  double *loc_matrix = loc_storage->loc_matrix;
  int nm = n * m;
  int mm = m * m;
  double *coef;
  int count = nm - mm * (i_main + 1);
  int source = i_main % p;
  double *coef_matrix_array = buf_;
  double *buf = buf_ + nm;

  if (source == rank)
    {
      int loc_col_coef = i_main / p;
      double *to_send = loc_matrix + loc_col_coef * nm + (i_main + 1) * mm;
      MPI_Bcast ( to_send, 
                  count, 
                  MPI_DOUBLE, 
                  source, 
                  comm
                );
      coef = to_send;
    }
  else
    {
      MPI_Bcast(
                coef_matrix_array, 
                count, 
                MPI_DOUBLE, 
                source, 
                comm
               );
      coef = coef_matrix_array;
    }


  SBC_gauss_subtr_all (
                       loc_storage,
                       coef,
                       buf,
                       i_main
                      );
}

void SBC_gauss_MPI_multi_row (
                              SBC_storage *loc_storage,
                              const int i_main,
                              MPI_Comm comm,
                              double *buf_,
                              int *permutation
                             )
{
  int rank = loc_storage->rank;
  int n = loc_storage->n;
  int m = loc_storage->m;
  int p = loc_storage->p; int l = loc_storage->l;
  int s = loc_storage->s;
  int has_s_col = loc_storage->has_s_col;
  int loc_size = loc_storage->total_size;
  double *loc_matrix = loc_storage->loc_matrix;
  int nm = n * m;
  int mm = m * m;
  int source = i_main % p;
  double *coef_matrix = buf_;
  double *buf = buf_ += mm;
  double *block;


  if (i_main == l)
    {
      if (has_s_col)
        {
          block = loc_matrix + loc_size - s * s;
          if (matrix_reverse (block, coef_matrix, permutation, s))
            MPI_Abort (comm, NO_SOLUTION);
          if (rank != 0)
            MPI_Send (coef_matrix, s * s, MPI_DOUBLE, 0, 0, comm);
          else
            {
              SBC_gauss_multi_row (
                                   loc_storage,
                                   coef_matrix,
                                   buf,
                                   i_main
                                  );
              return ;

            }
        }
      if (s != 0 && rank == 0)
        {
          MPI_Status st;
          MPI_Recv (coef_matrix, s * s, MPI_DOUBLE, source, 0, comm, &st); 

          SBC_gauss_multi_row (
                               loc_storage,
                               coef_matrix,
                               buf,
                               i_main
                              );

        }
      return ;
    }

  if (rank == source)
    {
      block = loc_matrix + (i_main / p) * nm + i_main * mm;
      if (matrix_reverse (block, coef_matrix, permutation, m))
        MPI_Abort (comm, NO_SOLUTION);
      MPI_Bcast (coef_matrix, mm, MPI_DOUBLE, source, comm);
    }
  else
    MPI_Bcast (coef_matrix, mm, MPI_DOUBLE, source, comm); 
  SBC_gauss_multi_row (
                       loc_storage,
                       coef_matrix,
                       buf,
                       i_main
                      ); 
}

void SBC_gauss_MPI_backward (
                             SBC_storage *loc_storage,
                             const int i_main,
                             MPI_Comm comm,
                             double *buf_
                            )
{
  int rank = loc_storage->rank;
  int n = loc_storage->n;
  int m = loc_storage->m;
  int p = loc_storage->p;
  int l = loc_storage->l;
  int s = loc_storage->s;
  int nm = n * m;
  double *loc_matrix = loc_storage->loc_matrix;
  int source = i_main % p;
  double *coef_matrix_array = buf_;
  double *buf = buf_ + nm;
  double *to_send;
  double *coef;

  int width = (i_main == l) ? s : m;

  to_send = loc_matrix + (i_main / p) * nm;

  if (width == 0)
    return;
  if (rank != 0 && source != rank)
    return;
  // rank == 0 || source == rank
  if (rank == source)
    {
      if (source != 0)
        {
          MPI_Send (to_send, width * i_main * m, MPI_DOUBLE, 0, 0, comm);
          return;
        }
      else 
        coef = to_send;
    }
  else //rank == 0 && source != 0
    {
      MPI_Status st;
      MPI_Recv (coef_matrix_array, width * i_main * m, MPI_DOUBLE, source, 0, comm, &st);
      coef = coef_matrix_array;
    }
  //rank == 0
  SBC_gauss_subtr_backward (loc_storage, i_main, coef, buf);
}

void SBC_gauss_subtr_backward (
                               SBC_storage *loc_storage,
                               const int i_main,
                               double *coef_col,
                               double *buf
                              )
{
  int m = loc_storage->m;
  int l = loc_storage->l;
  int s = loc_storage->s;
  double *b = loc_storage->b;
  double *minuend = b;
  double *subtrahend_ = buf;
  double *coef = coef_col;

  int width = (i_main == l) ? s : m;

  double *subtrahend = b + i_main * m;

  int iter = width * m;

  if (width == 0)
    {
      printf ("in width == 0\n");
      return;
    }

  for (int i = 0; i < i_main; i++)
    {
      matrix_product  (coef, subtrahend, subtrahend_, m, width, 1);
      matrix_subtr_rw (minuend, subtrahend_, m, 1);
      minuend += m;
      coef += iter;
    }
}

void SBC_MPI_multi_vector (
                           SBC_storage *loc_storage,
                           double *x,
                           double *b_save,
                           MPI_Comm comm,
                           double *buf
                          )
{
  int rank = loc_storage->rank;
  int n = loc_storage->n;
  int m = loc_storage->m;
  int p = loc_storage->p;
  int l = loc_storage->l;
  int s = loc_storage->s;
  int number_of_m_cols = loc_storage->number_of_m_cols;
  int has_s_col = loc_storage->has_s_col;
  int mm = m * m;
  int sm = s * m;
  int pm = p * m;
  double *loc_matrix = loc_storage->loc_matrix;
  double *block = loc_matrix;
  double *x_iter = x + rank * m;
  double *b_iter = b_save;

  MPI_Bcast (x, n, MPI_DOUBLE, 0, comm);

  for (int i = 0; i < n; i++)
    b_save[i] = 0;

  for (int i = 0; i < number_of_m_cols; i++)
    {
      b_iter = b_save;
      for (int j = 0; j < l; j++)
        {
          matrix_product (block, x_iter, buf, m, m, 1); 

          for (int k = 0; k < m; k ++)
            b_iter[k] += buf[k]; 

          b_iter += m;
          block += mm;
        } 
      matrix_product (block, x_iter, buf, s, m , 1);

      for (int k = 0; k < s; k ++)
        b_iter[k] += buf[k];

      block += sm; 
      x_iter += pm;
    }

  if (has_s_col)
    {
      b_iter = b_save;
      for (int j = 0; j < l; j++)
        {
          matrix_product (block, x_iter, buf, m, s, 1); 

          for (int k = 0; k < m; k ++)
            b_iter[k] += buf[k]; 

          b_iter += m;
          block += sm;
        } 
      matrix_product (block, x_iter, buf, s, s, 1);

      for (int k = 0; k < s; k ++)
        b_iter[k] += buf[k];
    }

  for (int i = 0; i < n; i++)
    {
      // printf ("[%d]: b_save[%d] = %lf\n", rank, i, b_save[i]);
    }
  MPI_Reduce (b_save, x, n, MPI_DOUBLE, MPI_SUM, 0, comm);
}


void SBC_destroy (
                  SBC_storage *loc_storage
                 )
{
  if (loc_storage->state >= SBC_STORAGE_ALLOCATED) 
    delete[] loc_storage->loc_matrix;
  if (loc_storage->b)
    delete[] loc_storage->b;
  loc_storage->state = SBC_STORAGE_DECLARED;
}

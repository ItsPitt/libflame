/*

    Copyright (C) 2014, The University of Texas at Austin
    Copyright (C) 2022, Advanced Micro Devices, Inc.

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_HIP

#include <flame/hiphopper.h>
#include "hip/hip_runtime_api.h"
#include "rocblas.h"
#include "rocsolver.h"

FLA_Error FLA_Apply_pivots_unb_external_hip( rocblas_handle handle, FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A, void* A_hip )
{

  FLA_Datatype datatype;
  int          n_A, cs_A;
  int          m_p;
  int          inc_p;
  int*         buff_p;
  int          k1_1, k2_1;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Apply_pivots_check( side, trans, p, A );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  /*if ( FLA_Obj_length( A ) > 1 )
  {
    fprintf( stderr,
             "HIP unb pivot apply operation not for objects with length >1. Is: %ld\n",
             FLA_Obj_length( A ) );

    return FLA_FAILURE;
  }*/

  datatype = FLA_Obj_datatype( A );

  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_p    = FLA_Obj_vector_inc( p );
  m_p      = FLA_Obj_vector_dim( p );

  buff_p   = FLA_INT_PTR( p );

  void* A_mat;
  if ( FLASH_Queue_get_malloc_managed_enabled_hip( ) )
  {
    A_mat = FLA_Obj_buffer_at_view( A );
  }
  else
  {
    A_mat = A_hip;
  }

  // Use one-based indices for LAPACK.
  k1_1     = 1;
  k2_1     = m_p;

  // Translate FLAME pivot indices to LAPACK-compatible indices.
  rocblas_int* pivots_lapack;
  hipMalloc( (void**) &pivots_lapack, sizeof( rocblas_int ) * m_p );

  hipError_t err_trans = trans_pivots_FLAME_to_LAPACK( handle, m_p, buff_p, pivots_lapack );
  if ( err_trans != hipSuccess )
  {
    fprintf( stderr,
             "Failure to translate LU pivots to LAPACK convention. err=%d\n",
             err_trans );
    return FLA_FAILURE;
  }

  switch ( datatype ){

  case FLA_FLOAT:
  {
    float* buff_A = ( float * ) A_mat; 

    rocsolver_slaswp( handle,
                      n_A,
                      buff_A, cs_A,
                      k1_1, 
                      k2_1,
                      pivots_lapack,
                      inc_p );
    break;
  }

  case FLA_DOUBLE:
  {
    double* buff_A = ( double * ) A_mat;

    rocsolver_dlaswp( handle,
                      n_A,
                      buff_A, cs_A,
                      k1_1, 
                      k2_1,
                      pivots_lapack,
                      inc_p );
    break;
  }

  case FLA_COMPLEX:
  {
    rocblas_float_complex* buff_A = ( rocblas_float_complex * ) A_mat;

    rocsolver_claswp( handle,
                      n_A,
                      buff_A, cs_A,
                      k1_1, 
                      k2_1,
                      pivots_lapack,
                      inc_p );
    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    rocblas_double_complex* buff_A = ( rocblas_double_complex * ) A_mat;

    rocsolver_zlaswp( handle,
                      n_A,
                      buff_A, cs_A,
                      k1_1, 
                      k2_1,
                      pivots_lapack,
                      inc_p );
    break;
  }

  }

  hipFree( pivots_lapack );

  return FLA_SUCCESS;
}

FLA_Error FLA_Apply_pivots_ln_unb_ext_hip( rocblas_handle handle, FLA_Obj p, FLA_Obj A, void* A_hip )
{
  return FLA_Apply_pivots_unb_external_hip( handle, FLA_LEFT, FLA_NO_TRANSPOSE, p, A, A_hip );
}

#endif

/*

    Copyright (C) 2014, The University of Texas at Austin
    Copyright (C) 2022, Advanced Micro Devices, Inc.

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_HIP

#include "rocblas.h"

FLA_Error FLA_Axpy_external_hip( rocblas_handle handle, FLA_Obj alpha, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip )
{
  FLA_Datatype datatype;
  int          m_B, n_B;
  int          ldim_A, inc_A;
  int          ldim_B, inc_B;
  int          i;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Axpy_check( alpha, A, B );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  ldim_A   = FLA_Obj_col_stride( A );
  inc_A    = 1;

  m_B      = FLA_Obj_length( B );
  n_B      = FLA_Obj_width( B );
  ldim_B   = FLA_Obj_col_stride( B );
  inc_B    = 1;

  switch ( datatype ){

  case FLA_FLOAT:
  {
    float* buff_alpha = ( float* ) FLA_FLOAT_PTR( alpha );
    float* buff_A_hip = ( float* ) A_hip;
    float* buff_B_hip = ( float* ) B_hip;

    for ( i = 0; i < n_B; i++ )
      rocblas_saxpy( handle,
                     m_B,
                     buff_alpha,
                     buff_A_hip + i * ldim_A, inc_A,
                     buff_B_hip + i * ldim_B, inc_B );

    break;
  }

  case FLA_DOUBLE:
  {
    double* buff_alpha = ( double* ) FLA_DOUBLE_PTR( alpha );
    double* buff_A_hip = ( double* ) A_hip;
    double* buff_B_hip = ( double* ) B_hip;

    for ( i = 0; i < n_B; i++ )
      rocblas_daxpy( handle,
                     m_B,
                     buff_alpha,
                     buff_A_hip + i * ldim_A, inc_A,
                     buff_B_hip + i * ldim_B, inc_B );

    break;
  }

  case FLA_COMPLEX:
  {
    rocblas_float_complex* buff_alpha = ( rocblas_float_complex* ) FLA_COMPLEX_PTR( alpha );
    rocblas_float_complex* buff_A_hip = ( rocblas_float_complex* ) A_hip;
    rocblas_float_complex* buff_B_hip = ( rocblas_float_complex* ) B_hip;

    for ( i = 0; i < n_B; i++ )
      rocblas_caxpy( handle,
                     m_B,
                     buff_alpha,
                     buff_A_hip + i * ldim_A, inc_A,
                     buff_B_hip + i * ldim_B, inc_B );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    rocblas_double_complex* buff_alpha = ( rocblas_double_complex* ) FLA_DOUBLE_COMPLEX_PTR( alpha );
    rocblas_double_complex* buff_A_hip = ( rocblas_double_complex* ) A_hip;
    rocblas_double_complex* buff_B_hip = ( rocblas_double_complex* ) B_hip;

    for ( i = 0; i < n_B; i++ )
      rocblas_zaxpy( handle,
                     m_B,
                     buff_alpha,
                     buff_A_hip + i * ldim_A, inc_A,
                     buff_B_hip + i * ldim_B, inc_B );

    break;
  }

  }
  
  return FLA_SUCCESS;
}

#endif

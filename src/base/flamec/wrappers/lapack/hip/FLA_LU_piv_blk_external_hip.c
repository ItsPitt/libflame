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

FLA_Error FLA_LU_piv_blk_external_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj p )
{
  FLA_Datatype datatype;
  dim_t        m_A, n_A, cs_A;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_LU_piv_check( A, p );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  rocblas_int* info;
  hipMalloc( (void**) &info, sizeof( rocblas_int ) );

  void* A_mat;
  void* p_vec;
  if ( FLASH_Queue_get_malloc_managed_enabled_hip( ) )
  {
    A_mat = FLA_Obj_buffer_at_view( A );
    p_vec = FLA_Obj_buffer_at_view( p );
  }
  else
  {
    A_mat = A_hip;
    p_vec = FLA_Obj_buffer_at_view( p ); // it's host malloc'd, so slow but usable
  }

  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A = ( float * ) A_mat;
    int   *buff_p = ( int   * ) p_vec;

    rocsolver_sgetrf( handle,
                      m_A,
                      n_A,
                      buff_A, cs_A,
                      buff_p,
                      info );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A = ( double * ) A_mat;
    int    *buff_p = ( int    * ) p_vec;

    rocsolver_dgetrf( handle,
                      m_A,
                      n_A,
                      buff_A, cs_A,
                      buff_p,
                      info );

    break;
  } 

  case FLA_COMPLEX:
  {
    rocblas_float_complex *buff_A = ( rocblas_float_complex * ) A_mat;
    int      *buff_p = ( int      * ) p_vec;

    rocsolver_cgetrf( handle,
                      m_A,
                      n_A,
                      buff_A, cs_A,
                      buff_p,
                      info );

    break;
  } 

  case FLA_DOUBLE_COMPLEX:
  {
    rocblas_double_complex *buff_A = ( rocblas_double_complex * ) A_mat;
    int      *buff_p = ( int      * ) p_vec;

    rocsolver_zgetrf( handle,
                      m_A,
                      n_A,
                      buff_A, cs_A,
                      buff_p,
                      info );

    break;
  } 

  }

  hipError_t shift_err = shift_pivots_LAPACK_to_FLAME( handle, FLA_Obj_length( p ), p_vec );
  if ( shift_err != hipSuccess )
  {
    fprintf( stderr,
             "Failure to shift LU pivots to FLAME convention. err=%d\n",
             shift_err );
    return FLA_FAILURE;
  }

  hipFree( info );

  return FLA_SUCCESS;
}

#endif

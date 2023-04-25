/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLASH

#include "FLASH_lapack2flash_util_defs.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_prototypes.h"

/*
   GELQF computes an LQ factorization of a M-by-N matrix A: A = L * Q.
*/

#define LAPACK_gelqf(prefix)                                            \
  int F77_ ## prefix ## gelqf( int* m,                                  \
                               int* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, int* lwork, \
                               int* info )

#define LAPACK_gelqf_body(prefix)                                       \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);                \
  int          min_m_n  = min( *m, *n );                                \
  FLA_Error    init_result;                                             \
  dim_t        blocksize = FLASH_get_preferred_blocksize();             \
                                                                        \
  /*This is only for testing, once this works correctly, remove this*/  \
  /*Set to 0 to "fix" one issue. Set to 1 to never use the FLA calls*/  \
  int always_use_flash = 0;                                             \
                                                                        \
  FLA_Init_safe( &init_result );                                        \
  if( *m <= *n || always_use_flash ){                                   \
    FLA_Obj      A_flat, A, t, T;                                       \
    /*I sent this up how I believe FLASH_LQ_UT() should be called*/     \
    /*There seems to be an issue while running the netlibs-test 3.11*/  \
    /*where if m > n, it fails and says */                              \
    /*"libflame: Detected inconsistent object dimensions." */           \
    /*Debug, remove this when fixed*/                                   \
    /*fprintf(stderr,"FLASH: m=%d n=%d ldim_A=%d\n", *m, *n, *ldim_A);*/\
    FLA_Obj_create_without_buffer( datatype, *m, *n, &A_flat );         \
    FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A_flat );               \
                                                                        \
    FLASH_Obj_create_without_buffer( datatype,                          \
                                     min_m_n,                           \
                                     1,                                 \
                                     FLASH_get_depth(),                 \
                                     &blocksize,                        \
                                     &t );                              \
    FLASH_Obj_attach_buffer( buff_t, 1, min_m_n, &t );                  \
    FLASH_Set( FLA_ZERO, t );                                           \
                                                                        \
    FLASH_LQ_UT_create_hier_matrices( A_flat,                           \
                                      FLASH_get_depth(),                \
                                      &blocksize,                       \
                                      &A,                               \
                                      &T );                             \
    FLASH_LQ_UT( A, T );                                                \
    FLA_LQ_UT_recover_tau( T, t );                                      \
    PREFIX2FLAME_INVERT_TAU(prefix,t);                                  \
                                                                        \
    FLA_Obj_free_without_buffer( &A_flat );                             \
    FLASH_Obj_free_without_buffer( &t );                                \
    FLASH_Obj_free( &T );                                               \
    FLASH_Obj_free( &A );                                               \
    /*This should possibly free without buffer, which might be a problem.*/     \
    /*If FLASH_LQ_UT_create_hier_matrices() is making a copy of the buffer,*/   \
    /*then this might not be properly manipulating the buffer data*/            \
    /*I'm unsure of the buffer data of A needs to be copied back to A_flat*/    \
    /*or if the buffer needs to be attached to A after it is made from A_flat*/ \
    /*In one of the netlibs-tests in 3.11, './EIG/xeigtsts < se2.in'. It seems*/\
    /*to loop for a long time until the test dies.*/                            \
    /*gdb shows it seg fault at sdrvst2stg_ , but I suspect that*/              \
    /*it's actually an issue of not updating the buffer correctly*/             \
  }else{                                                                \
    FLA_Obj      A, t, T;                                               \
    /*If the issue where m>n breaking this gets fixed,*/                \
    /*this block that defaults to FLA_ can be removed*/                 \
    /*Debug, remove this when fixed*/                                   \
    /*fprintf(stderr, "FLA: m=%d n=%d ldim_A=%d\n", *m, *n, *ldim_A);*/ \
    FLA_Obj_create_without_buffer( datatype, *m, *n, &A );              \
    FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                    \
                                                                        \
    FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &t );          \
    FLA_Obj_attach_buffer( buff_t, 1, min_m_n, &t );                    \
                                                                        \
    FLA_Set( FLA_ZERO, t );                                             \
                                                                        \
    FLA_LQ_UT_create_T( A, &T );                                        \
    FLA_LQ_UT( A, T );                                                  \
    FLA_LQ_UT_recover_tau( T, t );                                      \
    PREFIX2FLAME_INVERT_TAU(prefix,t);                                  \
                                                                        \
    FLA_Obj_free_without_buffer( &A );                                  \
    FLA_Obj_free_without_buffer( &t );                                  \
    FLA_Obj_free( &T );                                                 \
  }                                                                     \
  FLA_Finalize_safe( init_result );                                     \
                                                                        \
  *info = 0;                                                            \
                                                                        \
  return 0;

LAPACK_gelqf(s)
{
    {
        LAPACK_RETURN_CHECK( sgelqf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_gelqf_body(s)
    }
}
LAPACK_gelqf(d)
{
    {
        LAPACK_RETURN_CHECK( dgelqf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_gelqf_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gelqf(c)
{
    {
        LAPACK_RETURN_CHECK( cgelqf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_gelqf_body(c)
    }
}
LAPACK_gelqf(z)
{
    {
        LAPACK_RETURN_CHECK( zgelqf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_gelqf_body(z)
    }
}
#endif

#define LAPACK_gelq2(prefix)                                            \
  int F77_ ## prefix ## gelq2( int* m,                                  \
                               int* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,   \
                               int* info )
LAPACK_gelq2(s)
{
    {
        LAPACK_RETURN_CHECK( sgelq2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_gelqf_body(s)
    }
}
LAPACK_gelq2(d)
{
    {
        LAPACK_RETURN_CHECK( dgelq2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_gelqf_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gelq2(c)
{
    {
        LAPACK_RETURN_CHECK( cgelq2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_gelqf_body(c)
    }
}
LAPACK_gelq2(z)
{
    {
        LAPACK_RETURN_CHECK( zgelq2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_gelqf_body(z)
    }
}
#endif


#endif

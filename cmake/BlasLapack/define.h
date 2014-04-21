/* Choose the OS  */
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    #define USE_WINDOWS
#else
    #define USE_LINUX
#endif

#ifdef USE_LINUX
    // Misc
    #define dfill  dfill_
    // Level 1 BLAS Routines
    #define dswap  dswap_
    #define dscal  dscal_
    #define dcopy  dcopy_
    #define daxpy  daxpy_
    #define  ddot  ddot_
    #define dasum  dasum_
    #define damax  damax_
    // Level 2 BLAS Routines
    #define dgemv  dgemv_
    #define zgemv  zgemv_
    #define  dger  dger_
    #define zgeru  zgeru_
    // Level 3 BLAS Routines
    #define dgemm  dgemm_
    #define zgemm  zgemm_
    // LAPACK Routines 
    #define  dgesv  dgesv_
    #define  dgtsv  dgtsv_
    #define  dgbmv  dgbmv_
    #define dgetrf  dgetrf_
    #define dgttrf  dgttrf_
    #define dgbtrf  dgbtrf_
    #define dgetrs  dgetrs_
    #define dgttrs  dgttrs_
    #define dgbtrs  dgbtrs_
    #define dgetri  dgetri_
    #define dgbsv   dgbsv_
#endif


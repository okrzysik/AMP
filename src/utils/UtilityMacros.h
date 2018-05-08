// This file contains useful macros including AMP_ERROR, AMP_WARNING, AMP_INSIST, AMP_ASSERT, etc.
#ifndef included_AMP_UtilityMacros
#define included_AMP_UtilityMacros

#include "AMP/utils/Logger.h"
#include "AMP/utils/Utilities.h"


/*! \defgroup Macros Set of utility macro functions used in AMP
 *  \details  These functions are a list of C++ macros that are used within AMP
 *     for common operations, including checking for errors and profiling.
 *  \addtogroup Macros
 *  @{
 */


/*! \def NULL_STATEMENT
 *  \brief    A null statement
 *  \details  A statement that does nothing, for insure++ make it something
 * more complex than a simple C null statement to avoid a warning.
 */
#ifdef __INSURE__
#define NULL_STATEMENT            \
    do {                          \
        if ( 0 )                  \
            int nullstatement = 0 \
    } while ( 0 )
#else
#define NULL_STATEMENT
#endif


/*! \def NULL_USE(variable)
 *  \brief    A null use of a variable
 *  \details  A null use of a variable, use to avoid GNU compiler warnings about unused variables.
 *  \param variable  Variable to pretend to use
 */
#undef NULL_USE
#define NULL_USE( variable )                    \
    do {                                        \
        if ( 0 ) {                              \
            auto static t = (char *) &variable; \
            t++;                                \
        }                                       \
    } while ( 0 )


// Get a ostream
#ifndef LACKS_SSTREAM
#define TBOXOSTREAM std::ostringstream
#else
#define TBOXOSTREAM std::ostrstream
#endif


/*! \def AMP_ERROR(MSG)
 *  \brief      Throw error
 *  \details    Throw an error exception from within any C++ source code.  The
 *     macro argument may be any standard ostream expression.  The file and
 *     line number of the abort are also printed.
 *  \param MSG  Error message to print
 */
#define AMP_ERROR( MSG )                                           \
    do {                                                           \
        TBOXOSTREAM tboxos;                                        \
        tboxos << MSG;                                             \
        AMP::Utilities::abort( tboxos.str(), __FILE__, __LINE__ ); \
    } while ( 0 )


/*! \def AMP_WARNING(MSG)
 *  \brief   Print a warning
 *  \details Print a warning without exit.  Print file and line number of the warning.
 *  \param MSG  Warning message to print
 */
#define AMP_WARNING( MSG )                                                          \
    do {                                                                            \
        TBOXOSTREAM tboxos;                                                         \
        tboxos << MSG << std::ends;                                                 \
        printf( "WARNING: %s\n   Warning called in %s on line %i\n",                \
                tboxos.str().c_str(),                                               \
                __FILE__,                                                           \
                __LINE__ );                                                         \
        AMP::Logger::getInstance()->logWarning( tboxos.str(), __FILE__, __LINE__ ); \
    } while ( 0 )


/*! \def AMP_DEBUG(MSG)
 *  \brief   Print a debug without exit.
 *  \details Print a debug without exit.  Print file and line number of the debug.
 *  \param MSG  Debug message to print
 */
#define AMP_DEBUG( MSG )                                                          \
    do {                                                                          \
        TBOXOSTREAM tboxos;                                                       \
        tboxos << MSG << std::ends;                                               \
        printf( "WARNING: %s\n   Warning called in %s on line %i\n",              \
                tboxos.str().c_str(),                                             \
                __FILE__,                                                         \
                __LINE__ );                                                       \
        AMP::Logger::getInstance()->logDebug( tboxos.str(), __FILE__, __LINE__ ); \
    } while ( 0 )


/*! \def AMP_ASSERT(EXP)
 *  \brief Assert error
 *  \details Throw an error exception from within any C++ source code if the
 *     given expression is not true.  This is a parallel-friendly version
 *     of assert.
 *     The file and line number of the abort are printed along with the stack trace (if availible).
 *  \param EXP  Expression to evaluate
 */
#define AMP_ASSERT( EXP )                                              \
    do {                                                               \
        if ( !( EXP ) ) {                                              \
            TBOXOSTREAM tboxos;                                        \
            tboxos << "Failed assertion: " << #EXP;                    \
            AMP::Utilities::abort( tboxos.str(), __FILE__, __LINE__ ); \
        }                                                              \
    } while ( 0 )


/*! \def AMP_INSIST(EXP,MSG)
 *  \brief Insist error
 *  \details Throw an error exception from within any C++ source code if the
 *     given expression is not true.  This will also print the given message.
 *     This is a parallel-friendly version of assert.
 *     The file and line number of the abort are printed along with the stack trace (if availible).
 *  \param EXP  Expression to evaluate
 *  \param MSG  Debug message to print
 */
#define AMP_INSIST( EXP, MSG )                                         \
    do {                                                               \
        if ( !( EXP ) ) {                                              \
            TBOXOSTREAM tboxos;                                        \
            tboxos << "Failed insist: " << #EXP << std::endl;          \
            tboxos << "Message: " << MSG << std::ends;                 \
            AMP::Utilities::abort( tboxos.str(), __FILE__, __LINE__ ); \
        }                                                              \
    } while ( 0 )


/**
 * Macro for use when assertions are to be included
 * only when debugging.
 */
/*! \def AMP_CHECK_ASSERT(EXP)
 *  \brief Assert error (debug only)
 *  \details Throw an error exception from within any C++ source code if the
 *     given expression is not true.  This only runs if DEBUG_CHECK_ASSERTIONS
 *     is enabled.  If enabled, this is the same as a call to AMP_ASSERT.
 *  \param EXP  Expression to evaluate
 */
#ifdef DEBUG_CHECK_ASSERTIONS
#define AMP_CHECK_ASSERT( EXP ) AMP_ASSERT( EXP )
#else
#define AMP_CHECK_ASSERT( EXP )
#endif


/*! \def TYPE_HASH(X)
 *  \brief Get a hash key from the class type
 *  \details Get a hash key from the class type.
 *      This requires the RTTI (Run-time type information) to be available.
 *  \param X  Class to use for the hash key
 */
#define TYPE_HASH( X ) AMP::Utilities::hash_char( typeid( X ).name() )


/*! \def __VA_NARG__(...)
 *  \brief Macros to return the number of arguments in __VA_ARGS__
 *  \details These macros will return the number of arguments in __VA_ARGS__
 */
// clang-format off
#define __VA_NARG__( ... ) ( __VA_NARG_( _0, ##__VA_ARGS__, __RSEQ_N() ) - 1 )
#define __VA_NARG_( ... ) __VA_ARG_N( __VA_ARGS__ )
#define __VA_ARG_N( _1,  _2,  _3,  _4,  _5,  _6,  _7,  _8,  _9,  _10,  _11,  _12,  _13,  _14,   \
    _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26, _27, _28, _29, _30, _31, _32,   \
    _33, _34, _35, _36, _37, _38, _39, _40, _41, _42, _43, _44, _45, _46, _47, _48, _49, _50,   \
    _51, _52, _53, _54, _55, _56, _57, _58, _59, _60, _61, _62, _63, N,  ... )  N
#define __RSEQ_N()                                                                              \
    63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, \
        40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, \
        18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
// clang-format on


/**
 * Throw an error exception from within any C++ source code.  This is
 * is similar to AMP_ERROR(), but is designed to be invoked after a
 * call to a PETSc library function.  In other words, it acts similarly
 * to the PETSc CHKERRQ(ierr) macro.
 */
#ifdef HAVE_PETSC
/*
 * In the following, "CHKERRCONTINUE(ierr);" will cause PETSc to print out
 * a stack trace that led to the error; this may be useful for debugging.
 */
#ifndef LACKS_SSTREAM
#define PETSC_AMP_ERROR( ierr )                                        \
    do {                                                               \
        if ( ierr ) {                                                  \
            std::ostringstream tboxos;                                 \
            AMP::Utilities::abort( tboxos.str(), __FILE__, __LINE__ ); \
        }                                                              \
    }                                                                  \
    }                                                                  \
    while ( 0 )
#else
#define PETSC_AMP_ERROR( ierr )                                        \
    do {                                                               \
        if ( ierr ) {                                                  \
            std::ostrstream tboxos;                                    \
            CHKERRCONTINUE( ierr );                                    \
            AMP::Utilities::abort( tboxos.str(), __FILE__, __LINE__ ); \
        }                                                              \
    }                                                                  \
    }                                                                  \
    while ( 0 )
#endif
#endif


/*! \def DISABLE_WARNINGS
 *  \brief Reenable warnings
 *  \details This will re-enable warnings after a call to DIASABLE_WARNINGS
 */
/*! \def ENABLE_WARNINGS
 *  \brief Supress all warnings
 *  \details This will start to supress all compile warnings.
 *      Be sure to follow with ENABLE_WARNINGS
 */
// clang-format off
#ifndef DISABLE_WARNINGS
#if defined( USING_MSVC )
    #define DISABLE_WARNINGS __pragma( warning( push, 0 ) )
    #define ENABLE_WARNINGS __pragma( warning( pop ) )
#elif defined( USING_CLANG )
    #define DISABLE_WARNINGS                                                            \
        _Pragma( "clang diagnostic push" ) _Pragma( "clang diagnostic ignored \"-Wall\"" ) \
        _Pragma( "clang diagnostic ignored \"-Wextra\"" )                               \
        _Pragma( "clang diagnostic ignored \"-Wunused-private-field\"" )                \
        _Pragma( "clang diagnostic ignored \"-Wdeprecated-declarations\"" )             \
        _Pragma( "clang diagnostic ignored \"-Winteger-overflow\"" )
    #define ENABLE_WARNINGS _Pragma( "clang diagnostic pop" )
#elif defined( USING_GCC )
    #define DISABLE_WARNINGS                                                            \
        _Pragma( "GCC diagnostic push" ) _Pragma( "GCC diagnostic ignored \"-Wall\"" )  \
        _Pragma( "GCC diagnostic ignored \"-Wextra\"" )                                 \
        _Pragma( "GCC diagnostic ignored \"-Wunused-local-typedefs\"" )                 \
        _Pragma( "GCC diagnostic ignored \"-Woverloaded-virtual\"" )                    \
        _Pragma( "GCC diagnostic ignored \"-Wunused-parameter\"" )                      \
        _Pragma( "GCC diagnostic ignored \"-Wdeprecated-declarations\"" )               \
        _Pragma( "GCC diagnostic ignored \"-Wvirtual-move-assign\"" )                   \
        _Pragma( "GCC diagnostic ignored \"-Wunused-function\"" )                       \
        _Pragma( "GCC diagnostic ignored \"-Woverflow\"" )                              \
        _Pragma( "GCC diagnostic ignored \"-Wunused-variable\"" )                       \
        _Pragma( "GCC diagnostic ignored \"-Wignored-qualifiers\"" )                    \
        _Pragma( "GCC diagnostic ignored \"-Wenum-compare\"" )                          \
        _Pragma( "GCC diagnostic ignored \"-Wterminate\"" )                             \
        _Pragma( "GCC diagnostic ignored \"-Wpragmas\"" )
#define ENABLE_WARNINGS _Pragma( "GCC diagnostic pop" )
#else
#define DISABLE_WARNINGS
#define ENABLE_WARNINGS
#endif
#endif
// clang-format on

/*! @} */

#endif

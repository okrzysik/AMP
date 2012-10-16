// This file contains useful macros including AMP_ERROR, AMP_WARNING, AMP_INSIST, AMP_ASSERT, etc.
#ifndef included_AMP_UtilityMacros
#define included_AMP_UtilityMacros

#include "Utilities.h"


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
    #define NULL_STATEMENT ({ if(0) int nullstatement=0 })
#else
    #define NULL_STATEMENT
#endif


/*! \def NULL_USE(variable)
 *  \brief    A null use of a variable
 *  \details  A null use of a variable, use to avoid GNU compiler warnings about unused variables.
 *  \param variable  Variable to pretend to use
 */
#define NULL_USE(variable) ({ \
    if(0) {char *temp = (char *)&variable; temp++;} \
})


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
#define AMP_ERROR(MSG) ({					                    \
    TBOXOSTREAM tboxos;					                        \
    tboxos << MSG << std::ends;					                \
    AMP::Utilities::abort(tboxos.str(), __FILE__, __LINE__);    \
})


/*! \def AMP_WARNING(MSG)
 *  \brief   Print a warning
 *  \details Print a warning without exit.  Print file and line number of the warning.
 *  \param MSG  Warning message to print
 */
#define AMP_WARNING(MSG) ({					                    \
    TBOXOSTREAM tboxos;					                        \
    tboxos << MSG << std::ends;					                \
    AMP::Logger::getInstance() -> logWarning(tboxos.str(), __FILE__, __LINE__);\
})


/*! \def AMP_DEBUG(MSG)
 *  \brief   Print a debug without exit.
 *  \details Print a debug without exit.  Print file and line number of the debug.
 *  \param MSG  Debug message to print
 */
#define AMP_DEBUG(MSG) ({					                    \
    TBOXOSTREAM tboxos;					                        \
    tboxos << MSG << std::ends;					                \
    AMP::Logger::getInstance() -> logDebug(tboxos.str(), __FILE__, __LINE__);\
})


/*! \def AMP_ASSERT(EXP)
 *  \brief Assert error
 *  \details Throw an error exception from within any C++ source code if the
 *     given expression is not true.  This is a parallel-friendly version
 *     of assert.
 *     The file and line number of the abort are printed along with the stack trace (if availible).
 *  \param EXP  Expression to evaluate
 */
#define AMP_ASSERT(EXP) ({                                          \
    if ( !(EXP) ) {                                                 \
        TBOXOSTREAM tboxos;					                        \
        tboxos << "Failed assertion: " << #EXP << std::ends;        \
        AMP::Utilities::abort(tboxos.str(), __FILE__, __LINE__);    \
    }                                                               \
})


/*! \def AMP_INSIST(EXP,MSG)
 *  \brief Insist error
 *  \details Throw an error exception from within any C++ source code if the
 *     given expression is not true.  This will also print the given message.
 *     This is a parallel-friendly version of assert.
 *     The file and line number of the abort are printed along with the stack trace (if availible).
 *  \param EXP  Expression to evaluate
 *  \param MSG  Debug message to print
 */
#define AMP_INSIST(EXP,MSG) ({                                    \
    if ( !(EXP) ) {                                                 \
        TBOXOSTREAM tboxos;					                        \
        tboxos << "Failed insist: " << #EXP << std::endl;           \
        tboxos << "Message: " << MSG << std::ends;                  \
        AMP::Utilities::abort(tboxos.str(), __FILE__, __LINE__);    \
    }                                                               \
})



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
    #define AMP_CHECK_ASSERT(EXP) AMP_ASSERT(EXP)
#else
    #define AMP_CHECK_ASSERT(EXP) 
#endif



/*! \def TYPE_HASH(X)
 *  \brief Get a hash key from the class type
 *  \details Get a hash key from the class type.  
 *      This requires the RTTI (Run-time type information) to be available.
 *  \param X  Class to use for the hash key
 */
#define TYPE_HASH(X)  AMP::Utilities::hash_char(typeid(X).name())




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
    #define PETSC_AMP_ERROR(ierr) ({						\
          if (ierr) {                                   				\
             std::ostringstream tboxos;							\
             AMP::Utilities::abort(tboxos.str(), __FILE__, __LINE__);	\
          } 									\
    })
    #else
    #define PETSC_AMP_ERROR(ierr) ({						\
          if (ierr) {                                   				\
             std::ostrstream tboxos;							\
             CHKERRCONTINUE(ierr); 							\
             AMP::Utilities::abort(tboxos.str(), __FILE__, __LINE__);	        \
          } 									\
    })
    #endif
#endif



/*! @} */


#endif

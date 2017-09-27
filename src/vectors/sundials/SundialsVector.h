#ifndef included_AMP_SundialsVector
#define included_AMP_SundialsVector


#include "vectors/Vector.h"

extern "C" {
#include "sundials/sundials_nvector.h"
}

namespace AMP {
namespace LinearAlgebra {

/**  \class  SundialsVectorParameters
 *  \brief  Parameters class to construct a Sundials N_Vector
 *
 *  Since the constructor of SundialsVector is protected, this class is of little
 *  use to any but Backplane developers.  Also, since SundialsVector
 *  is almost trivial, so is the parameters class.
 */
class SundialsVectorParameters
{
};

/**
 *  \class  SundialsVector
 *  \brief  SundialsVector is a bridge between AMP::LinearAlgebra::Vector and
 *  the Sundials N_Vector data structure.
 *
 *  A SundialsVector has a N_Vector data structure.  Given an
 *  AMP::LinearAlgebra::Vector, this class can create a Sundials view without
 *  copying the data.  As such, this class serves three purposes:
 *  -# Provides a Sundials N_Vector for derived classes to use, fill, manage, etc.
 *  -# Provides an interface for accessing this Sundials N_Vector independent of derived classes
 *  -# Provides a static method for creating a Sundials view of an AMP Vector.
 */
class SundialsVector
{
protected:
    /**
     *  \brief  Sundials NVector wrapping the data in the Vector
     *
     *  However this is created, the N_Vector holds a pointer to the data used in the Vector
     *  in such a manner that is consistent with Sundials.
     */
    N_Vector d_n_vector;

    /**
     *  \brief  Construct a SundialsVector
     *
     *  This can only be called by a derived class or the static function below.  There is
     *  no need to create this vector directly since it is virtual.
     */
    SundialsVector();

public:
    /**
      *  \brief  Obtain a Sundials N_Vector for use in Sundials routines
      *
      *  This function is used to get a Sundials vector.  The following idiom should
      *  be used since it fails gracefully.  In this function, a view
      *  may be created before the NVector is extracted.
      \code
      double DoSundialsMin( Vector::shared_ptr &in )
      {
        // Create a Sundials N_Vector if necessary
        Vector::shared_ptr  in_sundials_view = SundialsVector::view( in );
        // Extract the N_Vector
        N_Vector in_nvector = dynamic_pointer_cast<SundialsVector>(in_sundials_view)->getNVector();
        return N_VMin( in_nvector );
      }
      \endcode
      */
    N_Vector &getNVector();

    /**
      *  \brief  Obtain a Sundials N_Vector for use in Sundials routines
      *
      *  This function is used to get a Sundials vector.  The following idiom should
      *  be used since it fails gracefully.  In this function, a view
      *  may be created before the NVector is extracted.
      \code
      double DoSundialsMin( Vector::shared_ptr &in )
      {
        // Create a Sundials N_Vector if necessary
        Vector::shared_ptr  in_sundials_view = SundialsVector::view( in );
        // Extract the N_Vector
        N_Vector in_nvector = dynamic_pointer_cast<SundialsVector>(in_sundials_view)->getNVector();
        return N_VMin( in_nvector );
      }
      \endcode
      */
    const N_Vector &getNVector() const;

    /**
     *  \param  AmpVector  a shared pointer to a Vector
     *  \brief  If one does not exist, create a Sundials wrapper for AmpVector.  Otherwise, return
     * AmpVector.
     *
     *  The function attempts to return a view with the least amount of work.
     *  IT WILL NEVER COPY DATA.
     *  - If AmpVector is already a SundialsVector, it is returned.
     *  - Else, if AmpVector is a ManagedVector, it is wrapped and returned
     *  - Else, if AmpVector can be used as a VectorEngine, a new ManagedVector
     *    is created and returned
     *  Otherwise, this function will throw an error.
     */
    static Vector::shared_ptr view( Vector::shared_ptr AmpVector );

    /**
     *  \param  AmpVector  a shared pointer to a Vector
     *  \brief  If one does not exist, create a const Sundials wrapper for AmpVector.  Otherwise,
     * return AmpVector.
     *
     *  The function attempts to return a view with the least amount of work.
     *  IT WILL NEVER COPY DATA.
     *  - If AmpVector is already a SundialsVector, it is returned.
     *  - Else, if AmpVector is a ManagedVector, it is wrapped and returned
     *  - Else, if AmpVector can be used as a VectorEngine, a new ManagedVector
     *    is created and returned.
     *  Otherwise, this function will throw an error.
     */
    static Vector::const_shared_ptr constView( Vector::const_shared_ptr AmpVector );

public:

    inline N_Vector& getNativeVec() { return getNVector(); }
    inline const N_Vector& getNativeVec() const { return getNVector(); }

};


} // namespace LinearAlgebra
} // namespace AMP

#include "SundialsVector.inline.h"


#endif

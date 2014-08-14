#ifndef included_AMP_SubsetVector
#define included_AMP_SubsetVector

#include "Vector.h"
#include <vector>

namespace AMP {
namespace LinearAlgebra {

/** \class SubsetVector
  * \brief This vector is a subset of an AMP Vector
  * \details
    Given an AMP Vector, this class will create a view of a subset of the
    vector.  For instance, if \f$\mathbf{a} = \{ a_0 a_1 a_2 \ldots a_n\}\f$,
    and \f$S\f$ is a set of non-negative integers \f$( s_0 s_1 \ldots s_m )\f$
    such that \f$0 \le s_i \le n\f$, then the SubsetVector will be the
    vector \f$\mathbf{a}_S = \{ a_{s_0} a_{s_1} a_{s_2} \ldots a_{s_m}\} \f$.

    This class provides a factory method called view:
    \code
      AMP::LinearAlgebra::Vector::shared_ptr  vec1;
      AMP::LinearAlgebra::Vector::shared_pt   vec2 = AMP::SubsetVector ( vec1 , subsetVar );
      AMP::LinearAlgebra::Vector::shared_ptr  vec3 = vec2->clone ( "subset2" );
    \code

    Since this is a view, any change to vec2 will be reflected on vec1 and
    vice versa.  vec2 is a sparse vector, with a mapping of new index to old.
    vec3, on the other hand, is a dense vector without an index.  If a lot
    of computation is necessary on the sparse vector, the data can be copied
    in and out:

    \code
      // Subset the vector to make a sparse vector.
      vec2 = AMP::SubsetVector ( vec1 , subsetVar )

      // Copy the sparse vector data to a dense vector
      vec3.copyVector( vec2 );

      // Perform whatever math
      performComputation ( vec3 );

      // Copy data back to vec2, and, consequently, vec1
      vec2.copyVector( vec3 );
    \endcode
  */
class  SubsetVector : public Vector
{

public:
    static Vector::shared_ptr   view ( Vector::shared_ptr , Variable::shared_ptr );
    static Vector::const_shared_ptr   view ( Vector::const_shared_ptr , Variable::shared_ptr );

    virtual std::string  type() const;

    using Vector::cloneVector;
    virtual Vector::shared_ptr  cloneVector ( Variable::shared_ptr ) const;
    virtual size_t numberOfDataBlocks () const;
    virtual size_t sizeOfDataBlock ( size_t i ) const;
    virtual void   swapVectors ( Vector &rhs );
    virtual void   aliasVector ( Vector &rhs );
    virtual size_t getLocalSize () const;
    virtual size_t getGlobalSize () const;
    virtual void   assemble () {}

    virtual void   addValuesByLocalID ( int , size_t * , const double * );
    virtual void   setValuesByLocalID ( int , size_t * , const double * );
    virtual void   getValuesByLocalID ( int , size_t * , double *vals ) const;
    virtual void   addLocalValuesByGlobalID ( int , size_t * , const double * );
    virtual void   setLocalValuesByGlobalID ( int , size_t * , const double * );
    virtual void   getLocalValuesByGlobalID ( int , size_t * , double * ) const ;
    virtual void   putRawData ( const double *in );
    virtual void   copyOutRawData ( double *out ) const;

private:
    SubsetVector  () {}
    void  computeIDMap ();

    void *getRawDataBlockAsVoid ( size_t i );
    const void *getRawDataBlockAsVoid ( size_t i ) const ;

    // Internal data
    Vector::shared_ptr      d_ViewVector;                   // Vector we subsetted for the view
    std::vector<size_t>     d_SubsetLocalIDToViewGlobalID;  // The list of global ID in the parent vector
    std::vector<size_t>     d_dataBlockSize;                // The size of the data blocks
    std::vector<double*>    d_dataBlockPtr;                 // The pointers to the data blocks

};


}
}


#endif


#ifndef included_AMP_DTK_AMPVectorHelpers
#define included_AMP_DTK_AMPVectorHelpers

#include <vectors/Vector.h>

#include <Tpetra_Vector.hpp>
#include <Teuchos_RCP.hpp>

namespace AMP {
namespace Operator {


/**
  * AMP Vector helpers for DTK.
*/
class DTKAMPVectorHelpers
{
public :

    /**
     * Constructor.
     */
    DTKAMPVectorHelpers() 
    { /* ... */ }

    //! Destructor
    ~DTKAMPVectorHelpers() 
    { /* ... */ }

    /*!
     * \brief Something
     */
    static Teuchos::RCP<Tpetra::Vector<double,int,std::size_t> >
    pullTpetraVectorFromAMPVector( const AMP::shared_ptr<const AMP::LinearAlgebra::Vector>& ampVector );

    /*!
     * \brief Something
     */
    static void
    pushTpetraVectorToAMPVector( const Tpetra::Vector<double,int,std::size_t>& tpetraVector,
                                 const AMP::shared_ptr<AMP::LinearAlgebra::Vector>& ampVector );
                                      

};


}
}

#endif

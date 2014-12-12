
#ifndef included_AMP_DTK_AMPVectorHelpers
#define included_AMP_DTK_AMPVectorHelpers

#include <vectors/Vector.h>

#include <Tpetra_Vector.hpp>
#include <Tpetra_Map.hpp>
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
     * \brief Given an AMP vector, create a compatible Tpetra vector and copy
     * the data from the AMP vector into the Tpetra vector.
     */
    static Teuchos::RCP<Tpetra::Vector<double,int,std::size_t> >
    pullTpetraVectorFromAMPVector( 
	const AMP::shared_ptr<const AMP::LinearAlgebra::Vector>& ampVector,
	Teuchos::RCP<const Tpetra::Map<int,std::size_t> > tpetra_map = Teuchos::null );

    /*!
     * \brief Given a Tpetra vector and a compatible AMP vector, copy the data
     * from the Tpetra vector into the AMP vector.
     */
    static void
    pushTpetraVectorToAMPVector( 
	const Tpetra::Vector<double,int,std::size_t>& tpetraVector,
	const AMP::shared_ptr<AMP::LinearAlgebra::Vector>& ampVector );
                                      
    /*!
     * \brief Given a DOFManager, create a compatible TpetraMap.
     */
    static Teuchos::RCP<const Tpetra::Map<int,std::size_t> >
    createTpetraMapFromAMPDOFManager(
	const AMP::shared_ptr<AMP::Discretization::DOFManager> dof_manager );
	
};


}
}

#endif

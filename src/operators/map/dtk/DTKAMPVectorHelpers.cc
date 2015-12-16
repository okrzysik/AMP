
#include "DTKAMPVectorHelpers.h"

namespace AMP {
namespace Operator {

//---------------------------------------------------------------------------//
Teuchos::RCP<Tpetra::Vector<double, int, std::size_t>>
DTKAMPVectorHelpers::pullTpetraVectorFromAMPVector(
    const AMP::shared_ptr<const AMP::LinearAlgebra::Vector> &ampVector,
    Teuchos::RCP<const Tpetra::Map<int, std::size_t>>
        tpetra_map )
{
    // Create a Tpetra map.
    if ( Teuchos::is_null( tpetra_map ) ) {
        AMP::shared_ptr<AMP::Discretization::DOFManager> dofManager = ampVector->getDOFManager();
        tpetra_map = createTpetraMapFromAMPDOFManager( dofManager );
    }

    // Create a Tpetra vector.
    Teuchos::RCP<Tpetra::Vector<double, int, std::size_t>> tpetraVector =
        Tpetra::createVector<double, int, std::size_t>( tpetra_map );

    // Copy the data from amp to tpetra
    ampVector->copyOutRawData( tpetraVector->get1dViewNonConst().getRawPtr() );
    return tpetraVector;
}

//---------------------------------------------------------------------------//
void DTKAMPVectorHelpers::pushTpetraVectorToAMPVector(
    const Tpetra::Vector<double, int, std::size_t> &tpetraVector,
    const AMP::shared_ptr<AMP::LinearAlgebra::Vector> &ampVector )
{
    ampVector->putRawData( tpetraVector.get1dView().getRawPtr() );
}

//---------------------------------------------------------------------------//
Teuchos::RCP<const Tpetra::Map<int, std::size_t>>
DTKAMPVectorHelpers::createTpetraMapFromAMPDOFManager(
    const AMP::shared_ptr<AMP::Discretization::DOFManager> dof_manager )
{
    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::rcp( new Teuchos::MpiComm<int>( dof_manager->getComm().getCommunicator() ) );

    // Iterate over the DOF manager elements and extract their DOF ids.
    AMP::Mesh::MeshIterator meshIterator = dof_manager->getIterator();
    const std::size_t n                  = meshIterator.size();
    Teuchos::Array<std::size_t> entityIDs( n );
    std::vector<std::size_t> dofIndices;
    meshIterator = meshIterator.begin();
    for ( std::size_t i = 0; i < n; ++i, ++meshIterator ) {
        dof_manager->getDOFs( meshIterator->globalID(), dofIndices );
        AMP_ASSERT( dofIndices.size() == 1 );
        entityIDs[i] = dofIndices[0];
    } // end for i
    AMP_ASSERT( meshIterator == meshIterator.end() );

    // Create a map from the DOF ids.
    return Tpetra::createNonContigMap<int, std::size_t>( entityIDs, comm );
}

//---------------------------------------------------------------------------//

} // end namespace Operator
} // end namespace AMP

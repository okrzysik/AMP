
#include "DTKAMPVectorHelpers.h"

#include <Tpetra_Map.hpp>

namespace AMP {
namespace Operator {


Teuchos::RCP<Tpetra::Vector<double,int,std::size_t> >
DTKAMPVectorHelpers::
pullTpetraVectorFromAMPVector( const AMP::shared_ptr<const AMP::LinearAlgebra::Vector>& ampVector )
{

    // Create a Tpetra map.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::rcp( new Teuchos::MpiComm<int>(ampVector->getComm().getCommunicator()) );

    AMP::shared_ptr<AMP::Discretization::DOFManager> dofManager = ampVector->getDOFManager();
    AMP::Mesh::MeshIterator meshIterator = dofManager->getIterator();
    const std::size_t n = meshIterator.size();
    Teuchos::Array<std::size_t> entityIDs(n);
    std::vector<std::size_t> dofIndices;
    meshIterator = meshIterator.begin();
    for (std::size_t i = 0; i < n; ++i, ++meshIterator) {
        dofManager->getDOFs(meshIterator->globalID(), dofIndices);
        AMP_ASSERT( dofIndices.size() == 1 );
        entityIDs[i] = dofIndices[0];
    } // end for i
    AMP_ASSERT( meshIterator == meshIterator.end() );
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > map =
        Tpetra::createNonContigMap<int,std::size_t>( entityIDs, comm );

    // Create a Tpetra vector.
    Teuchos::RCP<Tpetra::Vector<double,int,std::size_t> > tpetraVector =
        Tpetra::createVector<double,int,std::size_t>( map );

    // Copy the data from amp to tpetra
    ampVector->copyOutRawData(tpetraVector->get1dViewNonConst().getRawPtr());
    return tpetraVector;
}

void
DTKAMPVectorHelpers::
pushTpetraVectorToAMPVector( const Tpetra::Vector<double,int,std::size_t>& tpetraVector,
                             const AMP::shared_ptr<AMP::LinearAlgebra::Vector>& ampVector )
{
    ampVector->putRawData(tpetraVector.get1dView().getRawPtr());
}

}
}

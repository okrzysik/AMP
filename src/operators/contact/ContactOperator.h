
#ifndef included_AMP_ContactOperator
#define included_AMP_ContactOperator

#include "AMP/operators/ConstraintsEliminationOperator.h"
#include "AMP/operators/contact/ContactOperatorParameters.h"

namespace AMP {
namespace Operator {

/**
  An abstract base class for representing a linear operator. This class
  stores the matrix representation of the linear operator. It provides
  an implementation of the apply() function.
  @see Operator
  */
class ContactOperator : public ConstraintsEliminationOperator
{

public:
    /**
      Constructor. This resets the matrix shared pointer.
      @param [in] params
      */
    ContactOperator( const std::shared_ptr<ContactOperatorParameters> &params )
        : ConstraintsEliminationOperator( params )
    {
        d_Mesh = ( params->d_Mesh );

        d_GlobalComm  = ( params->d_GlobalComm );
        d_DOFsPerNode = ( params->d_DOFsPerNode );
        d_DOFManager  = ( params->d_DOFManager );

        d_MasterMeshID     = ( params->d_MasterMeshID );
        d_SlaveMeshID      = ( params->d_SlaveMeshID );
        d_MasterBoundaryID = ( params->d_MasterBoundaryID );
        d_SlaveBoundaryID  = ( params->d_SlaveBoundaryID );

        d_MasterMechanicsMaterialModel = ( params->d_MasterMechanicsMaterialModel );
    }

    AMP::Mesh::MeshID getMasterMeshID() const { return d_MasterMeshID; }

    AMP::Mesh::MeshID getSlaveMeshID() const { return d_SlaveMeshID; }

    /**
      @return The local number of constrained DOFs.
      */
    size_t numLocalConstraints() { return d_SlaveIndices.size(); }

    /**
      @return The global number of constrained DOFs.
      */
    size_t numGlobalConstraints() { return d_GlobalComm.sumReduce( d_SlaveIndices.size() ); }

    void getActiveSet( std::vector<AMP::Mesh::MeshElementID> const *&activeSet ) const
    {
        activeSet = &d_ActiveSet;
    }

    std::vector<AMP::Mesh::MeshElementID> const &getActiveSet() const { return d_ActiveSet; }

    virtual void initialize() = 0;

    virtual size_t updateActiveSet( AMP::LinearAlgebra::Vector::shared_ptr displacementFieldVector,
                                    bool skipDisplaceMesh ) = 0;

protected:
    AMP::AMP_MPI d_GlobalComm;
    AMP::Discretization::DOFManager::shared_ptr d_DOFManager;
    size_t d_DOFsPerNode;

    AMP::Mesh::Mesh::shared_ptr d_Mesh;

    AMP::Mesh::MeshID d_MasterMeshID;
    AMP::Mesh::MeshID d_SlaveMeshID;

    int d_MasterBoundaryID;
    int d_SlaveBoundaryID;

    std::vector<AMP::Mesh::MeshElementID> d_InactiveSet;
    std::vector<AMP::Mesh::MeshElementID> d_ActiveSet;

    std::shared_ptr<AMP::Operator::MechanicsMaterialModel> d_MasterMechanicsMaterialModel;

private:
};
} // namespace Operator
} // namespace AMP

#endif

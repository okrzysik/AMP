
#ifndef included_AMP_ContactOperatorParameters
#define included_AMP_ContactOperatorParameters

#include <discretization/DOF_Manager.h>
#include <operators/OperatorParameters.h>
#include <operators/mechanics/MechanicsMaterialModel.h>

namespace AMP {
namespace Operator {

class ContactOperatorParameters : public OperatorParameters
{
public:
    ContactOperatorParameters( const AMP::shared_ptr<AMP::Database> &db ) : OperatorParameters( db )
    {
        AMP_INSIST( d_db->keyExists( "MasterBoundaryID" ), "key not found" );
        d_MasterBoundaryID = d_db->getInteger( "MasterBoundaryID" );
        AMP_INSIST( d_db->keyExists( "SlaveBoundaryID" ), "key not found" );
        d_SlaveBoundaryID = d_db->getInteger( "SlaveBoundaryID" );
    }

    void reset() // because d_Mesh may point to NULL
    {
        std::vector<AMP::Mesh::MeshID> meshIDs = d_Mesh->getBaseMeshIDs();
        AMP_INSIST( d_db->keyExists( "MasterMeshIndex" ), "key not found" );
        d_MasterMeshID = meshIDs[d_db->getInteger( "MasterMeshIndex" )];
        AMP_INSIST( d_db->keyExists( "SlaveMeshIndex" ), "key not found" );
        d_SlaveMeshID = meshIDs[d_db->getInteger( "SlaveMeshIndex" )];
    }

    virtual ~ContactOperatorParameters() {}

    AMP::AMP_MPI d_GlobalComm;

    size_t d_DOFsPerNode;
    AMP::Discretization::DOFManager::shared_ptr d_DOFManager;

    AMP::Mesh::MeshID d_MasterMeshID;
    AMP::Mesh::MeshID d_SlaveMeshID;

    int d_MasterBoundaryID;
    int d_SlaveBoundaryID;

    AMP::shared_ptr<AMP::Operator::MechanicsMaterialModel> d_MasterMechanicsMaterialModel;

protected:
private:
};
} // namespace Operator
} // namespace AMP

#endif

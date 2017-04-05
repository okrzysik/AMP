//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MoabMapOperator.cc
 * \author Steven Hamilton
 * \brief  Member definitions for MoabMapOperator
 */
//---------------------------------------------------------------------------//

#include "MoabMapOperator.h"
#include "cell_hex8.h"
#include "discretization/simpleDOF_Manager.h"
#include "elem.h"

namespace AMP {
namespace Operator {

//---------------------------------------------------------------------------//
/*!
 *\brief Constructor
 */
//---------------------------------------------------------------------------//
MoabMapOperator::MoabMapOperator( const SP_MoabMapParams &params )
    : Base( params ), d_params( params ), d_moab( params->d_moabOp )
{
    // Get parameters from DB
    d_mapVar = params->d_db->getString( "MoabMapVariable" );

    // Get mesh manager
    d_meshMgr = params->d_mesh;

    // Get Moab Interface
    d_moabInterface = d_moab->getMoabInterface();

    // Interpolate to nodes or Gauss points?
    if ( params->d_db->getString( "InterpolateToType" ).compare( "Vertex" ) == 0 ) {
        d_interpType = NODES;
    } else if ( params->d_db->getString( "InterpolateToType" ).compare( "GaussPoint" ) == 0 ) {
        d_interpType = GAUSS_POINTS;
    } else
        AMP_ERROR( "InterpolateToType must be Vertex or GaussPoint" );
}

//---------------------------------------------------------------------------//
/*!
 *\brief Map variable on Moab mesh onto GPs of AMP mesh
 *
 * The apply function takes a variable distribution from a previously executed
 * Moab-based calculation (which is stored in a MoabBasedOperator object)
 * and maps it onto the Gauss points of a mesh adapter.  The vectors f and r
 * will be untouched (and can simply be NULL), the vector r will be populated
 * with the variable values, and the constants a and b will not be used.  A mesh
 * adapter must be set in the MoabMapOperatorParameters object prior to
 * calling apply.
 */
//---------------------------------------------------------------------------//
void MoabMapOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr f,
                             AMP::LinearAlgebra::Vector::const_shared_ptr u,
                             AMP::LinearAlgebra::Vector::shared_ptr r,
                             double a,
                             double b )
{
    AMP_INSIST( r, "Vector r must not be null" );
    AMP_INSIST( d_meshMgr, "Must have Mesh Adapter" );
    AMP_INSIST( d_moab, "Must have Moab Operator" );

    // Build Moab Coupler
    buildMoabCoupler();

    // Create vector to hold coordinates
    std::vector<double> allCoords;

    // Loop over meshes
    std::vector<AMP::Mesh::MeshID> meshIDs = d_meshMgr->getBaseMeshIDs();
    for ( size_t meshIndex = 0; meshIndex < meshIDs.size(); meshIndex++ ) {
        // this is an accessor to all the mesh info.
        AMP::Mesh::Mesh::shared_ptr currentMesh = d_meshMgr->Subset( meshIDs[meshIndex] );
        if ( currentMesh.get() == NULL )
            continue;

        std::vector<double> theseCoords;
        switch ( d_interpType ) {
        case NODES: {
            // Get nodes coords
            getNodeCoords( currentMesh, theseCoords );
            break;
        }
        case GAUSS_POINTS: {
            // Get GP coords for this mesh
            getGPCoords( currentMesh, theseCoords );
            break;
        }
        }
        // Add new coordinates to list
        allCoords.insert( allCoords.end(), theseCoords.begin(), theseCoords.end() );

        AMP::plog << "Found " << theseCoords.size() / 3 << " coordinates on this mesh" << std::endl;
    }

    AMP::plog << "Found " << allCoords.size() / 3 << " coordinates on all meshes" << std::endl;

    // Gives coordinates to Coupler
    unsigned int numCoords = allCoords.size() / 3;
    AMP_ASSERT( numCoords == r->getLocalSize() );
    // double relTol=1.0e-10, absTol=1.0e-10;
    d_coupler->locate_points( &allCoords[0], numCoords );
    // d_coupler->locate_points( &allCoords[0], numCoords,
    //                           relTol,       absTol );

    // Interpolate
    Vec_Dbl outputVar( numCoords, 0.0 );
    d_coupler->interpolate( moab::Coupler::LINEAR_FE, d_mapVar, &outputVar[0] );

    // This block was here for debugging, not ready to abandon it yet
    /*
    double minx= 100.0;
    double maxx=-100.0;
    double miny= 100.0;
    double maxy=-100.0;
    double minz= 100.0;
    double maxz=-100.0;
    double minp= 100.0;
    double maxp=-100.0;
    AMP::pout << "Interpolated values" << std::endl;
    for( unsigned int i=0; i<numCoords; ++i )
    {
        AMP::plog << allCoords[3*i] << " " << allCoords[3*i+1] << " " << allCoords[3*i+2] << " " <<
    outputVar[i] <<
    std::endl;
        minx = std::min( minx, allCoords[3*i] );
        maxx = std::max( maxx, allCoords[3*i] );
        miny = std::min( miny, allCoords[3*i+1] );
        maxy = std::max( maxy, allCoords[3*i+1] );
        minz = std::min( minz, allCoords[3*i+2] );
        maxz = std::max( maxz, allCoords[3*i+2] );
        minp = std::min( minp, outputVar[i] );
        maxp = std::max( maxp, outputVar[i] );
    }
    AMP::plog << "x Range: " << minx << " " << maxx << std::endl;
    AMP::plog << "y Range: " << miny << " " << maxy << std::endl;
    AMP::plog << "z Range: " << minz << " " << maxz << std::endl;
    AMP::plog << "p Range: " << minp << " " << maxp << std::endl;
    */

    // Copy values into r
    std::vector<size_t> myIndices( numCoords );
    for ( unsigned int i = 0; i < numCoords; ++i ) {
        myIndices[i] = i;
    }
    r->setValuesByLocalID( numCoords, &myIndices[0], &outputVar[0] );

    // Make consistent
    r->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
}

//---------------------------------------------------------------------------//
/*!
 *\brief Get vector of Gauss points for single mesh
 */
//---------------------------------------------------------------------------//
void MoabMapOperator::getGPCoords( AMP::Mesh::Mesh::shared_ptr &mesh, Vec_Dbl &xyz )
{
    AMP_INSIST( mesh, "Must have a mesh" );
    AMP_INSIST( d_interpType == GAUSS_POINTS, "Wrong interpolation type" );

    // Create Gauss point DOF manager
    size_t DOFsPerElement    = 8;
    int gaussPointGhostWidth = 0;
    bool split               = true;
    AMP::Discretization::DOFManager::shared_ptr gaussPointDofMap =
        AMP::Discretization::simpleDOFManager::create(
            mesh, AMP::Mesh::Volume, gaussPointGhostWidth, DOFsPerElement, split );

    // Get size of Gauss-point vectors
    // We're explicitly assuming every element has 8 Gauss points
    unsigned int numGauss = DOFsPerElement * mesh->numLocalElements( AMP::Mesh::Volume );

    // Resize vector
    xyz.resize( 3 * numGauss, 0.0 );

    // Create a variable for the coordinates
    AMP::LinearAlgebra::Variable::shared_ptr gpVariable(
        new AMP::LinearAlgebra::Variable( "coords" ) );

    // Convert from distance in m (AMP) to cm (Moab)
    // This should probably be specified on an input database
    //  rather than hard-coded here because I would guess that
    //  Moab meshes aren't all going to be in cm.
    double m_to_cm = 100.0;

    // Build Volume Integral Operator
    SP_VolIntOp volIntOp;
    buildVolumeIntOp( volIntOp, mesh );
    AMP_ASSERT( volIntOp );

    // Get FE Base from volume integral operator
    SP_FEBase fe_ptr = volIntOp->getSourceElement()->getFEBase();

    // Extract coordinates of each Gauss point
    unsigned int zeroGhostWidth  = 0;
    AMP::Mesh::MeshIterator elem = mesh->getIterator( AMP::Mesh::Volume, zeroGhostWidth );
    int elem_ctr = 0, gp_ctr = 0;
    for ( ; elem != elem.end(); elem++ ) {
        std::vector<AMP::Mesh::MeshElement> currNodes;
        currNodes = elem->getElements( AMP::Mesh::Vertex );
        ::Elem *currElemPtr;
        currElemPtr = new ::Hex8;
        for ( size_t j = 0; j < currNodes.size(); j++ ) {
            std::vector<double> pt     = currNodes[j].coord();
            currElemPtr->set_node( j ) = new ::Node( pt[0], pt[1], pt[2], j );
        } // end for j

        // Initialize FEBase for this object
        fe_ptr->reinit( currElemPtr );

        AMP_ASSERT( fe_ptr );

        std::vector<::Point> this_xyz = fe_ptr->get_xyz();

        // Loop over all Gauss-points on the element.
        for ( unsigned int i = 0; i < DOFsPerElement; i++ ) {
            // Assign coordinates
            xyz[3 * gp_ctr]     = this_xyz[i]( 0 ) * m_to_cm;
            xyz[3 * gp_ctr + 1] = this_xyz[i]( 1 ) * m_to_cm;
            xyz[3 * gp_ctr + 2] = this_xyz[i]( 2 ) * m_to_cm;
            gp_ctr++;
        }

        for ( size_t j = 0; j < currElemPtr->n_nodes(); j++ ) {
            delete ( currElemPtr->get_node( j ) );
            currElemPtr->set_node( j ) = NULL;
        } // end for j
        delete currElemPtr;
        currElemPtr = NULL;

        elem_ctr++;
    }
}

//---------------------------------------------------------------------------//
/*!
 *\brief Get vector of node coordinates for single mesh
 */
//---------------------------------------------------------------------------//
void MoabMapOperator::getNodeCoords( AMP::Mesh::Mesh::shared_ptr &mesh, Vec_Dbl &xyz )
{
    AMP_INSIST( mesh, "Must have Mesh Adapter" );
    AMP_INSIST( d_interpType == NODES, "Wrong interpolation type" );

    // Get size of nodal vectors
    unsigned int numNodes = mesh->numLocalElements( AMP::Mesh::Vertex );

    // Resize vector
    xyz.resize( 3 * numNodes, 0.0 );

    // Create Gauss point variable
    AMP::LinearAlgebra::Variable::shared_ptr gpVariable(
        new AMP::LinearAlgebra::Variable( "coords" ) );

    // Convert from distance in m (AMP) to cm (Moab)
    double m_to_cm = 100.0;

    // Extract coordinates of each node
    unsigned int zeroGhostWidth  = 0;
    AMP::Mesh::MeshIterator node = mesh->getIterator( AMP::Mesh::Vertex, zeroGhostWidth );
    int node_ctr                 = 0;
    for ( ; node != node.end(); node++ ) {
        xyz[3 * node_ctr]     = ( node->coord() )[0] * m_to_cm;
        xyz[3 * node_ctr + 1] = ( node->coord() )[1] * m_to_cm;
        xyz[3 * node_ctr + 2] = ( node->coord() )[2] * m_to_cm;
        node_ctr++;
    }
}

//---------------------------------------------------------------------------//
/*!
 *\brief Build volume integral operator
 */
//---------------------------------------------------------------------------//
void MoabMapOperator::buildVolumeIntOp( SP_VolIntOp &volIntOp, AMP::Mesh::Mesh::shared_ptr &mesh )
{
    using AMP::Operator::OperatorBuilder;

    std::string name = "VolumeIntegral";
    d_params->d_db->putDatabase( name );

    // Create volume database
    SP_Database volume_db = d_params->d_db->getDatabase( name );
    volume_db->putString( "name", "VolumeIntegralOperator" );
    volume_db->putString( "InputVariableType", "IntegrationPointScalar" );
    volume_db->putInteger( "Number_Active_Variables", 1 );
    volume_db->putInteger( "Number_Auxillary_Variables", 0 );
    volume_db->putBool( "Constant_Source", 1 );
    volume_db->putString( "OutputVariable", "VolumeIntegrated" );
    volume_db->putInteger( "print_info_level", 1 );
    volume_db->putDatabase( "ActiveInputVariables" );
    volume_db->putDatabase( "SourceElement" );

    // Source db
    SP_Database source_db = volume_db->getDatabase( "SourceElement" );
    source_db->putString( "name", "SourceNonlinearElement" );

    // Active variable db
    SP_Database act_db;
    act_db = volume_db->getDatabase( "ActiveInputVariables" );

    // Define active variable as Specific Power
    std::string interfaceVarName = "SpecificPowerInWattsPerGram";
    act_db->putString( "ActiveVariable_0", interfaceVarName );

    // Global DB
    SP_InpDatabase global_db = AMP::dynamic_pointer_cast<InpDatabase>( d_params->d_db );

    // We just need a dummy Element Physics Model
    SP_ElemPhysModel emptyModel;

    // Create the operator
    volIntOp = AMP::dynamic_pointer_cast<VolIntOp>(
        OperatorBuilder::createOperator( mesh, name, global_db, emptyModel ) );

    AMP_ASSERT( volIntOp );
}

//---------------------------------------------------------------------------//
/*!
 *\brief Build Moab Coupler
 */
//---------------------------------------------------------------------------//
void MoabMapOperator::buildMoabCoupler()
{
    // Get ParallelComm from Interface
    std::vector<moab::ParallelComm *> pcomm_vec;
    AMP::plog << "Getting ParallelComm" << std::endl;
    moab::ParallelComm::get_all_pcomm( d_moabInterface, pcomm_vec );

    AMP::plog << "Retrieved " << pcomm_vec.size() << " communicators" << std::endl;

    // Make sure we got exactly one parallel comm
    AMP_INSIST( pcomm_vec.size() == 1, "Must have exactly one Moab ParallelComm" );
    moab::ParallelComm *moabParComm = pcomm_vec[0];

    // Get source elements
    moab::Range srcElems;
    moab::ErrorCode moabError = moabParComm->get_part_entities( srcElems, 3 );
    AMP_ASSERT( moabError == moab::MB_SUCCESS );

    // Build Coupler
    int couplerID = 0;
    AMP::plog << "Calling Coupler constructor" << std::endl;
    d_coupler = AMP::shared_ptr<moab::Coupler>(
        new moab::Coupler( d_moabInterface, moabParComm, srcElems, couplerID ) );
}


} // namespace Operator
} // namespace AMP

//---------------------------------------------------------------------------//
//       end of MoabMapOperator.cc
//---------------------------------------------------------------------------//

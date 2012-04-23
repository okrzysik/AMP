//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MoabMapOperator.cc
 * \author Steven Hamilton
 * \brief  Member definitions for MoabMapOperator
 */
//---------------------------------------------------------------------------//

#include "MoabMapOperator.h"

namespace AMP {
namespace Operator {

//---------------------------------------------------------------------------//
/*!
 *\brief Constructor
 */
//---------------------------------------------------------------------------//
MoabMapOperator::MoabMapOperator( const SP_MoabMapParams &params )
    : Base(     params )
    , d_params( params )
    , d_moab( params->d_moabOp )
{
    // Get parameters from DB
    d_mapVar = params->d_db->getString("MoabMapVariable");

    // Get mesh manager
    d_meshMgr = params->d_mesh;

    // Get Moab Interface
    d_moabInterface = d_moab->getMoabInterface();

    // Interpolate to nodes or Gauss points?
    AMP::pout << "Interpolate type is " << params->d_db->getString("InterpolateToType") << std::endl;
    if( params->d_db->getString("InterpolateToType").compare("Vertex")==0 )
    {
        AMP::pout << "Interpolation type is nodes" << std::endl;
        d_interpType = NODES;
    }
    else if( params->d_db->getString("InterpolateToType").compare("GaussPoint")==0 )
    {
        AMP::pout << "Interpolation type is Gauss points" << std::endl;
        d_interpType = GAUSS_POINTS;
    }
    else
        AMP_ERROR("InterpolateToType must be Vertex or GaussPoint");

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
void MoabMapOperator::apply( const SP_Vector &f,
                             const SP_Vector &u,
                                   SP_Vector &r,
                                   double     a,
                                   double     b )
{
    AMP_INSIST(r,"Vector r must not be null");
    AMP_INSIST(d_meshMgr,"Must have Mesh Adapter"); 
    AMP_INSIST(d_moab,"Must have Moab Operator");

    // Build Moab Coupler
    AMP::pout << "Building Moab Coupler" << std::endl;
    buildMoabCoupler();
    AMP::pout << "Have Moab Coupler" << std::endl;

    // Create vector to hold coordinates
    std::vector<double> allCoords;

    // Loop over meshes
    AMP::Mesh::MeshManager::MeshIterator currentMesh;
    for( currentMesh  = d_meshMgr->beginMeshes();
         currentMesh != d_meshMgr->endMeshes();
         currentMesh++ )
    {
        AMP::pout << "Getting coords for mesh" << std::endl;
        AMP::pout << "Interpolate enum is " << d_interpType << std::endl;

        std::vector<double> theseCoords;
        switch( d_interpType )
        {
            case NODES:
            {
                // Get nodes coords
                AMP::pout << "Getting node coords" << std::endl;
                getNodeCoords( *currentMesh, theseCoords );
                break;
            }
            case GAUSS_POINTS:
            {
                // Get GP coords for this mesh
                AMP::pout << "Getting GP coords" << std::endl;
                getGPCoords( *currentMesh, theseCoords );
                break;
            }
        }
        // Add new coordinates to list
        allCoords.insert( allCoords.end(), 
                          theseCoords.begin(), 
                          theseCoords.end() );

        AMP::plog << "Found " << theseCoords.size()/3 << " coordinates on this mesh" << std::endl;
    }

    AMP::plog << "Found " << allCoords.size()/3 << " coordinates on all meshes" << std::endl;

    // Gives coordinates to Coupler
    unsigned int numCoords = allCoords.size() / 3;
    AMP_ASSERT( numCoords == r->getLocalSize() );
    d_coupler->locate_points( &allCoords[0], numCoords );

    // Interpolate
    Vec_Dbl outputVar(numCoords,0.0);
    d_coupler->interpolate( moab::Coupler::LINEAR_FE, 
                            d_mapVar, 
                            &outputVar[0] );

    AMP::pout << "Interpolated values" << std::endl;
    for( unsigned int i=0; i<numCoords; ++i )
    {
        AMP::plog << allCoords[3*i] << " " << allCoords[3*i+1] << " " << allCoords[3*i+2] << " " << outputVar[i] << std::endl;
    }

    // Copy values into r
    std::copy( outputVar.begin(), outputVar.end(), r->begin() );

    /*
    AMP::Vector::MultiVector::vector_iterator currentVector;
    for( currentVector  = r->beginVector();
         currentVector != r->endVector();
         currentVector++ )
    {
    }
    */


}

//---------------------------------------------------------------------------//
/*!
 *\brief Get vector of Gauss points for single mesh
 */
//---------------------------------------------------------------------------//
void MoabMapOperator::getGPCoords( SP_Mesh &mesh, Vec_Dbl &xyz )
{
    AMP_INSIST(mesh,"Must have Mesh Adapter"); 
    AMP_INSIST(d_interpType==GAUSS_POINTS,"Wrong interpolation type");

    // Get size of Gauss-point vectors
    // We're explicitly assuming every element has 8 Gauss points
    unsigned int numGauss = 8*mesh->numLocalElements();

    // Resize vector
    xyz.resize(3*numGauss,0.0);

    // Create Gauss point variable
    SP_Variable gpVariable(new HexGPVar("coords", mesh));

    // Convert from distance in m (AMP) to cm (Moab)
    double m_to_cm = 100.0;

    // Build Volume Integral Operator
    SP_VolIntOp volIntOp;
    buildVolumeIntOp( volIntOp, mesh );
    AMP_ASSERT(volIntOp);

    // Get FE Base from volume integral operator
    SP_FEBase fe_ptr = volIntOp->getSourceElement()->getFEBase();
    
    // Extract coordinates of each Gauss point
    AMP::Mesh::MeshAdapter::ElementIterator elem = mesh->beginElement();
    int elem_ctr=0, gp_ctr=0;
    for( ; elem != mesh->endElement();
           elem++ )
    {
        // Initialize FEBase for this object
        fe_ptr->reinit( &(elem->getElem()) );

        AMP_ASSERT( fe_ptr );

        std::vector< ::Point > this_xyz = fe_ptr->get_xyz();

        // Loop over all Gauss-points on the element.
        for( unsigned int i = 0; i < 8; i++ ) 
        {
            // Assign coordinates
            xyz[3*gp_ctr]   = this_xyz[i](0) * m_to_cm;
            xyz[3*gp_ctr+1] = this_xyz[i](1) * m_to_cm;
            xyz[3*gp_ctr+2] = this_xyz[i](2) * m_to_cm;
            gp_ctr++;
        }
        elem_ctr++;
    }

}

//---------------------------------------------------------------------------//
/*!
 *\brief Get vector of node coordinates for single mesh
 */
//---------------------------------------------------------------------------//
void MoabMapOperator::getNodeCoords( SP_Mesh &mesh, Vec_Dbl &xyz )
{
    AMP_INSIST(mesh,"Must have Mesh Adapter"); 
    AMP_INSIST(d_interpType==NODES,"Wrong interpolation type" );

    // Get size of nodal vectors
    unsigned int numNodes = mesh->numLocalNodes();

    // Resize vector
    xyz.resize(3*numNodes,0.0);

    // Create Gauss point variable
    SP_Variable gpVariable(new HexGPVar("coords", mesh));

    // Convert from distance in m (AMP) to cm (Moab)
    double m_to_cm = 100.0;

    // Extract coordinates of each node
    AMP::Mesh::MeshAdapter::NodeIterator node = mesh->beginNode();
    int node_ctr=0;
    for( ; node != mesh->endNode(); node++ )
    {
        xyz[3*node_ctr]   = node->x() * m_to_cm;
        xyz[3*node_ctr+1] = node->y() * m_to_cm;
        xyz[3*node_ctr+2] = node->z() * m_to_cm;
        node_ctr++;
    }

}

//---------------------------------------------------------------------------//
/*!
 *\brief Get GP coordinates for single mesh
 */
//---------------------------------------------------------------------------//
/*
void MoabMapOperator::getGPCoords( SP_Mesh &mesh, std::vector<double> &xyz )
{
    AMP_INSIST(mesh,"Must have mesh adapter");

    // Get size of Gauss-point vectors
    // We're explicitly assuming every element has 8 Gauss points
    unsigned int numGauss = 8*mesh->numLocalElements();

    // Create Gauss point variable
    SP_Variable gpVariable(new HexGPVar("coords", mesh));

    // Build Volume Integral Operator
    SP_VolIntOp volIntOp;
    buildVolumeIntOp( volIntOp, mesh );
    AMP_ASSERT(volIntOp);

    // Get FE Base from volume integral operator
    SP_FEBase fe_ptr = volIntOp->getSourceElement()->getFEBase();

    // Now copy powers back to AMP vector
    // If GP Vectors are guaranteed to be stored as
    //  index = gp + (8 * element) then this could
    //  be a straight copy, otherwise we have to loop
    //  over elements and use the DOF Map
    AMP::Mesh::MeshAdapter::ElementIterator elem = mesh->beginElement();
    int elem_ctr = 0;
    int gp_ctr   = 0;
    for( elem  = mesh->beginElement(); 
         elem != mesh->endElement();
         elem++ )
    {
        // Get DOF map for this element
        AMP::Mesh::DOFMap::shared_ptr dof_map = 
            mesh->getDOFMap( gpVariable );

        AMP_ASSERT( dof_map );

        std::vector<unsigned int> ndx;
        std::vector<unsigned int> empty;
        dof_map->getDOFs( *elem, ndx, empty );

        // Initialize FEBase for this object
        fe_ptr->reinit( &(elem->getElem()) );

        AMP_ASSERT( fe_ptr );

        // Loop over all Gauss-points on the element.
        for( unsigned int i = 0; i < 8; i++ ) 
        {
            int  offset = ndx[i];

            // Assign power
            r->setValueByGlobalID( offset, powers[gp_ctr] );
            gp_ctr++;
        }
        elem_ctr++;
    }

}
*/


//---------------------------------------------------------------------------//
/*!
 *\brief Build volume integral operator
 */
//---------------------------------------------------------------------------//
void MoabMapOperator::buildVolumeIntOp( SP_VolIntOp &volIntOp,
                                        SP_Mesh     &mesh )
{
    using AMP::Operator::OperatorBuilder;

    std::string name = "VolumeIntegral";
    d_params->d_db->putDatabase(name);

    // Create volume database
    SP_Database volume_db = d_params->d_db->getDatabase(name);
    volume_db->putString("name","VolumeIntegralOperator");
    volume_db->putString("InputVariableType","IntegrationPointScalar");
    volume_db->putInteger("Number_Active_Variables",1);
    volume_db->putInteger("Number_Auxillary_Variables",0);
    volume_db->putBool("Constant_Source",1);
    volume_db->putString("OutputVariable","VolumeIntegrated");
    volume_db->putInteger("print_info_level",1);
    volume_db->putDatabase("ActiveInputVariables");
    volume_db->putDatabase("SourceElement");

    // Source db
    SP_Database source_db = volume_db->getDatabase("SourceElement");
    source_db->putString("name", "SourceNonlinearElement");

    // Active variable db
    SP_Database act_db;
    act_db = volume_db->getDatabase("ActiveInputVariables");

    // Define active variable as Specific Power
    std::string interfaceVarName = "SpecificPowerInWattsPerGram";
    act_db->putString("ActiveVariable_0",interfaceVarName);

    // Global DB
    SP_InpDatabase global_db = 
        boost::dynamic_pointer_cast<InpDatabase>(d_params->d_db);

    // We just need a dummy Element Physics Model
    SP_ElemPhysModel emptyModel;

    // Create the operator
    volIntOp = boost::dynamic_pointer_cast<VolIntOp>
        (OperatorBuilder::createOperator(mesh, 
                                         name, 
                                         global_db, 
                                         emptyModel));

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
    std::vector< moab::ParallelComm* > pcomm_vec;
    AMP::plog << "Getting ParallelComm" << std::endl;
    moab::ParallelComm::get_all_pcomm( d_moabInterface, pcomm_vec );

    AMP::plog << "Retrieved " << pcomm_vec.size() << " communicators" << std::endl;

    // Make sure we got exactly one parallel comm
    AMP_INSIST( pcomm_vec.size() == 1,"Must have exactly one Moab ParallelComm" );
    moab::ParallelComm *moabParComm = pcomm_vec[0];

    // Get source elements
    moab::Range srcElems;
    moab::ErrorCode moabError = moabParComm->get_part_entities( srcElems, 3 );
    AMP_ASSERT( moabError == moab::MB_SUCCESS );

    // Build Coupler
    int couplerID = 0;
    AMP::plog << "Calling Coupler constructor" << std::endl;
    d_coupler = boost::shared_ptr< moab::Coupler >( 
            new moab::Coupler( d_moabInterface,
                               moabParComm,
                               srcElems,
                               couplerID ) );

}


} // namespace Operator
} // namespace AMP

//---------------------------------------------------------------------------//
//       end of MoabMapOperator.cc
//---------------------------------------------------------------------------//

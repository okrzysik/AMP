#include "AMP/operators/subchannel/FlowFrapconOperator.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"

#include <string>


namespace AMP {
namespace Operator {


FlowFrapconOperator::FlowFrapconOperator(
    std::shared_ptr<const FlowFrapconOperatorParameters> params )
    : Operator( params ), d_boundaryId( 0 )
{
    std::string inpVar = params->d_db->getString( "InputVariable" );
    d_inpVariable.reset( new AMP::LinearAlgebra::Variable( inpVar ) );

    std::string outVar = params->d_db->getString( "OutputVariable" );
    d_outVariable.reset( new AMP::LinearAlgebra::Variable( outVar ) );

    reset( params );
}

void FlowFrapconOperator::reset( std::shared_ptr<const OperatorParameters> params )
{
    auto myparams = std::dynamic_pointer_cast<const FlowFrapconOperatorParameters>( params );

    AMP_INSIST( ( ( myparams.get() ) != nullptr ), "NULL parameters" );
    AMP_INSIST( ( ( ( myparams->d_db ).get() ) != nullptr ), "NULL database" );

    AMP_INSIST( myparams->d_db->keyExists( "numpoints" ), "Key ''numpoints'' is missing!" );
    d_numpoints = myparams->d_db->getScalar<int>( "numpoints" );

    AMP_INSIST( myparams->d_db->keyExists( "Channel_Diameter" ), "Missing key: Channel_Dia" );
    d_De = myparams->d_db->getScalar<double>( "Channel_Diameter" );

    AMP_INSIST( myparams->d_db->keyExists( "Heat_Capacity" ), "Missing key: Heat_Capacity" );
    Cp = myparams->d_db->getScalar<double>( "Heat_Capacity" );

    AMP_INSIST( myparams->d_db->keyExists( "Mass_Flux" ), "Missing key: Mass_Flux" );
    d_G = myparams->d_db->getScalar<double>( "Mass_Flux" );

    AMP_INSIST( myparams->d_db->keyExists( "Temp_Inlet" ), "Missing key: Temp_In" );
    d_Tin = myparams->d_db->getWithDefault<double>( "Temp_Inlet", 300. );

    AMP_INSIST( myparams->d_db->keyExists( "Conductivity" ), "Missing key: Kconductivity" );
    d_K = myparams->d_db->getScalar<double>( "Conductivity" );

    AMP_INSIST( myparams->d_db->keyExists( "Reynolds" ), "Missing key: Reynolds" );
    d_Re = myparams->d_db->getScalar<double>( "Reynolds" );

    AMP_INSIST( myparams->d_db->keyExists( "Prandtl" ), "Missing key: Prandtl" );
    d_Pr = myparams->d_db->getScalar<double>( "Prandtl" );
}


// This is an in-place apply
void FlowFrapconOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                 AMP::LinearAlgebra::Vector::shared_ptr r )
{

    // AMP::Mesh::DOFMap::shared_ptr dof_map = d_MeshAdapter->getDOFMap(d_inpVariable);

    AMP_INSIST( ( ( r.get() ) != nullptr ), "NULL Residual Vector" );
    AMP_INSIST( ( ( u.get() ) != nullptr ), "NULL Solution Vector" );

    if ( !zPoints.empty() ) {
        d_numpoints = zPoints.size();
    }

    std::vector<double> box = d_Mesh->getBoundingBox();
    const double min_z      = box[4];
    const double max_z      = box[5];
    const double del_z      = ( max_z - min_z ) / d_numpoints;

    // std::cout << "Extreme Min Point in z = " << min_z << std::endl;
    // std::cout << "Extreme Max Point in z = " << max_z << std::endl;

    // Subset the vectors
    AMP::LinearAlgebra::Vector::const_shared_ptr flowInputVec = subsetInputVector( u );
    AMP::LinearAlgebra::Vector::shared_ptr outputVec          = subsetOutputVector( r );
    zPoints.resize( d_numpoints );

    // set the inlet flow temperature value
    double T1        = flowInputVec->getValueByLocalID( 0 );
    size_t idx       = 0;
    const double val = T1 - d_Tin;
    outputVec->setValuesByLocalID( 1, &idx, &val );

    zPoints[0] = min_z;
    for ( int j = 1; j < d_numpoints; j++ ) {
        zPoints[j] = zPoints[j - 1] + del_z;
    }

    // Iterate through the flow boundary
    for ( size_t i = 1; i < (size_t) d_numpoints; i++ ) {

        double cur_node, next_node;

        cur_node  = zPoints[i - 1];
        next_node = zPoints[i];

        double Heff, he_z, T_b_i, T_b_im1, T_c_i;
        double R_b = 0;

        T_c_i   = d_cladVec->getValueByLocalID( i );
        T_b_i   = flowInputVec->getValueByLocalID( i );
        T_b_im1 = flowInputVec->getValueByLocalID( i - 1 );

        Heff = ( 0.023 * d_K / d_De ) * pow( d_Re, 0.8 ) * pow( d_Pr, 0.4 );
        //       Cp   = getHeatCapacity(T_b_i);
        he_z = next_node - cur_node;

        R_b = T_b_i - T_b_im1 - ( ( 4 * Heff * ( T_c_i - T_b_i ) ) / ( Cp * d_G * d_De ) ) * he_z;

        outputVec->setValuesByLocalID( 1, &i, &R_b );

    } // end for i
}


std::shared_ptr<OperatorParameters>
FlowFrapconOperator::getJacobianParameters( AMP::LinearAlgebra::Vector::const_shared_ptr u_in )
{
    auto tmp_db = std::make_shared<AMP::Database>( "Dummy" );

    tmp_db->putScalar( "name", "FlowFrapconOperator" );
    tmp_db->putScalar( "numpoints", d_numpoints );
    tmp_db->putScalar( "Channel_Diameter", d_De );
    tmp_db->putScalar( "Mass_Flux", d_G );
    tmp_db->putScalar( "Heat_Capacity", Cp );
    tmp_db->putScalar( "Temp_Inlet", d_Tin );
    tmp_db->putScalar( "Conductivity", d_K );
    tmp_db->putScalar( "Reynolds", d_Re );
    tmp_db->putScalar( "Prandtl", d_Pr );

    auto outParams              = std::make_shared<FlowFrapconJacobianParameters>( tmp_db );
    auto u                      = std::const_pointer_cast<AMP::LinearAlgebra::Vector>( u_in );
    outParams->d_frozenSolution = subsetInputVector( u );
    return outParams;
}


AMP::LinearAlgebra::Vector::shared_ptr
FlowFrapconOperator::subsetOutputVector( AMP::LinearAlgebra::Vector::shared_ptr vec )
{
    AMP::LinearAlgebra::Variable::shared_ptr var = getInputVariable();
    // Subset the vectors, they are simple vectors and we need to subset for the current comm
    // instead of the mesh
    if ( d_Mesh ) {
        AMP::LinearAlgebra::VS_Comm commSelector( d_Mesh->getComm() );
        AMP::LinearAlgebra::Vector::shared_ptr commVec =
            vec->select( commSelector, var->getName() );
        return commVec->subsetVectorForVariable( var );
    } else {
        return vec->subsetVectorForVariable( var );
    }
}


AMP::LinearAlgebra::Vector::shared_ptr
FlowFrapconOperator::subsetInputVector( AMP::LinearAlgebra::Vector::shared_ptr vec )
{
    AMP::LinearAlgebra::Variable::shared_ptr var = getInputVariable();
    // Subset the vectors, they are simple vectors and we need to subset for the current comm
    // instead of the mesh
    if ( d_Mesh ) {
        AMP::LinearAlgebra::VS_Comm commSelector( d_Mesh->getComm() );
        AMP::LinearAlgebra::Vector::shared_ptr commVec =
            vec->select( commSelector, var->getName() );
        return commVec->subsetVectorForVariable( var );
    } else {
        return vec->subsetVectorForVariable( var );
    }
}


AMP::LinearAlgebra::Vector::const_shared_ptr
FlowFrapconOperator::subsetOutputVector( AMP::LinearAlgebra::Vector::const_shared_ptr vec )
{
    AMP::LinearAlgebra::Variable::shared_ptr var = getInputVariable();
    // Subset the vectors, they are simple vectors and we need to subset for the current comm
    // instead of the mesh
    if ( d_Mesh ) {
        AMP::LinearAlgebra::VS_Comm commSelector( d_Mesh->getComm() );
        AMP::LinearAlgebra::Vector::const_shared_ptr commVec =
            vec->constSelect( commSelector, var->getName() );
        return commVec->constSubsetVectorForVariable( var );
    } else {
        return vec->constSubsetVectorForVariable( var );
    }
}


AMP::LinearAlgebra::Vector::const_shared_ptr
FlowFrapconOperator::subsetInputVector( AMP::LinearAlgebra::Vector::const_shared_ptr vec )
{
    AMP::LinearAlgebra::Variable::shared_ptr var = getInputVariable();
    // Subset the vectors, they are simple vectors and we need to subset for the current comm
    // instead of the mesh
    if ( d_Mesh ) {
        AMP::LinearAlgebra::VS_Comm commSelector( d_Mesh->getComm() );
        AMP::LinearAlgebra::Vector::const_shared_ptr commVec =
            vec->constSelect( commSelector, var->getName() );
        return commVec->constSubsetVectorForVariable( var );
    } else {
        return vec->constSubsetVectorForVariable( var );
    }
}
} // namespace Operator
} // namespace AMP

#include "AMP/operators/subchannel/FlowFrapconJacobian.h"
#include "AMP/vectors/VectorSelector.h"

#include <string>


namespace AMP::Operator {


FlowFrapconJacobian::FlowFrapconJacobian( std::shared_ptr<const OperatorParameters> in_params )
    : Operator( in_params ), dCp( 0 )
{
    auto params = std::dynamic_pointer_cast<const FlowFrapconJacobianParameters>( in_params );
    std::string inpVar = params->d_db->getString( "InputVariable" );
    d_inpVariable.reset( new AMP::LinearAlgebra::Variable( inpVar ) );

    std::string outVar = params->d_db->getString( "OutputVariable" );
    d_outVariable.reset( new AMP::LinearAlgebra::Variable( outVar ) );

    d_SimpleVariable.reset( new AMP::LinearAlgebra::Variable( "FlowInternal" ) );

    reset( params );
}


std::shared_ptr<AMP::LinearAlgebra::Variable>
FlowFrapconJacobian::createInputVariable( const std::string &name, int varId )
{
    (void) varId;
    return d_inpVariable->clone( name );
}


std::shared_ptr<AMP::LinearAlgebra::Variable>
FlowFrapconJacobian::createOutputVariable( const std::string &name, int varId )
{
    (void) varId;
    return d_outVariable->clone( name );
}


std::shared_ptr<AMP::LinearAlgebra::Variable> FlowFrapconJacobian::getInputVariable()
{
    return d_inpVariable;
}


std::shared_ptr<AMP::LinearAlgebra::Variable> FlowFrapconJacobian::getOutputVariable()
{
    return d_outVariable;
}


void FlowFrapconJacobian::reset( std::shared_ptr<const OperatorParameters> params )
{
    AMP_ASSERT( params );

    auto myparams = std::dynamic_pointer_cast<const FlowFrapconJacobianParameters>( params );

    AMP_INSIST( ( ( myparams.get() ) != nullptr ), "NULL parameters" );
    AMP_INSIST( ( ( ( myparams->d_db ).get() ) != nullptr ), "NULL database" );

    bool skipParams = params->d_db->getWithDefault<bool>( "skip_params", false );

    if ( !skipParams ) {
        AMP_INSIST( myparams->d_db->keyExists( "numpoints" ), "Key ''numpoints'' is missing!" );
        d_numpoints = myparams->d_db->getScalar<int>( "numpoints" );

        AMP_INSIST( myparams->d_db->keyExists( "Channel_Diameter" ), "Missing key: Channel_Dia" );
        d_De = myparams->d_db->getScalar<double>( "Channel_Diameter" );

        AMP_INSIST( myparams->d_db->keyExists( "Heat_Capacity" ), "Missing key: Heat_Capacity" );
        Cp = myparams->d_db->getScalar<double>( "Heat_Capacity" );

        AMP_INSIST( myparams->d_db->keyExists( "Mass_Flux" ), "Missing key: Mass_Flux" );
        d_G = myparams->d_db->getScalar<double>( "Mass_Flux" );

        AMP_INSIST( myparams->d_db->keyExists( "Temp_Inlet" ), "Missing key: Temp_In" );
        d_Tin = myparams->d_db->getScalar<double>( "Temp_Inlet" );

        AMP_INSIST( myparams->d_db->keyExists( "Conductivity" ), "Missing key: Kconductivity" );
        d_K = myparams->d_db->getScalar<double>( "Conductivity" );

        AMP_INSIST( myparams->d_db->keyExists( "Reynolds" ), "Missing key: Reynolds" );
        d_Re = myparams->d_db->getScalar<double>( "Reynolds" );

        AMP_INSIST( myparams->d_db->keyExists( "Prandtl" ), "Missing key: Prandtl" );
        d_Pr = myparams->d_db->getScalar<double>( "Prandtl" );
    }

    if ( ( myparams->d_frozenSolution.get() ) != nullptr ) {
        d_frozenVec = myparams->d_frozenSolution;
    }
}


// This is an in-place apply
void FlowFrapconJacobian::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                 AMP::LinearAlgebra::Vector::shared_ptr r )
{

    // AMP::Mesh::DOFMap::shared_ptr dof_map = d_MeshAdapter->getDOFMap(d_inpVariable);

    AMP_INSIST( ( ( r.get() ) != nullptr ), "NULL Residual Vector" );
    AMP_INSIST( ( ( u.get() ) != nullptr ), "NULL Solution Vector" );

    std::vector<double> box = d_Mesh->getBoundingBox();
    const double min_z      = box[4];
    const double max_z      = box[5];
    const double del_z      = ( max_z - min_z ) / d_numpoints;

    //    std::cout << "Extreme Min Point in z = " << min_z << std::endl;
    //    std::cout << "Extreme Max Point in z = " << max_z << std::endl;


    AMP::LinearAlgebra::Vector::const_shared_ptr flowInputVec = subsetInputVector( u );

    AMP::LinearAlgebra::Vector::shared_ptr outputVec = subsetOutputVector( r );

    AMP_INSIST( ( d_frozenVec ), "Null Frozen Vector inside Jacobian" );

    if ( !zPoints.empty() ) {
        d_numpoints = zPoints.size();
    }

    zPoints.resize( d_numpoints );

    // set the inlet flow temperature value
    // flowInputVec->setValueByLocalID(0, 0.0);
    size_t idx       = 0;
    const double val = flowInputVec->getValueByLocalID( idx );
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

        double Heff, he_z, dT_b_i, dT_b_im1, T_c_i, T_b_i;
        double JT_b = 0;

        T_c_i = d_cladVec->getValueByLocalID( i );
        T_b_i = d_frozenVec->getValueByLocalID( i );

        dT_b_i   = flowInputVec->getValueByLocalID( i );
        dT_b_im1 = flowInputVec->getValueByLocalID( i - 1 );

        Heff = ( 0.023 * d_K / d_De ) * std::pow( d_Re, 0.8 ) * std::pow( d_Pr, 0.4 );
        //                  Cp   = getHeatCapacity(T_b_i);
        //                  dCp  = getHeatCapacityGradient(T_b_i);
        dCp  = 0.0;
        he_z = next_node - cur_node;

        JT_b =
            dT_b_i - dT_b_im1 -
            ( ( 4 * Heff * ( T_c_i - T_b_i ) ) / ( Cp * Cp * d_G * d_De ) ) * he_z * dCp * dT_b_i +
            ( ( 4 * Heff ) / ( Cp * d_G * d_De ) ) * he_z * dT_b_i;

        outputVec->setValuesByLocalID( 1, &i, &JT_b );

    } // end for i
}


AMP::LinearAlgebra::Vector::shared_ptr
FlowFrapconJacobian::subsetOutputVector( AMP::LinearAlgebra::Vector::shared_ptr vec )
{
    std::shared_ptr<AMP::LinearAlgebra::Variable> var = getInputVariable();
    // Subset the vectors, they are simple vectors and we need to subset for the current comm
    // instead of the mesh
    if ( d_Mesh ) {
        AMP::LinearAlgebra::VS_Comm commSelector( d_Mesh->getComm() );
        auto commVec = vec->select( commSelector );
        return commVec->subsetVectorForVariable( var );
    } else {
        return vec->subsetVectorForVariable( var );
    }
}


AMP::LinearAlgebra::Vector::shared_ptr
FlowFrapconJacobian::subsetInputVector( AMP::LinearAlgebra::Vector::shared_ptr vec )
{
    std::shared_ptr<AMP::LinearAlgebra::Variable> var = getInputVariable();
    // Subset the vectors, they are simple vectors and we need to subset for the current comm
    // instead of the mesh
    if ( d_Mesh ) {
        AMP::LinearAlgebra::VS_Comm commSelector( d_Mesh->getComm() );
        auto commVec = vec->select( commSelector );
        return commVec->subsetVectorForVariable( var );
    } else {
        return vec->subsetVectorForVariable( var );
    }
}


AMP::LinearAlgebra::Vector::const_shared_ptr
FlowFrapconJacobian::subsetOutputVector( AMP::LinearAlgebra::Vector::const_shared_ptr vec )
{
    std::shared_ptr<AMP::LinearAlgebra::Variable> var = getInputVariable();
    // Subset the vectors, they are simple vectors and we need to subset for the current comm
    // instead of the mesh
    if ( d_Mesh ) {
        AMP::LinearAlgebra::VS_Comm commSelector( d_Mesh->getComm() );
        auto commVec = vec->select( commSelector );
        return commVec->subsetVectorForVariable( var );
    } else {
        return vec->subsetVectorForVariable( var );
    }
}


AMP::LinearAlgebra::Vector::const_shared_ptr
FlowFrapconJacobian::subsetInputVector( AMP::LinearAlgebra::Vector::const_shared_ptr vec )
{
    std::shared_ptr<AMP::LinearAlgebra::Variable> var = getInputVariable();
    // Subset the vectors, they are simple vectors and we need to subset for the current comm
    // instead of the mesh
    if ( d_Mesh ) {
        AMP::LinearAlgebra::VS_Comm commSelector( d_Mesh->getComm() );
        auto commVec = vec->select( commSelector );
        return commVec->subsetVectorForVariable( var );
    } else {
        return vec->subsetVectorForVariable( var );
    }
}
} // namespace AMP::Operator

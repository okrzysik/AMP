
#include "MassLinearFEOperator.h"

#include "AMP/utils/Utilities.h"


namespace AMP::Operator {


std::shared_ptr<AMP::LinearAlgebra::Variable> MassLinearFEOperator::getInputVariable()
{
    return d_inpVariable;
}


std::shared_ptr<AMP::LinearAlgebra::Variable> MassLinearFEOperator::getOutputVariable()
{
    return d_outVariable;
}


MassLinearFEOperator::MassLinearFEOperator(
    std::shared_ptr<const MassLinearFEOperatorParameters> params )
    : LinearFEOperator( params )
{
    AMP_INSIST( params, "NULL parameter" );

    d_massLinElem = std::dynamic_pointer_cast<MassLinearElement>( d_elemOp );

    AMP_INSIST( ( ( d_massLinElem.get() ) != nullptr ),
                "d_elemOp is not of type MassLinearElement" );

    d_densityModel = params->d_densityModel;

    d_useConstantTemperature = params->d_db->keyExists( "FixedTemperature" );

    d_useConstantConcentration = params->d_db->keyExists( "FixedConcentration" );

    d_useConstantBurnup = params->d_db->keyExists( "FixedBurnup" );

    d_constantTemperatureValue = params->d_db->getWithDefault<double>( "FixedTemperature", 273.0 );

    d_constantConcentrationValue =
        params->d_db->getWithDefault<double>( "FixedConcentration", 0.0 );

    d_constantBurnupValue = params->d_db->getWithDefault<double>( "FixedBurnup", 0.0 );
    // d_inpVariable.reset(new AMP::Mesh::NodalScalarVariable("inpVar"));
    // d_outVariable.reset(new AMP::Mesh::NodalScalarVariable("outVar"));
    std::string inpVar = params->d_db->getString( "InputVariable" );
    d_inpVariable      = std::make_shared<AMP::LinearAlgebra::Variable>( inpVar );

    std::string outVar = params->d_db->getString( "OutputVariable" );
    d_outVariable      = std::make_shared<AMP::LinearAlgebra::Variable>( outVar );

    reset( params );
}


void MassLinearFEOperator::preAssembly( std::shared_ptr<const AMP::Operator::OperatorParameters> )
{
    d_matrix->zero();

    d_densityModel->preLinearAssembly();
}


void MassLinearFEOperator::postAssembly()
{
    d_densityModel->postLinearAssembly();

    d_matrix->makeConsistent();
}


void MassLinearFEOperator::preElementOperation( const AMP::Mesh::MeshElement &elem )
{
    d_currNodes                 = elem.getElements( AMP::Mesh::GeomType::Vertex );
    unsigned int num_local_dofs = d_currNodes.size();

    d_elementMassMatrix.resize( num_local_dofs );
    for ( unsigned int r = 0; r < num_local_dofs; r++ ) {
        d_elementMassMatrix[r].resize( num_local_dofs );
        for ( unsigned int c = 0; c < num_local_dofs; c++ ) {
            d_elementMassMatrix[r][c] = 0;
        }
    }

    std::vector<double> localTemperature( num_local_dofs );
    std::vector<double> localConcentration( num_local_dofs );
    std::vector<double> localBurnup( num_local_dofs );
    std::vector<size_t> dofs;

    if ( d_useConstantTemperature ) {
        for ( size_t r = 0; r < d_currNodes.size(); r++ ) {
            localTemperature[r] = d_constantTemperatureValue;
        }
    } else {
        std::shared_ptr<AMP::Discretization::DOFManager> DOF = d_temperature->getDOFManager();
        for ( size_t r = 0; r < d_currNodes.size(); r++ ) {
            DOF->getDOFs( d_currNodes[r].globalID(), dofs );
            AMP_ASSERT( dofs.size() == 1 );
            localTemperature[r] = d_temperature->getValueByGlobalID( dofs[0] );
        }
    }

    if ( d_useConstantConcentration ) {
        for ( size_t r = 0; r < d_currNodes.size(); r++ ) {
            localConcentration[r] = d_constantConcentrationValue;
        }
    } else {
        std::shared_ptr<AMP::Discretization::DOFManager> DOF = d_concentration->getDOFManager();
        for ( size_t r = 0; r < d_currNodes.size(); r++ ) {
            DOF->getDOFs( d_currNodes[r].globalID(), dofs );
            AMP_ASSERT( dofs.size() == 1 );
            localConcentration[r] = d_concentration->getValueByGlobalID( dofs[0] );
        }
    }

    if ( d_useConstantBurnup ) {
        for ( size_t r = 0; r < d_currNodes.size(); r++ ) {
            localBurnup[r] = d_constantBurnupValue;
        }
    } else {
        std::shared_ptr<AMP::Discretization::DOFManager> DOF = d_burnup->getDOFManager();
        for ( size_t r = 0; r < d_currNodes.size(); r++ ) {
            DOF->getDOFs( d_currNodes[r].globalID(), dofs );
            AMP_ASSERT( dofs.size() == 1 );
            localBurnup[r] = d_burnup->getValueByGlobalID( dofs[0] );
        }
    }

    createCurrentLibMeshElement();

    d_massLinElem->initializeForCurrentElement( d_currElemPtr, d_densityModel );

    d_massLinElem->setElementMassMatrix( d_elementMassMatrix );

    d_massLinElem->setElementVectors( localTemperature, localConcentration, localBurnup );
}


void MassLinearFEOperator::postElementOperation()
{
    std::vector<size_t> d_dofIndices( d_currNodes.size() ), dofs( 1 );
    for ( size_t i = 0; i < d_currNodes.size(); i++ ) {
        d_inDofMap->getDOFs( d_currNodes[i].globalID(), dofs );
        AMP_ASSERT( dofs.size() == 1 );
        d_dofIndices[i] = dofs[0];
    }

    for ( size_t r = 0; r < d_dofIndices.size(); r++ ) {
        for ( size_t c = 0; c < d_dofIndices.size(); c++ ) {
            d_matrix->addValueByGlobalID(
                d_dofIndices[r], d_dofIndices[c], d_elementMassMatrix[r][c] );
        }
    }

    destroyCurrentLibMeshElement();
}
} // namespace AMP::Operator

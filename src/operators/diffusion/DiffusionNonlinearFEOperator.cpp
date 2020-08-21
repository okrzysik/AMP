#include "AMP/operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "AMP/materials/Material.h"
#include "AMP/operators/ElementOperationParameters.h"
#include "AMP/operators/diffusion/DiffusionConstants.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperatorParameters.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"
#include "ProfilerApp.h"

#include <cstring>
#include <iostream>

namespace AMP {
namespace Operator {


AMP::LinearAlgebra::Variable::shared_ptr
DiffusionNonlinearFEOperator::createInputVariable( const std::string &name, int varId )
{
    if ( varId == -1 ) {
        return d_inpVariables->cloneVariable( name );
    } else {
        return ( d_inpVariables->getVariable( varId ) )->cloneVariable( name );
    }
}


AMP::LinearAlgebra::Variable::shared_ptr
DiffusionNonlinearFEOperator::createOutputVariable( const std::string &name, int varId )
{
    (void) varId;
    return d_outVariable->cloneVariable( name );
}


AMP::LinearAlgebra::Variable::shared_ptr DiffusionNonlinearFEOperator::getInputVariable()
{
    return d_inpVariables;
}


AMP::LinearAlgebra::Variable::shared_ptr DiffusionNonlinearFEOperator::getOutputVariable()
{
    return d_outVariable;
}

unsigned int DiffusionNonlinearFEOperator::numberOfDOFMaps() { return 1; }

AMP::LinearAlgebra::Variable::shared_ptr
DiffusionNonlinearFEOperator::getVariableForDOFMap( unsigned int id )
{
    (void) id;
    return d_inpVariables->getVariable( d_PrincipalVariable );
}


unsigned int DiffusionNonlinearFEOperator::getPrincipalVariableId() { return d_PrincipalVariable; }


std::vector<unsigned int> DiffusionNonlinearFEOperator::getNonPrincipalVariableIds()
{
    std::vector<unsigned int> ids;
    for ( size_t i = 0; i < Diffusion::NUMBER_VARIABLES; i++ ) {
        if ( i != d_PrincipalVariable and d_isActive[i] )
            ids.push_back( i );
    }
    return ids;
}


std::shared_ptr<DiffusionTransportModel> DiffusionNonlinearFEOperator::getTransportModel()
{
    return d_transportModel;
}


std::vector<AMP::LinearAlgebra::Vector::shared_ptr> DiffusionNonlinearFEOperator::getFrozen()
{
    return d_Frozen;
}


void DiffusionNonlinearFEOperator::setVector( unsigned int id,
                                              AMP::LinearAlgebra::Vector::shared_ptr frozenVec )
{
    AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
    // AMP::LinearAlgebra::Vector::shared_ptr meshSubsetVec = frozenVec->select(meshSelector,
    // d_inpVariables->getVariable(id)->getName());
    AMP::LinearAlgebra::Vector::shared_ptr meshSubsetVec =
        frozenVec->select( meshSelector, frozenVec->getVariable()->getName() );
    d_Frozen[id] = meshSubsetVec->subsetVectorForVariable( d_inpVariables->getVariable( id ) );
    ( d_Frozen[id] )->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
}


DiffusionNonlinearFEOperator::DiffusionNonlinearFEOperator(
    const std::shared_ptr<DiffusionNonlinearFEOperatorParameters> &params )
    : NonlinearFEOperator( params ), d_Frozen( Diffusion::NUMBER_VARIABLES )
{
    AMP_INSIST( ( ( params.get() ) != nullptr ), "NULL parameter!" );

    d_diffNonlinElem = std::dynamic_pointer_cast<DiffusionNonlinearElement>( d_elemOp );

    AMP_INSIST( ( ( d_diffNonlinElem.get() ) != nullptr ),
                "d_elemOp is not of type DiffusionNonlinearElement" );

    d_transportModel = params->d_transportModel;

    d_isActive.resize( Diffusion::NUMBER_VARIABLES );
    d_isFrozen.resize( Diffusion::NUMBER_VARIABLES );
    d_inVec.resize( Diffusion::NUMBER_VARIABLES );

    std::shared_ptr<AMP::Database> activeVariables_db =
        params->d_db->getDatabase( "ActiveInputVariables" );

    for ( size_t var = 0; var < Diffusion::NUMBER_VARIABLES; var++ ) {
        std::string namespec = activeVariables_db->getWithDefault<std::string>(
            Diffusion::names[var], "not_specified" );
        bool isthere           = namespec != "not_specified";
        d_isActive[var]        = isthere;
        std::string frozenName = "Freeze" + Diffusion::names[var];
        d_isFrozen[var]        = params->d_db->getWithDefault( frozenName, false );
        if ( d_isFrozen[var] )
            AMP_INSIST( d_isActive[var], "can not freeze a variable unless it is active" );
    }

    d_numberActive = 0;
    d_numberFrozen = 0;
    for ( unsigned int i = 0; i < Diffusion::NUMBER_VARIABLES; i++ ) {
        if ( d_isActive[i] )
            d_numberActive++;
        if ( d_isFrozen[i] )
            d_numberFrozen++;
    }

    AMP_INSIST( params->d_db->keyExists( "PrincipalVariable" ), "must specify PrincipalVariable" );
    int value = params->d_db->getScalar<int>( "PrincipalVariable" );
    AMP_INSIST( value >= 0, "can not specify negative PrincipalVariable" );
    AMP_INSIST( static_cast<unsigned int>( value ) < Diffusion::NUMBER_VARIABLES,
                "PrincipalVariable is too large" );
    d_PrincipalVariable = static_cast<unsigned int>( value );

    AMP_INSIST( d_isActive[d_PrincipalVariable], "must have Principal_Variable active" );
    d_diffNonlinElem->setPrincipalVariable( d_PrincipalVariable );

    resetFrozen( params );

    d_inpVariables.reset( new AMP::LinearAlgebra::MultiVariable( "InputVariables" ) );
    for ( unsigned int i = 0; i < Diffusion::NUMBER_VARIABLES; i++ ) {
        AMP::LinearAlgebra::Variable::shared_ptr dummyVar;
        d_inpVariables->add( dummyVar );
    } // end for i

    for ( unsigned int var = 0; var < Diffusion::NUMBER_VARIABLES; var++ ) {
        if ( d_isActive[var] ) {
            std::string name = activeVariables_db->getString( Diffusion::names[var] );
            AMP::LinearAlgebra::Variable::shared_ptr dummyVar(
                new AMP::LinearAlgebra::Variable( name ) );
            d_inpVariables->setVariable( var, dummyVar );
            if ( d_isFrozen[var] ) {
                d_inVec[var] = d_Frozen[var];
                if ( d_inVec[var] != nullptr )
                    AMP_ASSERT( d_inVec[var]->getUpdateStatus() ==
                                AMP::LinearAlgebra::VectorData::UpdateState::UNCHANGED );
            }
        } else {
            AMP::LinearAlgebra::Variable::shared_ptr dummyVar;
            d_inpVariables->add( dummyVar );
        }
    }

    d_outVariable.reset(
        new AMP::LinearAlgebra::Variable( params->d_db->getString( "OutputVariable" ) ) );

    init( params );
}


void DiffusionNonlinearFEOperator::preAssembly( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                                AMP::LinearAlgebra::Vector::shared_ptr r )
{
    // PROFILE_START("preAssembly",2);
    AMP_INSIST( ( u != nullptr ), "NULL Input Vector!" );
    AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
    AMP::LinearAlgebra::Vector::const_shared_ptr u_meshVec =
        u->constSelect( meshSelector, "u_mesh" );

    if ( d_iDebugPrintInfoLevel > 7 )
        AMP::pout << "DiffusionNonlinearFEOperator::preAssembly, entering" << std::endl;

    for ( unsigned int var = 0; var < Diffusion::NUMBER_VARIABLES; var++ ) {
        if ( d_isActive[var] ) {
            if ( d_isFrozen[var] ) {
                d_Frozen[var]->makeConsistent(
                    AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
                d_inVec[var] = d_Frozen[var];
            } else {
                AMP::LinearAlgebra::Variable::shared_ptr tvar = d_inpVariables->getVariable( var );
                d_inVec[var] = u_meshVec->constSubsetVectorForVariable( tvar );
                AMP_ASSERT( d_inVec[var] );
            }

            AMP_ASSERT( d_inVec[var] != nullptr );
            AMP_ASSERT( d_inVec[var]->getUpdateStatus() ==
                        AMP::LinearAlgebra::VectorData::UpdateState::UNCHANGED );
            if ( d_iDebugPrintInfoLevel > 5 )
                std::cout << "Max Value inside preAssembly: " << d_inVec[var]->max() << std::endl;
        }
    }

    d_outVec = subsetOutputVector( r );
    d_outVec->zero();

    d_transportModel->preNonlinearAssembly();

    if ( d_iDebugPrintInfoLevel > 7 )
        AMP::pout << "DiffusionNonlinearFEOperator::preAssembly, leaving" << std::endl;
    // PROFILE_STOP("preAssembly",2);
}


void DiffusionNonlinearFEOperator::postAssembly()
{
    if ( d_iDebugPrintInfoLevel > 7 )
        AMP::pout << "DiffusionNonlinearFEOperator::postAssembly, entering" << std::endl;

    d_transportModel->postNonlinearAssembly();
    d_outVec->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_ADD );

    if ( d_iDebugPrintInfoLevel > 7 )
        AMP::pout << "DiffusionNonlinearFEOperator::postAssembly, leaving" << std::endl;
}


void DiffusionNonlinearFEOperator::preElementOperation(


    const AMP::Mesh::MeshElement &elem )
{

    // PROFILE_START("preElementOperation",2);
    if ( d_iDebugPrintInfoLevel > 7 )
        AMP::pout << "DiffusionNonlinearFEOperator::preElementOperation, entering" << std::endl;

    std::vector<std::vector<double>> elementInputVectors( Diffusion::NUMBER_VARIABLES );

    d_currNodes = elem.getElements( AMP::Mesh::GeomType::Vertex );
    std::vector<AMP::Mesh::MeshElementID> ids( d_currNodes.size() );
    for ( size_t i = 0; i < d_currNodes.size(); i++ )
        ids[i] = d_currNodes[i].globalID();

    std::vector<size_t> dofs( d_currNodes.size() );
    for ( unsigned int var = 0; var < Diffusion::NUMBER_VARIABLES; var++ ) {
        if ( d_isActive[var] ) {
            AMP::Discretization::DOFManager::shared_ptr DOF = ( d_inVec[var] )->getDOFManager();
            DOF->getDOFs( ids, dofs );
            AMP_ASSERT( dofs.size() == d_currNodes.size() );
            elementInputVectors[var].resize( dofs.size() );
            ( d_inVec[var] )
                ->getValuesByGlobalID( dofs.size(), &dofs[0], &elementInputVectors[var][0] );
        }
    }

    d_elementOutputVector.resize( d_currNodes.size() );
    for ( unsigned int i = 0; i < d_currNodes.size(); i++ )
        d_elementOutputVector[i] = 0.0;

    d_diffNonlinElem->setElementVectors( elementInputVectors, d_elementOutputVector );

    d_diffNonlinElem->initializeForCurrentElement( d_currElemPtrs[d_currElemIdx],
                                                   d_transportModel );

    if ( d_iDebugPrintInfoLevel > 7 )
        AMP::pout << "DiffusionNonlinearFEOperator::preElementOperation, leaving" << std::endl;

    // PROFILE_STOP("preElementOperation",2);
}


void DiffusionNonlinearFEOperator::postElementOperation()
{
    // PROFILE_START("postElementOperation",2);
    if ( d_iDebugPrintInfoLevel > 7 )
        AMP::pout << "DiffusionNonlinearFEOperator::postElementOperation, entering" << std::endl;

    std::vector<AMP::Mesh::MeshElementID> ids( d_currNodes.size() );
    for ( size_t i = 0; i < d_currNodes.size(); i++ )
        ids[i] = d_currNodes[i].globalID();

    AMP::Discretization::DOFManager::shared_ptr DOF = d_outVec->getDOFManager();
    std::vector<size_t> dofs( d_currNodes.size() );
    DOF->getDOFs( ids, dofs );
    AMP_ASSERT( dofs.size() == d_currNodes.size() );

    d_outVec->addValuesByGlobalID( dofs.size(), &dofs[0], &d_elementOutputVector[0] );

    if ( d_iDebugPrintInfoLevel > 7 )
        AMP::pout << "DiffusionNonlinearFEOperator::postElementOperation, leaving" << std::endl;
    // PROFILE_STOP("postElementOperation",2);
}


void DiffusionNonlinearFEOperator::init(
    const std::shared_ptr<DiffusionNonlinearFEOperatorParameters> &params )
{
    if ( d_iDebugPrintInfoLevel > 7 )
        AMP::pout << "DiffusionNonlinearFEOperator::init, entering" << std::endl;

    (void) params;
    AMP::Mesh::MeshIterator el     = d_Mesh->getIterator( AMP::Mesh::GeomType::Volume, 0 );
    AMP::Mesh::MeshIterator end_el = el.end();

    for ( d_currElemIdx = 0; el != end_el; ++el, ++d_currElemIdx ) {
        d_currNodes = el->getElements( AMP::Mesh::GeomType::Vertex );
        d_diffNonlinElem->initializeForCurrentElement( d_currElemPtrs[d_currElemIdx],
                                                       d_transportModel );
        d_diffNonlinElem->initTransportModel();
    } // end for el
    d_currElemIdx = static_cast<unsigned int>( -1 );

    if ( d_iDebugPrintInfoLevel > 7 )
        AMP::pout << "DiffusionNonlinearFEOperator::init, leaving" << std::endl;
}


void DiffusionNonlinearFEOperator::reset( const std::shared_ptr<OperatorParameters> &params )
{
    std::shared_ptr<DiffusionNonlinearFEOperatorParameters> dnlparams_sp =
        std::dynamic_pointer_cast<DiffusionNonlinearFEOperatorParameters, OperatorParameters>(
            params );

    if ( d_PrincipalVariable == Diffusion::TEMPERATURE )
        d_inVec[d_PrincipalVariable] = dnlparams_sp->d_FrozenTemperature;
    if ( d_PrincipalVariable == Diffusion::CONCENTRATION )
        d_inVec[d_PrincipalVariable] = dnlparams_sp->d_FrozenConcentration;
    if ( d_PrincipalVariable == Diffusion::BURNUP )
        d_inVec[d_PrincipalVariable] = dnlparams_sp->d_FrozenBurnup;
    AMP_ASSERT( d_inVec[d_PrincipalVariable]->getUpdateStatus() ==
                AMP::LinearAlgebra::VectorData::UpdateState::UNCHANGED );

    resetFrozen( dnlparams_sp );
    for ( unsigned int var = 0; var < Diffusion::NUMBER_VARIABLES; var++ ) {
        if ( d_isActive[var] ) {
            if ( d_isFrozen[var] ) {
                d_inVec[var] = d_Frozen[var];
                AMP_ASSERT( d_inVec[var]->getUpdateStatus() ==
                            AMP::LinearAlgebra::VectorData::UpdateState::UNCHANGED );
            }
        }
    }
}


std::shared_ptr<OperatorParameters> DiffusionNonlinearFEOperator::getJacobianParameters(
    AMP::LinearAlgebra::Vector::const_shared_ptr u )
{
    //    AMP::LinearAlgebra::Vector::shared_ptr u  =
    //    std::const_pointer_cast<AMP::LinearAlgebra::Vector>(u_in);
    std::shared_ptr<AMP::Database> tmp_db( new AMP::Database( "Dummy" ) );
    AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
    auto u_meshVec = u->constSelect( meshSelector, "u_mesh" );

    // set up a database for the linear operator params
    tmp_db->putScalar( "name", "DiffusionLinearFEOperator" );
    tmp_db->putScalar( "InputVariable", Diffusion::names[d_PrincipalVariable] );
    tmp_db->putScalar( "OutputVariable", d_outVariable->getName() );
    tmp_db->putScalar( "FixedTemperature", d_isActive[Diffusion::TEMPERATURE] ? false : true );
    tmp_db->putScalar( "FixedConcentration", d_isActive[Diffusion::CONCENTRATION] ? false : true );
    tmp_db->putScalar( "FixedBurnup", d_isActive[Diffusion::BURNUP] ? false : true );

    // create the linear operator params
    std::shared_ptr<DiffusionLinearFEOperatorParameters> outParams(
        new DiffusionLinearFEOperatorParameters( tmp_db ) );

    // create the linear element object
    std::shared_ptr<AMP::Database> elem_db( new AMP::Database( "Dummy" ) );
    tmp_db->putScalar( "TransportAtGaussPoints", d_diffNonlinElem->getTransportAtGauss() );
    std::shared_ptr<ElementOperationParameters> eparams(
        new ElementOperationParameters( elem_db ) );
    std::shared_ptr<DiffusionLinearElement> linearElement( new DiffusionLinearElement( eparams ) );

    // add miscellaneous to output parameters
    outParams->d_transportModel = d_transportModel;
    outParams->d_elemOp         = linearElement;
    outParams->d_Mesh           = d_Mesh;

    // add variables to parameters
    if ( d_isActive[Diffusion::TEMPERATURE] ) {
        if ( d_isFrozen[Diffusion::TEMPERATURE] ) {
            outParams->d_temperature = d_Frozen[Diffusion::TEMPERATURE];
        } else {
            auto tvar        = d_inpVariables->getVariable( Diffusion::TEMPERATURE );
            auto temperature = u_meshVec->constSubsetVectorForVariable( tvar );
            outParams->d_temperature =
                std::const_pointer_cast<AMP::LinearAlgebra::Vector>( temperature );
            outParams->d_temperature->makeConsistent(
                AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
        }
    }

    if ( d_isActive[Diffusion::CONCENTRATION] ) {
        if ( d_isFrozen[Diffusion::CONCENTRATION] ) {
            outParams->d_concentration = d_Frozen[Diffusion::CONCENTRATION];
        } else {
            auto cvar          = d_inpVariables->getVariable( Diffusion::CONCENTRATION );
            auto concentration = u_meshVec->constSubsetVectorForVariable( cvar );
            outParams->d_concentration =
                std::const_pointer_cast<AMP::LinearAlgebra::Vector>( concentration );
            outParams->d_concentration->makeConsistent(
                AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
        }
    }

    if ( d_isActive[Diffusion::BURNUP] ) {
        if ( d_isFrozen[Diffusion::BURNUP] ) {
            outParams->d_burnup = d_Frozen[Diffusion::BURNUP];
        } else {
            auto bvar           = d_inpVariables->getVariable( Diffusion::BURNUP );
            auto burnup         = u_meshVec->constSubsetVectorForVariable( bvar );
            outParams->d_burnup = std::const_pointer_cast<AMP::LinearAlgebra::Vector>( burnup );
            outParams->d_burnup->makeConsistent(
                AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
        }
    }

    return outParams;
}


void DiffusionNonlinearFEOperator::resetFrozen(
    const std::shared_ptr<DiffusionNonlinearFEOperatorParameters> params )
{
    using namespace Diffusion;
    for ( size_t var = 0; var < Diffusion::NUMBER_VARIABLES; var++ ) {
        if ( d_isActive[var] and d_isFrozen[var] ) {
            if ( var == TEMPERATURE )
                d_Frozen[var] = params->d_FrozenTemperature;
            if ( var == CONCENTRATION )
                d_Frozen[var] = params->d_FrozenConcentration;
            if ( var == BURNUP )
                d_Frozen[var] = params->d_FrozenBurnup;
        }
    }
}


bool DiffusionNonlinearFEOperator::isValidInput( AMP::LinearAlgebra::Vector::shared_ptr &u )
{
    auto property                  = d_transportModel->getProperty();
    std::vector<std::string> names = property->get_arguments();
    size_t nnames                  = names.size();
    std::string argname;
    bool found = false;
    if ( d_PrincipalVariable == Diffusion::TEMPERATURE ) {
        for ( size_t i = 0; i < nnames; i++ ) {
            if ( names[i] == "temperature" ) {
                argname = "temperature";
                found   = true;
                break;
            }
        }
    }
    if ( d_PrincipalVariable == Diffusion::CONCENTRATION ) {
        for ( size_t i = 0; i < nnames; i++ ) {
            if ( names[i] == "concentration" ) {
                argname = "concentration";
                found   = true;
                break;
            }
        }
    }

    bool result = true;
    AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
    AMP::LinearAlgebra::Vector::shared_ptr u_meshVec = u->select( meshSelector, "u_mesh" );
    if ( found ) {
        AMP::LinearAlgebra::Vector::shared_ptr uinp = u_meshVec->subsetVectorForVariable(
            d_inpVariables->getVariable( d_PrincipalVariable ) );
        std::vector<double> vals( uinp->getLocalSize() );
        size_t nit = 0;
        for ( auto &elem : *uinp ) {
            vals[nit] = elem;
            nit++;
        }
        result = property->in_range( argname, vals );
    }

    return result;
}
} // namespace Operator
} // namespace AMP

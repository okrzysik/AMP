#include "AMP/operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/materials/Material.h"
#include "AMP/operators/ElementOperationParameters.h"
#include "AMP/operators/diffusion/DiffusionConstants.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperatorParameters.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorSelector.h"
#include "ProfilerApp.h"

#include <cstring>
#include <iostream>

namespace AMP::Operator {


static int getVariableID( const std::string &name )
{
    int id = -1;
    for ( size_t i = 0; i < Diffusion::NUMBER_VARIABLES; i++ ) {
        if ( name == Diffusion::names[i] )
            id = i;
    }
    if ( id == 100000 )
        AMP_ERROR( "Unknown variable: " + name );
    return id;
}


std::shared_ptr<AMP::LinearAlgebra::Variable>
DiffusionNonlinearFEOperator::createInputVariable( const std::string &name, int varId )
{
    if ( varId == -1 ) {
        return d_inpVariables->cloneVariable( name );
    } else {
        return d_inpVariables->getVariable( varId )->cloneVariable( name );
    }
}


std::shared_ptr<AMP::LinearAlgebra::Variable>
DiffusionNonlinearFEOperator::createOutputVariable( const std::string &name, int varId )
{
    (void) varId;
    return d_outVariable->cloneVariable( name );
}


std::shared_ptr<AMP::LinearAlgebra::Variable> DiffusionNonlinearFEOperator::getInputVariable()
{
    return d_inpVariables;
}


std::shared_ptr<AMP::LinearAlgebra::Variable> DiffusionNonlinearFEOperator::getOutputVariable()
{
    return d_outVariable;
}

unsigned int DiffusionNonlinearFEOperator::numberOfDOFMaps() { return 1; }


std::string DiffusionNonlinearFEOperator::getPrincipalVariable() { return d_PrincipalVariable; }


std::shared_ptr<DiffusionTransportModel> DiffusionNonlinearFEOperator::getTransportModel()
{
    return d_transportModel;
}


std::vector<AMP::LinearAlgebra::Vector::shared_ptr> DiffusionNonlinearFEOperator::getFrozen()
{
    return d_Frozen;
}


void DiffusionNonlinearFEOperator::setVector( const std::string &name,
                                              AMP::LinearAlgebra::Vector::shared_ptr frozenVec )
{
    bool set = false;
    for ( size_t id = 0; id < Diffusion::NUMBER_VARIABLES; id++ ) {
        if ( Diffusion::names[id] == name ) {
            AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
            auto meshSubsetVec = frozenVec->select( meshSelector, frozenVec->getName() );
            d_Frozen[id] =
                meshSubsetVec->subsetVectorForVariable( d_inpVariables->getVariable( id ) );
            d_Frozen[id]->makeConsistent(
                AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
            set = true;
        }
    }
    if ( !set )
        AMP_ERROR( "Variable " + name + " not found" );
}


DiffusionNonlinearFEOperator::DiffusionNonlinearFEOperator(
    std::shared_ptr<const DiffusionNonlinearFEOperatorParameters> params )
    : NonlinearFEOperator( params ), d_Frozen( Diffusion::NUMBER_VARIABLES )
{
    AMP_INSIST( params, "NULL parameter!" );

    d_diffNonlinElem = std::dynamic_pointer_cast<DiffusionNonlinearElement>( d_elemOp );

    AMP_INSIST( d_diffNonlinElem, "d_elemOp is not of type DiffusionNonlinearElement" );

    d_transportModel = params->d_transportModel;

    d_isActive.resize( Diffusion::NUMBER_VARIABLES );
    d_isFrozen.resize( Diffusion::NUMBER_VARIABLES );
    d_inVec.resize( Diffusion::NUMBER_VARIABLES );

    auto activeVariables_db = params->d_db->getDatabase( "ActiveInputVariables" );

    for ( size_t var = 0; var < Diffusion::NUMBER_VARIABLES; var++ ) {
        auto namespec   = activeVariables_db->getWithDefault<std::string>( Diffusion::names[var],
                                                                         "not_specified" );
        bool isthere    = namespec != "not_specified";
        d_isActive[var] = isthere;
        std::string frozenName = "Freeze" + Diffusion::names[var];
        d_isFrozen[var]        = params->d_db->getWithDefault<bool>( frozenName, false );
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

    d_PrincipalVariable     = params->d_db->getString( "PrincipalVariable" );
    int principalVariableId = getVariableID( d_PrincipalVariable );
    AMP_INSIST( d_isActive[principalVariableId], "must have Principal_Variable active" );
    d_diffNonlinElem->setPrincipalVariable( d_PrincipalVariable );

    resetFrozen( params );

    d_inpVariables.reset( new AMP::LinearAlgebra::MultiVariable( "InputVariables" ) );
    for ( unsigned int i = 0; i < Diffusion::NUMBER_VARIABLES; i++ ) {
        std::shared_ptr<AMP::LinearAlgebra::Variable> dummyVar;
        d_inpVariables->add( dummyVar );
    } // end for i

    for ( unsigned int var = 0; var < Diffusion::NUMBER_VARIABLES; var++ ) {
        if ( d_isActive[var] ) {
            std::string name = activeVariables_db->getString( Diffusion::names[var] );
            auto dummyVar    = std::make_shared<AMP::LinearAlgebra::Variable>( name );
            d_inpVariables->setVariable( var, dummyVar );
            if ( d_isFrozen[var] ) {
                d_inVec[var] = d_Frozen[var];
                if ( d_inVec[var] )
                    AMP_ASSERT( d_inVec[var]->getUpdateStatus() ==
                                AMP::LinearAlgebra::VectorData::UpdateState::UNCHANGED );
            }
        } else {
            std::shared_ptr<AMP::LinearAlgebra::Variable> dummyVar;
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
    AMP_INSIST( u, "NULL Input Vector!" );
    AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
    auto u_meshVec = u->select( meshSelector, "u_mesh" );

    if ( d_iDebugPrintInfoLevel > 7 )
        AMP::pout << "DiffusionNonlinearFEOperator::preAssembly, entering" << std::endl;

    for ( unsigned int var = 0; var < Diffusion::NUMBER_VARIABLES; var++ ) {
        if ( d_isActive[var] ) {
            if ( d_isFrozen[var] ) {
                d_Frozen[var]->makeConsistent(
                    AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
                d_inVec[var] = d_Frozen[var];
            } else {
                auto tvar    = d_inpVariables->getVariable( var );
                d_inVec[var] = u_meshVec->subsetVectorForVariable( tvar );
                if ( !d_inVec[var] )
                    AMP_ERROR( "Unable to subset for " + tvar->getName() );
            }

            AMP_ASSERT( d_inVec[var] );
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


void DiffusionNonlinearFEOperator::preElementOperation( const AMP::Mesh::MeshElement &elem )
{
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
            auto DOF = ( d_inVec[var] )->getDOFManager();
            DOF->getDOFs( ids, dofs );
            AMP_ASSERT( dofs.size() == d_currNodes.size() );
            elementInputVectors[var].resize( dofs.size() );
            d_inVec[var]->getValuesByGlobalID(
                dofs.size(), &dofs[0], &elementInputVectors[var][0] );
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
}


void DiffusionNonlinearFEOperator::postElementOperation()
{
    if ( d_iDebugPrintInfoLevel > 7 )
        AMP::pout << "DiffusionNonlinearFEOperator::postElementOperation, entering" << std::endl;

    std::vector<AMP::Mesh::MeshElementID> ids( d_currNodes.size() );
    for ( size_t i = 0; i < d_currNodes.size(); i++ )
        ids[i] = d_currNodes[i].globalID();

    auto DOF = d_outVec->getDOFManager();
    std::vector<size_t> dofs( d_currNodes.size() );
    DOF->getDOFs( ids, dofs );
    AMP_ASSERT( dofs.size() == d_currNodes.size() );

    d_outVec->addValuesByGlobalID( dofs.size(), &dofs[0], &d_elementOutputVector[0] );

    if ( d_iDebugPrintInfoLevel > 7 )
        AMP::pout << "DiffusionNonlinearFEOperator::postElementOperation, leaving" << std::endl;
}


void DiffusionNonlinearFEOperator::init(
    std::shared_ptr<const DiffusionNonlinearFEOperatorParameters> )
{
    if ( d_iDebugPrintInfoLevel > 7 )
        AMP::pout << "DiffusionNonlinearFEOperator::init, entering" << std::endl;

    auto el     = d_Mesh->getIterator( AMP::Mesh::GeomType::Cell, 0 );
    auto end_el = el.end();
    for ( d_currElemIdx = 0; el != end_el; ++el, ++d_currElemIdx ) {
        d_currNodes = el->getElements( AMP::Mesh::GeomType::Vertex );
        d_diffNonlinElem->initializeForCurrentElement( d_currElemPtrs[d_currElemIdx],
                                                       d_transportModel );
        d_diffNonlinElem->initTransportModel();
    }
    d_currElemIdx = static_cast<unsigned int>( -1 );

    if ( d_iDebugPrintInfoLevel > 7 )
        AMP::pout << "DiffusionNonlinearFEOperator::init, leaving" << std::endl;
}


void DiffusionNonlinearFEOperator::reset( std::shared_ptr<const OperatorParameters> inParams )
{
    auto params =
        std::dynamic_pointer_cast<const DiffusionNonlinearFEOperatorParameters>( inParams );

    for ( auto [name, vec] : params->d_FrozenVecs ) {
        if ( name == d_PrincipalVariable ) {
            int id      = getVariableID( d_PrincipalVariable );
            d_inVec[id] = vec;
            AMP_ASSERT( vec->getUpdateStatus() ==
                        AMP::LinearAlgebra::VectorData::UpdateState::UNCHANGED );
        }
    }

    resetFrozen( params );
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
    auto db = std::make_shared<AMP::Database>( "Dummy" );
    AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
    auto u_meshVec = u->select( meshSelector, "u_mesh" );

    // set up a database for the linear operator params
    db->putScalar( "name", "DiffusionLinearFEOperator" );
    db->putScalar( "InputVariable", d_PrincipalVariable );
    db->putScalar( "OutputVariable", d_outVariable->getName() );
    db->putScalar( "FixedTemperature", d_isActive[Diffusion::TEMPERATURE] ? false : true );
    db->putScalar( "FixedConcentration", d_isActive[Diffusion::CONCENTRATION] ? false : true );
    db->putScalar( "FixedBurnup", d_isActive[Diffusion::BURNUP] ? false : true );

    // create the linear operator params
    auto outParams = std::make_shared<DiffusionLinearFEOperatorParameters>( db );

    // create the linear element object
    auto elem_db = std::make_shared<AMP::Database>( "Dummy" );
    db->putScalar( "TransportAtGaussPoints", d_diffNonlinElem->getTransportAtGauss() );
    auto eparams       = std::make_shared<ElementOperationParameters>( elem_db );
    auto linearElement = std::make_shared<DiffusionLinearElement>( eparams );

    // add miscellaneous to output parameters
    auto createDOFManager = []( std::shared_ptr<AMP::Mesh::Mesh> mesh ) {
        return AMP::Discretization::simpleDOFManager::create(
            mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    };
    outParams->d_transportModel = d_transportModel;
    outParams->d_elemOp         = linearElement;
    outParams->d_Mesh           = d_Mesh;
    outParams->d_inDofMap       = createDOFManager( d_Mesh );
    outParams->d_outDofMap      = outParams->d_inDofMap;

    // add variables to parameters
    for ( size_t i = 0; i < d_isActive.size(); i++ ) {
        if ( d_isActive[i] ) {
            auto name = Diffusion::names[i];
            if ( d_isFrozen[i] ) {
                outParams->d_inputVecs[name] = d_Frozen[i];
            } else {
                auto var = d_inpVariables->getVariable( i );
                auto vec = std::const_pointer_cast<AMP::LinearAlgebra::Vector>(
                    u_meshVec->subsetVectorForVariable( var ) );
                vec->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
                outParams->d_inputVecs[name] = vec;
            }
        }
    }

    return outParams;
}


void DiffusionNonlinearFEOperator::resetFrozen(
    std::shared_ptr<const DiffusionNonlinearFEOperatorParameters> params )
{
    for ( size_t var = 0; var < Diffusion::NUMBER_VARIABLES; var++ ) {
        if ( d_isActive[var] && d_isFrozen[var] ) {
            d_Frozen[var].reset();
            auto it = params->d_FrozenVecs.find( Diffusion::names[var] );
            if ( it != params->d_FrozenVecs.end() )
                d_Frozen[var] = it->second;
        }
    }
}


bool DiffusionNonlinearFEOperator::isValidInput( AMP::LinearAlgebra::Vector::const_shared_ptr u )
{
    auto property = d_transportModel->getProperty();
    auto names    = property->get_arguments();
    size_t nnames = names.size();
    std::string argname;
    bool found = false;
    for ( size_t i = 0; i < nnames; i++ ) {
        if ( names[i] == d_PrincipalVariable ) {
            argname = d_PrincipalVariable;
            found   = true;
            break;
        }
    }

    bool result = true;
    AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
    auto u_meshVec = u->select( meshSelector, "u_mesh" );
    if ( found ) {
        auto uinp = u_meshVec->subsetVectorForVariable( d_PrincipalVariable );
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


std::vector<std::string> DiffusionNonlinearFEOperator::getNonPrincipalVariableIds()
{
    std::vector<std::string> vars;
    for ( size_t i = 0; i < Diffusion::NUMBER_VARIABLES; i++ ) {
        auto var = Diffusion::names[i];
        if ( d_isActive[i] && var != d_PrincipalVariable )
            vars.push_back( var );
    }
    return vars;
}


} // namespace AMP::Operator

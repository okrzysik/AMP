#include "solvers/trilinos/muelu/TrilinosMueLuSolver.h"

DISABLE_WARNINGS
// Teuchos
#include <Teuchos_RCP.hpp>

// Xpetra
#include <Xpetra_EpetraVector.hpp>
// MueLu
#include "MueLu.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_Ifpack2Smoother.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_RAPFactory.hpp"
#include <MueLu_CreateEpetraPreconditioner.hpp>

ENABLE_WARNINGS

#include "ProfilerApp.h"
#include "matrices/Matrix.h"
#include "matrices/trilinos/EpetraMatrix.h"
#include "operators/LinearOperator.h"
#include "utils/Utilities.h"
#include "vectors/DataChangeFirer.h"
#include "vectors/trilinos/epetra/EpetraVector.h"

namespace AMP {
namespace Solver {


    

/****************************************************************
* Constructors / Destructor                                     *
****************************************************************/
TrilinosMueLuSolver::TrilinosMueLuSolver()
{
    d_bCreationPhase = true;
}
TrilinosMueLuSolver::TrilinosMueLuSolver( AMP::shared_ptr<SolverStrategyParameters> parameters )
    : SolverStrategy( parameters )
{
    AMP_ASSERT( parameters.get() != nullptr );
    initialize( parameters );
}
TrilinosMueLuSolver::~TrilinosMueLuSolver()
{
    d_mueluSolver.reset();
    d_matrix.reset(); // Need to keep a copy of the matrix alive until after the solver is destroyed
}

void TrilinosMueLuSolver::initialize( AMP::shared_ptr<SolverStrategyParameters> const parameters )
{
    getFromInput( parameters->d_db );

    if ( d_pOperator.get() != nullptr ) {

        registerOperator( d_pOperator );

    }
}

Teuchos::RCP<MueLu::TentativePFactory<SC, LO, GO, NO>> TrilinosMueLuSolver::getTentativePFactory( void )
{
    return Teuchos::rcp( new MueLu::TentativePFactory<SC, LO, GO, NO>() );
}

Teuchos::RCP<MueLu::SaPFactory<SC, LO, GO, NO>> TrilinosMueLuSolver::getSaPFactory( void )
{
    return Teuchos::rcp( new MueLu::SaPFactory<SC, LO, GO, NO>() );
}
    
Teuchos::RCP<MueLu::TransPFactory<SC, LO, GO, NO>> TrilinosMueLuSolver::getRFactory( void )
{
    return Teuchos::rcp( new MueLu::TransPFactory<SC, LO, GO, NO>() );
}

Teuchos::RCP<MueLu::SmootherFactory<SC, LO, GO, NO>> TrilinosMueLuSolver::getCoarseSolverFactory( void )
{
    auto coarseSolverPrototype = Teuchos::rcp( new MueLu::DirectSolver<SC, LO, GO, NO>() );
    return Teuchos::rcp( new MueLu::SmootherFactory<SC, LO, GO, NO>(coarseSolverPrototype, Teuchos::null) );
}

Teuchos::RCP<MueLu::SmootherFactory<SC, LO, GO, NO>> TrilinosMueLuSolver::getSmootherFactory( void )
{
    std::string ifpackType; 

    Teuchos::RCP<MueLu::SmootherPrototype<SC, LO, GO, NO>> smootherPrototype;
    
    if ( d_smoother_type == "TrilinosSmoother" ) {
        ifpackType = "RELAXATION";
        smootherPrototype = Teuchos::rcp(new MueLu::TrilinosSmoother<SC, LO, GO, NO>(ifpackType, d_MueLuParameterList));
    } else if ( d_smoother_type == "IfpackSmoother" ) {
        ifpackType = "point relaxation stand-alone"; 
        smootherPrototype = Teuchos::rcp(new MueLu::IfpackSmoother<NO>(ifpackType, d_MueLuParameterList));
    } else if ( d_smoother_type == "Ifpack2Smoother" ) {
        ifpackType = "point relaxation stand-alone"; 
        smootherPrototype = Teuchos::rcp(new MueLu::Ifpack2Smoother<SC, LO, GO, NO>(ifpackType, d_MueLuParameterList));
    }

    return Teuchos::rcp(new MueLu::SmootherFactory<SC, LO, GO, NO>(smootherPrototype));
}

Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> TrilinosMueLuSolver::getXpetraMatrix( AMP::shared_ptr<AMP::Operator::LinearOperator> & op )
{
    // wrap in a Xpetra matrix
    auto ampMatrix       = op->getMatrix();
    auto epetraMatrix    = AMP::dynamic_pointer_cast<AMP::LinearAlgebra::EpetraMatrix> ( AMP::LinearAlgebra::EpetraMatrix::createView( ampMatrix ) );
    auto epA             = Teuchos::rcpFromRef( epetraMatrix->getEpetra_CrsMatrix() );
    Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> exA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix( epA ) );
    auto crsWrapMat      = Teuchos::rcp( new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>( exA ) ); 
    auto xA    = Teuchos::rcp_dynamic_cast<Xpetra::Matrix<SC, LO, GO, NO>> ( crsWrapMat );

    return xA;
}

void TrilinosMueLuSolver::buildHierarchyFromDefaults( void )
{

    // Transfer operators
    auto TentativePFact       = Teuchos::rcp( new MueLu::TentativePFactory<SC, LO, GO, NO>() );
    auto SaPFact              = Teuchos::rcp( new MueLu::SaPFactory<SC, LO, GO, NO>() );
    auto RFact                = Teuchos::rcp( new MueLu::TransPFactory<SC, LO, GO, NO>());
    // coarsest solver.
    auto coarseSolverPrototype = Teuchos::rcp( new MueLu::DirectSolver<SC, LO, GO, NO>() );
    auto   coarseSolverFact    = Teuchos::rcp( new MueLu::SmootherFactory<SC, LO, GO, NO>(coarseSolverPrototype, Teuchos::null) );
                
    d_factoryManager.SetFactory("Ptent", TentativePFact);
    d_factoryManager.SetFactory("P",     SaPFact);
    d_factoryManager.SetFactory("R",     RFact);
    //    d_factoryManager.SetFactory("Smoother", Teuchos::null);      //skips smoother setup
    d_factoryManager.SetFactory("CoarseSolver", coarseSolverFact);
                
    // extract the Xpetra matrix from AMP
    auto linearOperator  = AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator> ( d_pOperator );           
    auto fineLevelA = getXpetraMatrix( linearOperator );
                
    d_mueluHierarchy     = Teuchos::rcp( new MueLu::Hierarchy<SC, LO, GO, NO >() );
    d_mueluHierarchy->SetDefaultVerbLevel(MueLu::Medium);
    auto finestMGLevel   = d_mueluHierarchy->GetLevel();
    finestMGLevel->Set( "A", fineLevelA );
                
    const int startLevel = 0;
    const int maxLevels  = 10;
    d_mueluHierarchy->Setup( d_factoryManager, startLevel, maxLevels );

}

void TrilinosMueLuSolver::buildHierarchyByLevel( void )
{

    d_mueluHierarchyManager = Teuchos::rcp(new MueLu::ParameterListInterpreter<SC, LO, GO, NO >( d_MueLuParameterList ) );

#if 0    
    d_mueluHierarchy     = Teuchos::rcp( new MueLu::Hierarchy<SC, LO, GO, NO >() );
#else
    d_mueluHierarchy     = d_mueluHierarchyManager->CreateHierarchy();
#endif

    d_mueluHierarchy->SetDefaultVerbLevel(MueLu::Medium);
    auto finestMGLevel   = d_mueluHierarchy->GetLevel();
    // extract the Xpetra matrix from AMP
    auto linearOperator  = AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator> ( d_pOperator );           
    auto fineLevelA = getXpetraMatrix( linearOperator );            
    finestMGLevel->Set( "A", fineLevelA );

    d_levelFactoryManager.resize( d_maxLevels );

#if 1
    for ( size_t i=0u; i<d_maxLevels; ++i ) {

        d_levelFactoryManager[i] = Teuchos::rcp( new MueLu::FactoryManager<SC,LO,GO,NO> );
        d_levelFactoryManager[i]->SetFactory("Ptent", getTentativePFactory() );
        d_levelFactoryManager[i]->SetFactory("P",     getSaPFactory() );
        d_levelFactoryManager[i]->SetFactory("R",     getRFactory() );
    }

    // seeing if smoothers can be initialized for each level
    for ( size_t i=0u; i<d_maxLevels; ++i ) {
        
        auto finerLevelManager = (i==0u) ? Teuchos::null : d_levelFactoryManager[i-1];
        auto currentLevelManager = d_levelFactoryManager[i];
        auto coarserLevelManager = ( i < d_maxLevels-1 ) ? d_levelFactoryManager[i+1] : Teuchos::null;

        // setup a level
        bool bIsLastLevel = d_mueluHierarchy->Setup( i,
                                                     finerLevelManager,
                                                     currentLevelManager,
                                                     coarserLevelManager );

        auto currentLevel = d_mueluHierarchy->GetLevel( i );

        // keep some of the factories
        auto PtentFactory = currentLevelManager->GetFactory("Ptent");
        currentLevel->Keep("Ptent",  PtentFactory.get() );

        auto PFactory = currentLevelManager->GetFactory("P");
        currentLevel->Keep("P",  PFactory.get() );

        auto RFactory = currentLevelManager->GetFactory("R");
        currentLevel->Keep("R",  RFactory.get() );

        currentLevelManager->SetFactory("Smoother", getSmootherFactory() );      
        currentLevelManager->SetFactory("CoarseSolver", getCoarseSolverFactory() );

        // setup a level
        bIsLastLevel = d_mueluHierarchy->Setup( i,
                                                finerLevelManager,
                                                currentLevelManager,
                                                coarserLevelManager );
        
        if ( bIsLastLevel ) break;
        
    }
#else
    for ( size_t i=0u; i<d_maxLevels; ++i ) {

        d_levelFactoryManager[i] = Teuchos::rcp( new MueLu::FactoryManager<SC,LO,GO,NO> );
        d_levelFactoryManager[i]->SetFactory("Ptent", getTentativePFactory() );
        d_levelFactoryManager[i]->SetFactory("P",     getSaPFactory() );
        d_levelFactoryManager[i]->SetFactory("R",     getRFactory() );
        // explicitly setting the factory to null appears to be the wrong thing to do
        // MueLu decides no smoothing is required on the level
        d_levelFactoryManager[i]->SetFactory("Smoother", getSmootherFactory() );      
        d_levelFactoryManager[i]->SetFactory("CoarseSolver", getCoarseSolverFactory() );
    }

    for ( size_t i=0u; i<d_maxLevels; ++i ) {
        
        auto finerLevelManager = (i==0u) ? Teuchos::null : d_levelFactoryManager[i-1];
        auto currentLevelManager = d_levelFactoryManager[i];
        auto coarserLevelManager = ( i < d_maxLevels-1 ) ? d_levelFactoryManager[i+1] : Teuchos::null;

        bool bIsLastLevel = d_mueluHierarchy->Setup( i,
                                                     finerLevelManager,
                                                     currentLevelManager,
                                                     coarserLevelManager );

        if ( bIsLastLevel ) break;
        
    }
#endif
}

void TrilinosMueLuSolver::getFromInput( const AMP::shared_ptr<AMP::Database> &db )
{
    d_bRobustMode = db->getBoolWithDefault( "ROBUST_MODE", false );
    d_bUseEpetra  = db->getBoolWithDefault( "USE_EPETRA", true );

#if 1
    d_build_hierarchy = db->getBoolWithDefault( "build_hierarchy", true );
    d_build_hierarchy_from_defaults = db->getBoolWithDefault( "build_hierarchy_from_defaults", false );
#else
    d_build_hierarchy = db->getBoolWithDefault( "build_hierarchy", false );
    d_build_hierarchy_from_defaults = db->getBoolWithDefault( "build_hierarchy_from_defaults", true );
#endif
    // general parameters
    d_MueLuParameterList.set( "verbosity", db->getStringWithDefault("verbosity", "medium") );
    d_MueLuParameterList.set( "problem: type", db->getStringWithDefault("problem_type", "unknown"));
    d_MueLuParameterList.set( "number of equations", db->getIntegerWithDefault("number_of_equations", 1));
    d_maxLevels = db->getIntegerWithDefault("max_levels", 10);
    d_MueLuParameterList.set( "max levels", (int) d_maxLevels);
    d_MueLuParameterList.set( "cycle type", db->getStringWithDefault("cycle_type", "V"));
    d_MueLuParameterList.set( "problem: symmetric", db->getBoolWithDefault("problem_symmetric", false));
    d_MueLuParameterList.set( "xml parameter file", db->getStringWithDefault("xml_parameter_file", ""));
    
    // smoothing and coarse solver options
    d_smoother_type = db->getStringWithDefault("smoother_type", "IfpackSmoother");
    
    d_MueLuParameterList.set( "smoother: pre or post", db->getStringWithDefault("smoother_pre_or_post", "both"));
    d_MueLuParameterList.set( "smoother: pre type", db->getStringWithDefault("smoother_pre_type", "RELAXATION"));
    d_MueLuParameterList.set( "smoother: post type", db->getStringWithDefault("smoother_post_type", "RELAXATION"));
    d_MueLuParameterList.set( "smoother: overlap", db->getIntegerWithDefault("smoother_overlap", 0));
    d_MueLuParameterList.set( "smoother: pre overlap", db->getIntegerWithDefault("smoother_pre_overlap", 0));
    d_MueLuParameterList.set( "smoother: post overlap", db->getIntegerWithDefault("smoother_post_overlap", 0));

#if 0    
    d_MueLuParameterList.set( "relaxation: type", db->getStringWithDefault("relaxation_type", "Gauss-Seidel") );
    d_MueLuParameterList.set( "relaxation: sweeps", db->getIntegerWithDefault("relaxation_sweeps", 1) );
    d_MueLuParameterList.set( "relaxation: damping factor", db->getDoubleWithDefault("relaxation_damping_factor", 1.0) );

    if ( db->keyExists("relaxation_zero_starting_solution") ) {     
        d_MueLuParameterList.set( "relaxation: zero starting solution", db->getBool("relaxation_zero_starting_solution") );
    }

    if ( db->keyExists("relaxation_backward_mode") ) { 
        d_MueLuParameterList.set( "relaxation: backward mode", db->getBool("relaxation_backward_mode"));
    }
    
    if ( db->keyExists("relaxation_use_l1") ) { 
        d_MueLuParameterList.set( "relaxation: use l1", db->getBool("relaxation_use_l1") );
    }
    
    if ( db->keyExists("relaxation_l1_eta") ) { 
        d_MueLuParameterList.set( "relaxation: l1 eta", db->getDouble("relaxation_l1_eta") );
    }
    
    if ( db->keyExists("relaxation_min_diagonal_value") ) { 
        d_MueLuParameterList.set( "relaxation: min diagonal value", db->getDouble("relaxation_min_diagonal_value") );
    }

    if ( db->keyExists("relaxation_fix_tiny_diagonal_entries") ) { 
        d_MueLuParameterList.set( "relaxation: fix tiny diagonal entries", db->getBool("relaxation_fix_tiny_diagonal_entries") );
    }

    if ( db->keyExists("relaxation_check_diagonal_entries") ) { 
        d_MueLuParameterList.set( "relaxation: check diagonal entries", db->getBool("relaxation_check_diagonal_entries") );
    }

    if ( db->keyExists("relaxation_local_smoothing_indices") ) { 
        d_MueLuParameterList.set( "relaxation: local smoothing indices", db->getIntegerWithDefault("relaxation_local_smoothing_indices") );
    }
#endif
    
    d_MueLuParameterList.set( "coarse: max size", db->getIntegerWithDefault("coarse_max_size", 2000));
    d_MueLuParameterList.set( "coarse: type", db->getStringWithDefault("coarse_type", "SuperLU"));
    d_MueLuParameterList.set( "coarse: overlap", db->getIntegerWithDefault("coarse_overlap", 0));

    // aggregation options: incomplete list
    d_MueLuParameterList.set( "aggregation: type", db->getStringWithDefault("aggregation_type", "uncoupled"));
    d_MueLuParameterList.set( "aggregation: ordering", db->getStringWithDefault("aggregation_ordering", "natural"));
    d_MueLuParameterList.set( "aggregation: drop scheme", db->getStringWithDefault("aggregation_drop_scheme", "classical"));
    d_MueLuParameterList.set( "aggregation: drop tol", db->getDoubleWithDefault("aggregation_drop_tol", 0.0));
    d_MueLuParameterList.set( "aggregation: min agg size", db->getIntegerWithDefault("aggregation_min_agg_size", 2));
    d_MueLuParameterList.set( "aggregation: max agg size", db->getIntegerWithDefault("aggregation_max_agg_size", -1));
    d_MueLuParameterList.set( "aggregation: brick x size", db->getIntegerWithDefault("aggregation_brick_x_size", 2));
    d_MueLuParameterList.set( "aggregation: brick y size", db->getIntegerWithDefault("aggregation_brick_y_size", 2));
    d_MueLuParameterList.set( "aggregation: brick z size", db->getIntegerWithDefault("aggregation_brick_z_size", 2));

    // rebalance options
    d_MueLuParameterList.set( "repartition: enable", db->getBoolWithDefault("repartition_enable", false));
    d_MueLuParameterList.set( "repartition: partitioner", db->getStringWithDefault("repartition_partitioner", "zoltan2"));
    d_MueLuParameterList.set( "repartition: start level", db->getIntegerWithDefault("repartition_start_level",2));
    d_MueLuParameterList.set( "repartition: min rows per proc", db->getIntegerWithDefault("repartition_min_rows_per_proc", 800));
    d_MueLuParameterList.set( "repartition: max imbalance", db->getDoubleWithDefault("repartition_max_imbalance", 1.2));
    d_MueLuParameterList.set( "repartition: remap parts", db->getBoolWithDefault("repartition_remap_parts", true));
    d_MueLuParameterList.set( "repartition: rebalance P and R", db->getBoolWithDefault("repartition_P_and_R",false));

    // mg algorithm options, incomplete list
    d_MueLuParameterList.set( "multigrid algorithm", db->getStringWithDefault("multigrid_algorithm", "sa"));
    d_MueLuParameterList.set( "sa: damping factor", db->getDoubleWithDefault("sa_damping_factor", 1.33));

    // mg reuse options: none at present

    // miscellaneous options 
    d_MueLuParameterList.set( "print initial parameters", db->getBoolWithDefault("print_initial_parameters", true));
    d_MueLuParameterList.set( "print unused parameters", db->getBoolWithDefault("print_unused_parameters", true));

    d_MueLuParameterList.set( "transpose: use implicit", db->getBoolWithDefault("transpose_use_implicit", false));

}

void TrilinosMueLuSolver::registerOperator( const AMP::shared_ptr<AMP::Operator::Operator> op )
{
    d_pOperator = op;
    AMP_INSIST( d_pOperator.get() != nullptr,
                "ERROR: TrilinosMueLuSolver::initialize() operator cannot be NULL" );

    if ( d_bUseEpetra ) {

        if( d_build_hierarchy ) {

            if( d_build_hierarchy_from_defaults ) {

                buildHierarchyFromDefaults();

            } else {

                buildHierarchyByLevel();
            }

        } else {
            
            auto linearOperator =
                AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>( d_pOperator );
            AMP_INSIST( linearOperator.get() != nullptr, "linearOperator cannot be NULL" );
            
            d_matrix = AMP::dynamic_pointer_cast<AMP::LinearAlgebra::EpetraMatrix>(
                                                                                   linearOperator->getMatrix() );
            AMP_INSIST( d_matrix.get() != nullptr, "d_matrix cannot be NULL" );
            
            // MueLu expects a Teuchos ref pointer
            Teuchos::RCP<Epetra_CrsMatrix> fineLevelMatrixPtr = Teuchos::rcpFromRef(d_matrix->getEpetra_CrsMatrix());
            auto mueluRCP = MueLu::CreateEpetraPreconditioner(fineLevelMatrixPtr, 
                                                              d_MueLuParameterList);
            
            d_mueluSolver.reset( mueluRCP.get() );
            mueluRCP.release();
        }
        
    } else {
        AMP_ERROR("Only Epetra interface currently supported");
    }
    d_bCreationPhase = false;
}


void TrilinosMueLuSolver::resetOperator(
    const AMP::shared_ptr<AMP::Operator::OperatorParameters> params )
{
    PROFILE_START( "resetOperator" );
    AMP_INSIST( ( d_pOperator.get() != nullptr ),
                "ERROR: TrilinosMueLuSolver::resetOperator() operator cannot be NULL" );
    d_pOperator->reset( params );
    reset( AMP::shared_ptr<SolverStrategyParameters>() );
    PROFILE_STOP( "resetOperator" );
}


void TrilinosMueLuSolver::reset( AMP::shared_ptr<SolverStrategyParameters> )
{
    PROFILE_START( "reset" );
    d_mueluSolver.reset();
    d_matrix.reset(); // Need to keep a copy of the matrix alive until after the solver is destroyed
    registerOperator( d_pOperator );
    PROFILE_STOP( "reset" );
}

void TrilinosMueLuSolver::solveWithHierarchy( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                                              AMP::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    PROFILE_START( "solveWithHierarchy" );

    if ( d_bUseEpetra ) {

        
        // These functions throw exceptions if this cannot be performed.
        Epetra_Vector &fVec = ( AMP::LinearAlgebra::EpetraVector::view( AMP::const_pointer_cast<AMP::LinearAlgebra::Vector>(f) ) )
            ->castTo<AMP::LinearAlgebra::EpetraVector>()
            .getEpetra_Vector();
        Epetra_Vector &uVec = ( AMP::LinearAlgebra::EpetraVector::view( u ) )
            ->castTo<AMP::LinearAlgebra::EpetraVector>()
            .getEpetra_Vector();

        auto rcp_u = Teuchos::rcpFromRef( uVec );
        auto rcp_f = Teuchos::rcpFromRef( fVec );

        // Epetra -> Xpetra
        auto xu = Teuchos::rcp( new Xpetra::EpetraVectorT <int, NO> ( rcp_u ) );
        auto xf = Teuchos::rcp( new Xpetra::EpetraVectorT <int, NO> ( rcp_f ) );

        d_mueluHierarchy->Iterate( *xf, *xu, d_iMaxIterations );
        
    } else {
        AMP_ERROR("Only Epetra interface currently supported");
    }        
    

    PROFILE_STOP( "solveWithHierarchy" );
    
}

void TrilinosMueLuSolver::solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                                 AMP::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    PROFILE_START( "solve" );

    AMP_ASSERT( f != nullptr );
    AMP_ASSERT( u != nullptr );

    // in this case we make the assumption we can access a EpetraMat for now
    AMP_INSIST( d_pOperator.get() != nullptr,
                "ERROR: TrilinosMueLuSolver::solve() operator cannot be NULL" );
    
    if ( d_bUseZeroInitialGuess ) {
        u->zero();
    }

    
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> r;
    
    bool computeResidual = false;
    if ( d_bRobustMode || ( d_iDebugPrintInfoLevel > 1 ) ) {
        computeResidual = true;
    }
    
    double initialResNorm = 0., finalResNorm = 0.;
    
    if ( computeResidual ) {
        r = f->cloneVector();
        d_pOperator->residual( f, u, r );
        initialResNorm = r->L2Norm();
        
        if ( d_iDebugPrintInfoLevel > 1 ) {
            AMP::pout << "TrilinosMueLuSolver::solve(), L2 norm of residual before solve "
                      << std::setprecision( 15 ) << initialResNorm << std::endl;
        }
    }
    
    if ( d_iDebugPrintInfoLevel > 2 ) {
        double solution_norm = u->L2Norm();
        AMP::pout << "TrilinosMueLuSolver : before solve solution norm: " << std::setprecision( 15 )
                  << solution_norm << std::endl;
    }
    
    if( d_build_hierarchy ) {

        solveWithHierarchy( f, u );

    } else {
        
        if ( d_bCreationPhase ) {
            if ( d_bUseEpetra ) {
                // MueLu expects a Teuchos ref pointer
                Teuchos::RCP<Epetra_CrsMatrix> fineLevelMatrixPtr = Teuchos::rcpFromRef(d_matrix->getEpetra_CrsMatrix());
                auto mueluRCP = MueLu::CreateEpetraPreconditioner(fineLevelMatrixPtr, 
                                                                  d_MueLuParameterList);
                d_mueluSolver.reset( mueluRCP.get() );
                mueluRCP.release();
            } else {
                AMP_ERROR("Only Epetra interface currently supported");
            }
            d_bCreationPhase = false;
        }
        
        if ( d_bUseEpetra ) {
            // These functions throw exceptions if this cannot be performed.
            const Epetra_Vector &fVec = ( AMP::LinearAlgebra::EpetraVector::constView( f ) )
                ->castTo<const AMP::LinearAlgebra::EpetraVector>()
                .getEpetra_Vector();
            Epetra_Vector &uVec = ( AMP::LinearAlgebra::EpetraVector::view( u ) )
                ->castTo<AMP::LinearAlgebra::EpetraVector>()
                .getEpetra_Vector();
            
            d_mueluSolver->ApplyInverse( fVec, uVec );
        } else {
            AMP_ERROR("Only Epetra interface supported");
        }
    }
    
    // Check for NaNs in the solution (no communication necessary)
    double localNorm = u->localL2Norm();
    AMP_INSIST( localNorm == localNorm, "NaNs detected in solution" );
    
    // we are forced to update the state of u here
    // as Epetra is not going to change the state of a managed vector
    // an example where this will and has caused problems is when the
    // vector is a petsc managed vector being passed back to PETSc
    if ( u->isA<AMP::LinearAlgebra::DataChangeFirer>() ) {
        u->castTo<AMP::LinearAlgebra::DataChangeFirer>().fireDataChange();
    }

    if ( d_iDebugPrintInfoLevel > 2 ) {
        double solution_norm = u->L2Norm();
        AMP::pout << "TrilinosMueLuSolver : after solve solution norm: " << std::setprecision( 15 )
                  << solution_norm << std::endl;
    }

    if ( computeResidual ) {
        d_pOperator->residual( f, u, r );
        finalResNorm = r->L2Norm();

        if ( d_iDebugPrintInfoLevel > 1 ) {
            AMP::pout << "TrilinosMueLuSolver::solve(), L2 norm of residual after solve "
                      << std::setprecision( 15 ) << finalResNorm << std::endl;
        }
    }

    PROFILE_STOP( "solve" );

    if ( d_bRobustMode ) {
        if ( finalResNorm > initialResNorm ) {
            AMP::pout << "Warning: ML was not able to reduce the residual. Using LU instead."
                      << std::endl;
            reSolveWithLU( f, u );
        }
    }
}


void TrilinosMueLuSolver::reSolveWithLU( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                                      AMP::shared_ptr<AMP::LinearAlgebra::Vector>
                                          u )
{
    PROFILE_START( "reSolveWithLU" );

    if ( !d_bUseEpetra ) {
        AMP_ERROR( "Robust mode can only be used with Epetra matrices." );
    }

    AMP::shared_ptr<AMP::Operator::LinearOperator> linearOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>( d_pOperator );
    AMP_INSIST( linearOperator.get() != nullptr, "linearOperator cannot be NULL" );

    d_matrix =
        AMP::dynamic_pointer_cast<AMP::LinearAlgebra::EpetraMatrix>( linearOperator->getMatrix() );
    AMP_INSIST( d_matrix.get() != nullptr, "d_matrix cannot be NULL" );

    auto tmpMueLuParameterList = d_MueLuParameterList;
    tmpMueLuParameterList.set( "verbosity", "medium" );
    tmpMueLuParameterList.set( "max levels", 1 );
    tmpMueLuParameterList.set( "coarse: type", "SuperLU");

    // MueLu expects a Teuchos ref pointer
    Teuchos::RCP<Epetra_CrsMatrix> fineLevelMatrixPtr = Teuchos::rcpFromRef(d_matrix->getEpetra_CrsMatrix());
    auto mueluRCP = MueLu::CreateEpetraPreconditioner(fineLevelMatrixPtr, 
                                                      tmpMueLuParameterList);
        
    d_mueluSolver.reset( mueluRCP.get() );
    mueluRCP.release();

    d_bCreationPhase = false;

    solve( f, u );

    mueluRCP = MueLu::CreateEpetraPreconditioner(fineLevelMatrixPtr, 
                                                 d_MueLuParameterList);
    
    d_mueluSolver.reset( mueluRCP.get() );
    mueluRCP.release();

    d_bCreationPhase = false;

    PROFILE_STOP( "reSolveWithLU" );
}

} // Solver
} // AMP

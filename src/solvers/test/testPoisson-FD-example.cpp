#include "AMP/IO/PIO.h"
#include "AMP/utils/AMPManager.h"

#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/Vector.h"

#include "AMP/discretization/boxMeshDOFManager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/mesh/MeshElement.h"
#include "AMP/mesh/structured/BoxMesh.h"

#include "AMP/matrices/MatrixBuilder.h"

#include "AMP/operators/Operator.h"
#include "AMP/operators/LinearOperator.h"
#include "AMP/operators/OperatorFactory.h"

#include "AMP/solvers/SolverFactory.h"
#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/solvers/SolverFactory.h"

#include <iostream>


// Object to hold column indices and associated data
struct colsDataPair {
    std::vector<size_t> cols;
    std::vector<double> data;
};

/* --------------------------------------
    Implementation of utility functions 
----------------------------------------- */

// Convert an element box to a node box. 
// Modified from src/mesh/test/test_BoxMeshIndex.cpp by removing the possibility of any of the grid dimensions being periodic.
AMP::Mesh::BoxMesh::Box getNodeBox( std::shared_ptr<AMP::Mesh::BoxMesh> mesh, AMP::Mesh::BoxMesh::Box box ) {
    auto global = mesh->getGlobalBox();
    for ( int d = 0; d < 3; d++ ) {
        if ( box.last[d] == global.last[d] )
            box.last[d]++;
    }
    return box;
};

// As above, but just gets a localNodeBox from the mesh
AMP::Mesh::BoxMesh::Box getLocalNodeBox( std::shared_ptr<AMP::Mesh::BoxMesh> mesh ) {
    auto local  = mesh->getLocalBox();
    auto global = mesh->getGlobalBox();
    for ( int d = 0; d < 3; d++ ) {
        if ( local.last[d] == global.last[d] )
            local.last[d]++;
    }
    return local;
};

// As above, but just gets a GlobalNodeBox from the mesh
AMP::Mesh::BoxMesh::Box getGlobalNodeBox( std::shared_ptr<AMP::Mesh::BoxMesh> mesh ) {
    auto global = mesh->getGlobalBox();
    for ( int d = 0; d < 3; d++ ) {
        global.last[d]++;
    }
    return global;
};


/* Fill CSR matrix with data from CSRData */
void fillMatWithLocalCSRData( std::shared_ptr<AMP::LinearAlgebra::Matrix> matrix,
                          std::shared_ptr<AMP::Discretization::DOFManager> DOFMan,
                          std::map<size_t, colsDataPair> CSRData ) {

    // Iterate through local rows in matrix
    for ( size_t dof = DOFMan->beginDOF(); dof != DOFMan->endDOF(); dof++ ) {
        size_t nrows = 1;
        size_t ncols = CSRData[dof].cols.size();

        // It seems that the cols I parse here don't have to be ordered the same as those that are stored internally.
        matrix->setValuesByGlobalID<double>( nrows, ncols, &dof, CSRData[dof].cols.data(), CSRData[dof].data.data() );
    }
}


// Parameters for CG, optionally using a preconditioner
std::unique_ptr<AMP::Database> getCGParameters( bool use_nested ) {
    auto db = std::make_unique<AMP::Database>( "LinearSolver" );
    db->putScalar<std::string>( "name", "CGSolver" );

    db->putScalar<bool>( "uses_preconditioner", use_nested );
    db->putScalar<double>( "absolute_tolerance", 1.0e-12 );
    db->putScalar<double>( "relative_tolerance", 1.0e-12 );
    db->putScalar<int>( "print_info_level", 1 );
    db->putScalar<int>( "max_iterations", 100 );
    // Later, a pc Database will be created with key "Preconditioner" 
    if ( use_nested )
        db->putScalar<std::string>( "pc_name", "Preconditioner" );
    return db;
}


// These parameters are modified from "SolverTestParameters.cpp"
// Boomer docs: https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html
std::unique_ptr<AMP::Database> getBoomerAMGParameters( bool as_preconditioner ) {
    auto db = std::make_unique<AMP::Database>( "LinearSolver" );
    db->putScalar<std::string>( "name", "BoomerAMGSolver" );
    
    // Ensure use of 1 iteration pre preconditioner application
    if (as_preconditioner) {
        db->putScalar<double>( "absolute_tolerance", 0.0 );
        db->putScalar<double>( "relative_tolerance", 0.0 );
        db->putScalar<int>( "max_iterations", 1 );
        db->putScalar<int>( "print_info_level", 0 );
    } else {
        db->putScalar<double>( "absolute_tolerance", 1.0e-12 );
        db->putScalar<double>( "relative_tolerance", 1.0e-12 );
        db->putScalar<int>( "max_iterations", 25 );
        db->putScalar<int>( "print_info_level", 1 );
    }

    db->putScalar<int>( "min_coarse_size", 5 );
    db->putScalar<int>( "relax_type", 6 ); // Hybrid GS/Jacobi
    db->putScalar<int>( "coarsen_type", 10 );
    db->putScalar<int>( "interp_type", 17 ); // distance-two
    db->putScalar<int>( "cycle_type", 1 );
    db->putScalar<int>( "relax_order", 0 );
    db->putScalar<double>( "strong_threshold", 0.5 );
    return db;
}


/* Create a SolverStrategyParameters object, as required to build a linear solver from the SolverFactory. input_db should contain a Database with key "solver_name", which contains a "name" key that is a valid class in AMP::Solver. 
    For similar examples, see: 
        "buildSolver()" in testSolversForUserMatrix.cpp 
        "buildSolver()" in SolverTestParameters.cpp */
std::shared_ptr<AMP::Solver::SolverStrategyParameters> 
getSolverStrategyParameters(const std::string &solver_name,
                            std::shared_ptr<AMP::Database> input_db,
                            const AMP::AMP_MPI &comm,
                            std::shared_ptr<AMP::LinearAlgebra::Vector> initialGuess,
                            std::shared_ptr<AMP::Operator::Operator> op )
{
    // Ensure input_db contains a database with key "solver_name" 
    AMP_INSIST( input_db->keyExists( solver_name ), "Key " + solver_name + " is missing!" );
    auto db = input_db->getDatabase( solver_name );

    // Ensure the "name" of the solver is specified
    AMP_INSIST( db->keyExists( "name" ), "input_db must contain a 'name' key inside its '" + solver_name + "' database" );
    auto name = db->getScalar<std::string>( "name" );

    auto solverStratParams = std::make_shared<AMP::Solver::SolverStrategyParameters>( db );

    solverStratParams->d_pOperator     = op;
    solverStratParams->d_comm          = comm;
    solverStratParams->d_pInitialGuess = initialGuess;
    solverStratParams->d_db            = db;
    solverStratParams->d_global_db     = input_db;

    return solverStratParams;
}


/* Given a solver defined by solverStratParams, augment it with a preconditioner if one is specified to be used in input_db[solver_name] */
void augmentWithPcSolver(const std::string &solver_name,
                        std::shared_ptr<AMP::Database> input_db,
                        std::shared_ptr<AMP::Solver::SolverStrategyParameters> solverStratParams ) {

    // First some basic error checking
    // Ensure input_db contains a database with key "solver_name" 
    AMP_INSIST( input_db->keyExists( solver_name ), "Key " + solver_name + " is missing!" );
    auto db = input_db->getDatabase( solver_name );
    // Ensure the "name" of the solver is specified
    AMP_INSIST( db->keyExists( "name" ), "input_db must contain a 'name' key inside its '" + solver_name + "' database" );
    auto name = db->getScalar<std::string>( "name" );

    // Construct preconditioner if one is to be used 
    auto uses_preconditioner = db->getWithDefault<bool>( "uses_preconditioner", false );
    if ( uses_preconditioner ) {
        auto pc_name = db->getWithDefault<std::string>( "pc_name", "Preconditioner" );
        AMP_INSIST( input_db->keyExists( pc_name ), "'" + solver_name + "' uses a pc with name '" + pc_name + "' but input_db does not contain a key for it" );

        // Get a solverStategyParameters object for the preconditioner
        auto pcSolverParams = getSolverStrategyParameters(pc_name, input_db, solverStratParams->d_comm, nullptr, solverStratParams->d_pOperator );

        // Set NestedSolver field for the outer solver as the preconditioner
        solverStratParams->d_pNestedSolver = AMP::Solver::SolverFactory::create( pcSolverParams );
    }
}


/* Create a linear solver for the linear operator AOp. 
    solverName can be one of: "CG", "AMG", "CG+AMG" */
static std::unique_ptr<AMP::Solver::SolverStrategy>
getLinearSolver(AMP::AMP_MPI comm,
                std::shared_ptr<AMP::Operator::LinearOperator> AOp,
                const AMP::Database &PDE_db) {

    auto solverName = PDE_db.getString( "solverName" );

    // Get solver-specific parameters. Parameters for the solver are stored as a Database in solverParams with key "LinearSolver"
    auto solverParams = std::make_shared<AMP::Database>();
    if ( solverName == "CG" ) {
        // CG params when used as a solver
        solverParams->putDatabase( "LinearSolver", getCGParameters( false ) ); 
    } else if ( solverName == "AMG" ) {
        // BoomerAMG params when used as a solver
        solverParams->putDatabase( "LinearSolver", getBoomerAMGParameters( false ) ); 
    } else if ( solverName == "CG+AMG" ) {
        // CG parameters, when used with a preconditioner
        solverParams->putDatabase( "LinearSolver", getCGParameters( true ) ); 
        // BoomerAMG params when used as a preconditioner stored as "Preconditioner" Database
        solverParams->putDatabase( "Preconditioner", getBoomerAMGParameters( true ) ); 
    } else {
        AMP_ERROR( "solverName '" + solverName + "' is not one of those supported: 'CG', 'AMG', or 'CG+AMG'" );
    }

    // Create solverStrategyParameters according to the "LinearSolver" Database in solverParams
    auto solverStratParams = getSolverStrategyParameters( "LinearSolver", solverParams, comm, nullptr, AOp );
    // Add preconditioner to outer solver if required
    augmentWithPcSolver( "LinearSolver", solverParams, solverStratParams );
    // Create linear solver
    auto solver = AMP::Solver::SolverFactory::create( solverStratParams );

    return solver;
}


/* Create a d-dimensional BoxMesh over [0,1]^d with n+1 points in each direction */
static std::shared_ptr<AMP::Mesh::BoxMesh> createBoxMesh( AMP::AMP_MPI comm, const AMP::Database &PDE_db )
{
    auto n   = PDE_db.getScalar<int>( "n" );
    auto dim = PDE_db.getScalar<int>( "dim" );

    auto mesh_db = std::make_shared<AMP::Database>( "Mesh" );
    mesh_db->putScalar<int>( "dim", dim );
    mesh_db->putScalar<std::string>( "MeshName", "AMP::cube" );
    mesh_db->putScalar<std::string>( "Generator", "cube" );
    if ( dim == 1 ) {
        mesh_db->putVector<int>( "Size", { n } ); // mesh has n+1 points
        mesh_db->putVector<double>( "Range", { 0.0, 1.0 } );
    } else if ( dim == 2 ) {
        mesh_db->putVector<int>( "Size", { n, n } ); // mesh has n+1 x n+1 points
        mesh_db->putVector<double>( "Range", { 0.0, 1.0, 0.0, 1.0 } );
    } else if ( dim == 3 ) {
        mesh_db->putVector<int>( "Size", { n, n, n } ); // mesh has n+1 x n+1 x n+1 points
        mesh_db->putVector<double>( "Range", { 0.0, 1.0, 0.0, 1.0, 0.0, 1.0 } );
    }

    // Create MeshParameters
    auto mesh_params = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    mesh_params->setComm( comm );
    
    // Create Mesh
    static std::shared_ptr<AMP::Mesh::BoxMesh> boxMesh = AMP::Mesh::BoxMesh::generate( mesh_params );

    // Print some statistics of mesh
    # if 0
    std::cout << "Mesh: " 
    << static_cast<int>(boxMesh->getDim()) << "-dimensional. " 
    << "Global size = " << boxMesh->numGlobalElements(AMP::Mesh::GeomType::Vertex) << ". " 
    << "Local size on rank " << boxMesh->getComm().getRank() << " (of " << boxMesh->getComm().getSize() << ") = " << boxMesh->numLocalElements(AMP::Mesh::GeomType::Vertex) 
    << std::endl;
    //std::cout << "mesh comm size = " << boxMesh->getComm().getSize() << std::endl;
    #endif
    
    #if 0
    // Bounding box of mesh w/ zero ghost
    auto globalBox = getNodeBox( boxMesh, boxMesh->getGlobalBox(0) );
    auto localBox  = getNodeBox( boxMesh, boxMesh->getLocalBox(0) );
    std::cout << " local box on p" << comm.getRank() << ": ";
    for ( int i = localBox.first[0]; i <= localBox.last[0]; i++ ) {
        std::cout << i << ",";
    }
    std::cout << std::endl;

    std::cout << "global box on p" << comm.getRank() << ": ";
    for ( int i = globalBox.first[0]; i <= globalBox.last[0]; i++ ) {
        std::cout << i << ",";
    }
    std::cout << std::endl;
    #endif

    #if 0
    AMP::pout << "boundary ids:" << std::endl;
    auto bndry_ids = boxMesh->getBoundaryIDs();
    for (auto id : bndry_ids) {
        AMP::pout << id << " ";
    }
    AMP::pout << std::endl;
    #endif

    return boxMesh;
}



/* Compute discrete norms of vector u */
std::vector<double> getDiscreteNorms(double h,  
                    std::shared_ptr<const AMP::LinearAlgebra::Vector> u) {
    // Compute norms
    double uL1Norm  = static_cast<double>( u->L1Norm()  ) * h*h;
    double uL2Norm  = static_cast<double>( u->L2Norm()  ) * h;
    double uMaxNorm = static_cast<double>( u->maxNorm() );

    std::vector<double> unorms = { uL1Norm, uL2Norm, uMaxNorm }; 
    return unorms;
}



/* ------------------------------------------------
    Class implementing a discrete Poisson problem 
------------------------------------------------- */
class PoissonOp : public AMP::Operator::LinearOperator {

private:
    std::shared_ptr<AMP::LinearAlgebra::Matrix> getLaplacianMatrix();

public:

    std::shared_ptr<AMP::Mesh::BoxMesh>              d_BoxMesh;
    std::shared_ptr<AMP::Database>                   d_db;
    std::shared_ptr<AMP::Discretization::DOFManager> d_DOFMan;
    AMP::Mesh::GeomType                              d_geomType = AMP::Mesh::GeomType::Vertex;

    // Constructor call's base class's constructor
    PoissonOp(std::shared_ptr<const AMP::Operator::OperatorParameters> params_) : 
            AMP::Operator::LinearOperator( params_ ) { 

        // Keep a pointer to my BoxMesh to save having to do this downcast repeatedly 
        d_BoxMesh = std::dynamic_pointer_cast<AMP::Mesh::BoxMesh>(this->getMesh());
        AMP_INSIST( d_BoxMesh, "Mesh must be a AMP::Mesh::BoxMesh" );

        // Set my database
        d_db = params_->d_db->getDatabase( "PDE_db" );

        // Set DOFManager
        this->set_DOFManager();
        AMP_INSIST(  d_DOFMan, "Requires non-null DOF" );

        // Get the matrix
        auto A = getLaplacianMatrix( );
        
        // Set linear operator's matrix
        this->setMatrix( A );

    }

    // Used by OperatorFactory to create a PoissonOp
    static std::unique_ptr<AMP::Operator::Operator> create( std::shared_ptr<AMP::Operator::OperatorParameters> params ) {  
        return std::make_unique<PoissonOp>( params ); };

    std::shared_ptr<AMP::LinearAlgebra::Vector> getRightVector() const override {
        auto tempVar = std::make_shared<AMP::LinearAlgebra::Variable> (" ");
        return AMP::LinearAlgebra::createVector( this->d_DOFMan, tempVar );
    }

    // Used to register this operator in a factory
    std::string type() const override {
        return "PoissonOp";
    }

    // Populate vectors with exact PDE solution and corresponing source term
    void fill_uexact_and_fsource( 
        std::shared_ptr<AMP::LinearAlgebra::Vector> uexact, 
        std::shared_ptr<AMP::LinearAlgebra::Vector> fsource );
    

private:
    /* Build and set d_DOFMan */
    void set_DOFManager() {
        int DOFsPerElement = 1; 
        int gcw  = 1; // Ghost-cell width. 
        d_DOFMan = AMP::Discretization::boxMeshDOFManager::create(this->getMesh(), d_geomType, gcw, DOFsPerElement);
    }

    // Exact solution and corresponding source term
    double exactSolution(double x);
    double exactSolution(double z, double y);
    double exactSolution(double x, double y, double z);
    double sourceTerm(double x);
    double sourceTerm(double x, double y);
    double sourceTerm(double x, double y, double z);
    
    // FD stencils
    std::vector<double> getStencil1D();
    std::vector<double> getStencil2D();
    std::vector<double> getStencil3D();

    // Map from row index to col+data in that row
    std::map<size_t, colsDataPair> get1DCSRData();
    std::map<size_t, colsDataPair> get2DCSRData();
    std::map<size_t, colsDataPair> get3DCSRData();

    // Map from grid index i (or i,j, or i,j,k) to a MeshElementIndex to a MeshElementId and then to the corresponding DOF
    size_t grid_inds_to_DOF( int i, int j = 0, int k = 0 ) {
        AMP::Mesh::BoxMesh::MeshElementIndex ind(
                        AMP::Mesh::GeomType::Vertex, 0, i, j, k );
        AMP::Mesh::MeshElementID id = d_BoxMesh->convert( ind );
        std::vector<size_t> dof;
        d_DOFMan->getDOFs(id, dof);
        return dof[0];
    };
}; 


// -u_xx = f
double PoissonOp::exactSolution(double x) {
    return sin(2.0 * M_PI * x);
}

// -u_xx = f
double PoissonOp::sourceTerm(double x) {
    return 4.0 * M_PI * M_PI * sin(2.0 * M_PI * x);
}

// -(alpha*u_xx + beta*u_yy + gamma*u_xy) = f (make sol non-symmetric in x and y to ensure the code doesn't mix up the two)
double PoissonOp::exactSolution(double x, double y) {
    double u = sin(2.0 * M_PI * x) * sin(4.0 * M_PI * y);
    return u;
}

// -(alpha*u_xx + beta*u_yy + gamma*u_xy) = f
double PoissonOp::sourceTerm(double x, double y) {
    
    // Unpack parameters
    auto theta   = d_db->getScalar<double>("theta");
    auto eps     = d_db->getScalar<double>("eps");
    double c     = cos(theta);
    double s     = sin(theta);
    // PDE coefficients
    double alpha = c*c + eps * s*s;
    double beta  = eps * c*c + s*s;
    double gamma = 2*(1 - eps) * c * s;

    double f = 4*M_PI*M_PI*(alpha*sin(2*M_PI*x)*sin(4*M_PI*y) + 4*beta*sin(2*M_PI*x)*sin(4*M_PI*y) - 2*gamma*cos(2*M_PI*x)*cos(4*M_PI*y));

    return f;
}

// -(u_xx + u_yy + epsilon*u_zz) = f
double PoissonOp::exactSolution(double x, double y, double z) {
    return std::sin(2*M_PI*x)*std::sin(4*M_PI*y)*std::sin(6*M_PI*z);
}

// -(u_xx + u_yy + epsilon*u_zz) = f
double PoissonOp::sourceTerm(double x, double y, double z) {
    auto epsilon = d_db->getScalar<double>("eps");
    return 4*std::pow(M_PI, 2)*(9*epsilon + 5)*std::sin(2*M_PI*x)*std::sin(4*M_PI*y)*std::sin(6*M_PI*z);
}


/* Get 3-point stencil for 1D Poisson. */
std::vector<double> PoissonOp::getStencil1D() {
    
    // Stencil coefficients
    double W  = -1.0; 
    double O  = +2.0; 
    double E  = -1.0;
    // Populate stencil
    std::vector<double> stencil = { O, W, E };

    // Introduce 1/h^2 scaling 
    auto h = d_db->getScalar<double>("h");
    for ( auto &s : stencil ) {
        s *= 1.0/(h*h);
    }

    return stencil;
}

/* Get 9-point stencil for upwind FD disretization of rotated anisotropic 2D Poisson. */
std::vector<double> PoissonOp::getStencil2D() {
    
    double eps   = this->d_db->getScalar<double>( "eps" );
    double theta = this->d_db->getScalar<double>( "theta" );

    AMP_INSIST( theta >= 0.0 && theta <= M_PI/2.0, "Upwind discretization only valid for theta in [0,pi/2]" );

    double c     = cos(theta);
    double s     = sin(theta);

    // PDE coefficients
    double alpha = c*c + eps * s*s;
    double beta  = eps * c*c + s*s;
    double gamma = 2*(1 - eps) * c * s;
    gamma *= 0.5; // gamma only ever appears multiplied by 1/2, except at O.
    double NW = 0.0;
    double N  =           -   beta +   gamma; 
    double NE =                    -   gamma; 
    double W  = -   alpha          +   gamma; 
    double O  = + 2*alpha + 2*beta - 2*gamma; 
    double E  = -   alpha          +   gamma;
    double SW =                    -   gamma;
    double S  =           -   beta +   gamma;
    double SE = 0.0; 

    // Populate stencil
    std::vector<double> stencil = { O, SW, S, SE, W, E, NW, N, NE };

    // Introduce 1/h^2 scaling 
    auto h = d_db->getScalar<double>("h");
    for ( auto &s : stencil ) {
        s *= 1.0/(h*h);
    }

    return stencil;
}

/* Get 7-point stencil for 3D Poisson. */
std::vector<double> PoissonOp::getStencil3D() {
    
    auto eps  = d_db->getScalar<double>("eps");

    double D  = -1.0*eps;
    double S  = -1.0;
    double W  = -1.0; 
    double O  = +4.0 + 2.0*eps; 
    double E  = -1.0;
    double N  = -1.0;
    double U  = -1.0*eps;

    // Populate stencil
    std::vector<double> stencil = { O, D, S, W, E, N, U };

    // Introduce 1/h^2 scaling 
    auto h = d_db->getScalar<double>("h");
    for ( auto &s : stencil ) {
        s *= 1.0/(h*h);
    }

    return stencil;
}


// CSR structure for identity row
colsDataPair localCSR_identity(size_t dof) { 
    std::vector<double> vals = { 1.0 };  // data
    std::vector<size_t> cols = { dof };  // Column index
    colsDataPair identity = { cols, vals }; 
    return identity;
};


/* Get CSR structure of 1D Laplacian */
std::map<size_t, colsDataPair> PoissonOp::get1DCSRData() {    

    // Get 3-point stencil
    auto stencil = getStencil1D( );

    // Get local grid index box w/ zero ghosts
    auto globalBox = getGlobalNodeBox( d_BoxMesh );

    // Segment localBox into DOFs on the global boundary and those on the interior
    auto localBoxInterior = getLocalNodeBox( d_BoxMesh );
    std::vector<int> localBoundaryIDs = {-1, -1};
    // Local west boundary is global west boundary
    if (localBoxInterior.first[0] == globalBox.first[0]) {
        localBoundaryIDs[0] = localBoxInterior.first[0];
        localBoxInterior.first[0] += 1;
    }
    // Local east boundary is global east boundary
    if (localBoxInterior.last[0] == globalBox.last[0]) {
        localBoundaryIDs[1] = localBoxInterior.last[0];
        localBoxInterior.last[0] -= 1;
    }

    // Create a map from the DOF to a pair a vectors
    // Map from a DOF to vector of col inds and associated data
    std::map<size_t, colsDataPair> localCSRData;

    // Set identity in boundary rows that I own
    for (auto i : localBoundaryIDs) {
        if (i != -1) {
            auto dof = grid_inds_to_DOF( i );
            localCSRData[dof] = localCSR_identity( dof );
        }
    }

    // Iterate over local interior box
    for (auto i = localBoxInterior.first[0]; i <= localBoxInterior.last[0]; i++) {
        
        // The current row
        size_t dof = grid_inds_to_DOF( i );
        // Copy of stencil
        std::vector<double> vals = stencil; 
        // Column indices, ordered consistently with the stencil
        std::vector<size_t> cols = { 
            dof,
            grid_inds_to_DOF( i-1 ),
            grid_inds_to_DOF( i+1 ) };
        localCSRData[dof] = { cols, vals };
    }  

    return localCSRData;
}

/* Get CSR structure of rotated anisotropic 2D Laplacian */
std::map<size_t, colsDataPair> PoissonOp::get2DCSRData() {    

    // Get 9-point stencil 
    auto stencil = getStencil2D( );
    
    // Get local grid index box w/ zero ghosts
    auto localBox  = getLocalNodeBox( d_BoxMesh );
    auto globalBox = getGlobalNodeBox( d_BoxMesh );

    // Create a map from the DOF to a pair a vectors
    // Map from a DOF to vector of col inds and associated data
    std::map<size_t, colsDataPair> localCSRData;

    // Iterate over local box
    for (auto j = localBox.first[1]; j <= localBox.last[1]; j++) {
        for (auto i = localBox.first[0]; i <= localBox.last[0]; i++) {
        
            // The current row
            size_t dof = grid_inds_to_DOF( i, j );

            // Set identity in boundary rows
            if (j == globalBox.first[1] || j == globalBox.last[1] || i == globalBox.first[0] || i == globalBox.last[0]) {
                localCSRData[dof] = localCSR_identity( dof );
                continue;
            }
            
            // Copy of stencil
            std::vector<double> vals = stencil; 
            // Column indices, ordered consistently with the stencil
            std::vector<size_t> cols = { 
                dof,
                grid_inds_to_DOF( i-1, j-1 ),
                grid_inds_to_DOF( i ,  j-1 ),
                grid_inds_to_DOF( i+1, j-1 ),
                grid_inds_to_DOF( i-1, j   ),
                grid_inds_to_DOF( i+1, j   ),
                grid_inds_to_DOF( i-1, j+1 ),
                grid_inds_to_DOF( i ,  j+1 ),
                grid_inds_to_DOF( i+1, j+1 ) };

            localCSRData[dof] = { cols, vals };
        }
    }   

    return localCSRData;
}

/* Get CSR structure of anisotropic 3D Laplacian */
std::map<size_t, colsDataPair> PoissonOp::get3DCSRData( ) {    

    // Get 7-point stencil
    auto stencil = getStencil3D( );

    // Get local grid index box w/ zero ghosts
    auto localBox  = getLocalNodeBox( d_BoxMesh );
    auto globalBox = getGlobalNodeBox( d_BoxMesh );

    // Create a map from the DOF to a pair a vectors
    // Map from a DOF to vector of col inds and associated data
    std::map<size_t, colsDataPair> localCSRData;

    // Iterate over local box
    for (auto k = localBox.first[2]; k <= localBox.last[2]; k++) {
        for (auto j = localBox.first[1]; j <= localBox.last[1]; j++) {
            for (auto i = localBox.first[0]; i <= localBox.last[0]; i++) {

                // The current row
                size_t dof = grid_inds_to_DOF( i, j, k );

                // Set identity in boundary rows
                if (k == globalBox.first[2] || k == globalBox.last[2] || j == globalBox.first[1] || j == globalBox.last[1] || i == globalBox.first[0] || i == globalBox.last[0]) {
                    localCSRData[dof] = localCSR_identity( dof );
                    continue;
                }

                // Copy of stencil
                std::vector<double> vals = stencil; 
                // Column indices, ordered consistently with the stencil
                // O, D, S, W, E, N, U
                std::vector<size_t> cols = { 
                    dof,
                    grid_inds_to_DOF( i,   j,   k-1 ),
                    grid_inds_to_DOF( i,   j-1, k   ),
                    grid_inds_to_DOF( i-1, j,   k   ),
                    grid_inds_to_DOF( i+1, j,   k   ),
                    grid_inds_to_DOF( i,   j+1, k   ),
                    grid_inds_to_DOF( i,   j,   k+1 ) };

                localCSRData[dof] = { cols, vals };
            }
        } 
    }

    return localCSRData;
}



/* Populate exact solution and RHS source-term vectors. 
    Note that zero values are put into boundary rows of the 
    RHS vector f since the problem is posed with zero Dirichlet BCs */
void PoissonOp::fill_uexact_and_fsource( 
        std::shared_ptr<AMP::LinearAlgebra::Vector> uexact, 
        std::shared_ptr<AMP::LinearAlgebra::Vector> fsource ) {

    auto it      = d_BoxMesh->getIterator(d_geomType); // Mesh iterator
    auto meshDim = d_BoxMesh->getDim(); // Dimension

    // Fill in exact solution and source term vectors
    for ( auto elem = it.begin(); elem != it.end(); elem++ ) {
        
        std::vector<double> u;
        std::vector<double> f;

        if ( meshDim == 1 ) {
            double x = ( elem->coord() )[0];
            u.push_back(exactSolution(x));
            
            if (!(elem->isOnSurface())) {
                f.push_back(sourceTerm( x ));
            } else {
                f.push_back(0.0);
            }

        } else if ( meshDim == 2 ) {
            double x = ( elem->coord() )[0];
            double y = ( elem->coord() )[1];
            u.push_back( exactSolution( x, y ) );

            if (!(elem->isOnSurface())) {
                f.push_back( sourceTerm( x, y ) );
            } else {
                f.push_back( 0.0 );
            }

        } else if (meshDim == 3) {
            double x = ( elem->coord() )[0];
            double y = ( elem->coord() )[1];
            double z = ( elem->coord() )[2];
            u.push_back( exactSolution( x, y, z ) );

            if (!(elem->isOnSurface())) {
                f.push_back( sourceTerm( x, y, z ) );
            } else {
                f.push_back( 0.0 );
            }
        }

        std::vector<size_t> i;
        d_DOFMan->getDOFs( elem->globalID(), i );
        uexact->setValuesByGlobalID( 1, &i[0], &u[0] );
        fsource->setValuesByGlobalID( 1, &i[0], &f[0] );
    }
    uexact->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
    fsource->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}


/* Return a constructed CSR matrix corresponding to the discretized Laplacian on the mesh */
std::shared_ptr<AMP::LinearAlgebra::Matrix> PoissonOp::getLaplacianMatrix( ) {

    auto inVec  = this->getRightVector();
    auto outVec = this->getRightVector();

    int meshDim = this->getMesh()->getDim();
    std::map<size_t, colsDataPair> localCSRData;
    if ( meshDim == 1 ) {
        localCSRData = this->get1DCSRData( );
    } else if ( meshDim == 2 ) {
        localCSRData = this->get2DCSRData( );
    } else if ( meshDim == 3 ) {
        localCSRData = this->get3DCSRData( );
    }

    // Create Lambda to return col inds from a given row ind
    auto getColumnIDs = [&](int row) { return localCSRData[row].cols; };
    // Create CSR matrix
    std::shared_ptr<AMP::LinearAlgebra::Matrix> A = AMP::LinearAlgebra::createMatrix( inVec, outVec, "CSRMatrix", getColumnIDs );
    fillMatWithLocalCSRData( A, d_DOFMan, localCSRData );
    
    // Finalize A
    A->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );
    // Print statistics of A
    #if 0
    size_t nGlobalRows = A->numGlobalRows();
    size_t nLocalRows  = A->numLocalRows();
    std::cout << "Matrix A: #global rows=" << nGlobalRows << ". #local rows on rank " << A->getComm().getRank() << " (of " << A->getComm().getSize() << ") = " << nLocalRows
              << std::endl;
    #endif
    return A;
}



void driver(AMP::AMP_MPI comm, 
            std::shared_ptr<AMP::Database> PDE_db ) {

    AMP_INSIST( PDE_db, "driver() requires non-null PDE_db" );

    /****************************************************************
    * Create a mesh                                                 *
    ****************************************************************/
    static std::shared_ptr<AMP::Mesh::BoxMesh> mesh = createBoxMesh( comm, *PDE_db );


    /****************************************************************
    * Create the LinearOperator PoissonOp over the mesh             *
    ****************************************************************/

    const auto OpDB = std::make_shared<AMP::Database>( "linearOperatorDB" );
    OpDB->putDatabase( "PDE_db", PDE_db->cloneDatabase() );
    OpDB->putScalar<int>( "print_info_level", 0 );
    OpDB->putScalar<std::string>( "name", "PoissonOp" ); // Operator factory requires the db contain a "name" key providing the name of the operator   
    auto OpParameters = std::make_shared<AMP::Operator::OperatorParameters>( OpDB );
    
    OpParameters->d_name = "PoissonOp";
    OpParameters->d_Mesh = mesh;

    // There are two different methods for creating a PoissonOp
    // This is method 1, using directly the PoissonOp constructor
    //auto AOp             = std::make_shared<PoissonOp>( OpParameters );    
    
    // This is method 2, using an operator factory. In this particular example, method 2 is more complicated, but in general having the factory infrastructure set up allows for more flexibility because then solvers, etc. can create a PoissonOp on their own.
    // Create an OperatorFactory and register PoissonOp in it 
    auto & operatorFactory = AMP::Operator::OperatorFactory::getFactory();
    operatorFactory.registerFactory( "PoissonOp", PoissonOp::create );
    std::shared_ptr<AMP::Operator::Operator> AOp_ = AMP::Operator::OperatorFactory::create( OpParameters );
    auto AOp = std::dynamic_pointer_cast<PoissonOp>( AOp_ );


    /****************************************************************
    * Set up relevant vectors over the mesh                         *
    ****************************************************************/
    // Create required vectors over the mesh
    auto unumVec    = AOp->getRightVector();
    auto uexactVec  = AOp->getRightVector();
    auto rexactVec  = AOp->getRightVector();
    auto fsourceVec = AOp->getRightVector();

    // Set exact solution and RHS vector
    AOp->fill_uexact_and_fsource( uexactVec, fsourceVec );

    // Initialize unum to random values
    unumVec->setRandomValues();
    unumVec->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
    

    /****************************************************************
    * Compute discrete residual norm on continuous solution (truncation error) *
    ****************************************************************/
    auto A = AOp->getMatrix();
    AMP::pout << "\nDiscrete residual of continuous solution: ";
    AOp->residual( fsourceVec, uexactVec, rexactVec );
    auto rnorms = getDiscreteNorms( PDE_db->getScalar<double>( "h" ), rexactVec );
    // Print residual norms
    AMP::pout << "||r|| = (" << rnorms[0] << ", " << rnorms[1] << ", " << rnorms[2] << ")" << std::endl;

    
    /****************************************************************
    * Construct linear solver of the LinearOperator and apply it    *
    ****************************************************************/
    // Get the linear solver for operator AOp
    auto mySolver = getLinearSolver( comm, AOp, *PDE_db );

    // Use zero initial iterate and apply solver
    mySolver->setZeroInitialGuess( true );
    mySolver->apply(fsourceVec, unumVec);
    
    // Compute disretization error
    AMP::pout << "\nDiscretization error post linear solve: "; 
    auto e = uexactVec->clone();
    e->axpy(-1.0, *unumVec, *uexactVec); 
    auto enorms = getDiscreteNorms( PDE_db->getScalar<double>( "h" ), e );
    AMP::pout << "||e|| = (" << enorms[0] << ", " << enorms[1] << ", " << enorms[2] << ")" << std::endl;

}
// end of driver()




/*  Input usage is: >> <dim> <n> <solverName>
    Where:
        dim        : 1, 2, or 3 is the number of spatial dimensions
        n          : There are n+1 mesh points in each dimension
        solverName : CG, or AMG, or CG+AMG is the linear solver

    e.g., >> mpirun -n 2 testPoisson-FD-example 1 16 CG+AMG

    In 1D, the PDE is -u_xx = f. Standard 3-point finite differences are used. 
    In 2D, the PDE is -u_xx -eps*u_yy = f but rotated an angle of theta degrees. 7-point upwind finite-differences are used.
    In 3D, the PDE is -u_xx -u_yy - eps*u_zz = f

*/
int main( int argc, char **argv )
{

    // Default constants in the 2D PDE.
    double eps   = 1e-2;
    double theta = 30.0 * (M_PI / 180.0);

    AMP::AMPManager::startup( argc, argv );

    // Create a global communicator
    AMP::AMP_MPI comm( AMP_COMM_WORLD );
    // int myRank   = comm.getRank();
    // int numRanks = comm.getSize();

    // Unpack inputs
    int dim = atoi(argv[1]);          // PDE dimension
    int n   = atoi(argv[2]);          // Grid size in each dimension
    std::string solverName = argv[3]; // Solver name
    
    // Create DB with PDE- and mesh-related quantities
    auto PDE_db = std::make_shared<AMP::Database>( "PDE" );
    PDE_db->putScalar( "n",   n );
    PDE_db->putScalar( "dim", dim );
    PDE_db->putScalar( "solverName", solverName );
    if ( dim == 2 ) {
        PDE_db->putScalar<double>( "eps",   eps );
        PDE_db->putScalar<double>( "theta", theta );
    } else if ( dim == 3 ) {
        PDE_db->putScalar<double>( "eps",   eps );
    }

    // Store the mesh size: There are n+1 points on the mesh, including the boundaries
    double h = 1.0 / n;
    PDE_db->putScalar( "h",   h );

    // Driver
    driver( comm, PDE_db );

    AMP::AMPManager::shutdown();

    return 0;
}
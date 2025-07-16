#include "AMP/IO/PIO.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/IO/AsciiWriter.h"

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

#include "AMP/solvers/testHelpers/testSolverHelpers.h"

#include "AMP/time_integrators/TimeIntegratorParameters.h"
#include "AMP/time_integrators/TimeIntegratorFactory.h"
#include "AMP/time_integrators/TimeIntegrator.h"
#include "AMP/time_integrators/TimeOperator.h"
#include "AMP/time_integrators/ImplicitIntegrator.h"

#include <iostream>
#include <iomanip>
#include <filesystem>


/* Representation of constant coefficients in front of derivatives in second-order PDE
In 1D:
    cxx*u_xx
In 2D:
    cxx*u_xx + cyy*u_yy + cxy*u_xy
In 3D:
    cxx*u_xx + cyy*u_yy + czz*u_zz + cxy*u_xy + cxz*u_xz + cyz*u_yz
*/
struct PDECoefficients {
    double xx = std::nan("");
    double yy = std::nan("");
    double zz = std::nan("");
    double xy = std::nan("");
    double xz = std::nan("");
    double yz = std::nan("");

    PDECoefficients() {};
    PDECoefficients( int dim, std::vector<double> c ) {
        if ( dim == 1 ) {
            xx = c[0];
        } else if ( dim == 2 ) {
            xx = c[0], yy = c[1], xy = c[2];
        } else if ( dim == 3 ) {
            xx = c[0], yy = c[1], zz = c[2], xy = c[3], xz = c[4], yz = c[5];
        }
    }
};

// Object to hold column indices and associated data for packing into a CSR matrix
struct colsDataPair {
    std::vector<size_t> cols;
    std::vector<double> data;
};



/* --------------------------------------
    Implementation of utility functions 
----------------------------------------- */

// CSR structure for identity row
colsDataPair localCSR_identity(size_t dof) { 
    std::vector<double> vals = { 1.0 };  // data
    std::vector<size_t> cols = { dof };  // Column index
    colsDataPair identity    = { cols, vals }; 
    return identity;
}

// Convert an element box to a node box. 
// Modified from src/mesh/test/test_BoxMeshIndex.cpp by removing the possibility of any of the grid dimensions being periodic.
AMP::Mesh::BoxMesh::Box getNodeBox( std::shared_ptr<AMP::Mesh::BoxMesh> mesh, AMP::Mesh::BoxMesh::Box box ) {
    auto global = mesh->getGlobalBox();
    for ( int d = 0; d < 3; d++ ) {
        if ( box.last[d] == global.last[d] )
            box.last[d]++;
    }
    return box;
}

// As above, but just gets a localNodeBox from the mesh
AMP::Mesh::BoxMesh::Box getLocalNodeBox( std::shared_ptr<AMP::Mesh::BoxMesh> mesh ) {
    auto local  = mesh->getLocalBox();
    auto global = mesh->getGlobalBox();
    for ( int d = 0; d < 3; d++ ) {
        if ( local.last[d] == global.last[d] )
            local.last[d]++;
    }
    return local;
}

// As above, but just gets a GlobalNodeBox from the mesh
AMP::Mesh::BoxMesh::Box getGlobalNodeBox( std::shared_ptr<AMP::Mesh::BoxMesh> mesh ) {
    auto global = mesh->getGlobalBox();
    for ( int d = 0; d < 3; d++ ) {
        global.last[d]++;
    }
    return global;
}


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


/* In 1D: We want the computational vertices at the mid points of each cell. 
The n cells look like:
    [0, 1]*h, [1, 2]*h, ..., [n-1, n]*h.
Since there are n cells in [0,1], the width of each must be 1/n
The computational vertices are the mid points of each cell, but note that we eliminate the first and last vertices 0.5*h and 1 - 0.5*h, so we actually have n-2 active vertices:
        h + 0.5*h, ..., 1 - (h + 0.5*). 

    */
std::shared_ptr<AMP::Mesh::BoxMesh> createBoxMesh( AMP::AMP_MPI comm, std::shared_ptr<AMP::Database> PDE_db )
{
    auto n   = PDE_db->getScalar<int>( "n" );
    auto dim = PDE_db->getScalar<int>( "dim" );
    double h = 1.0 / n;
    PDE_db->putScalar( "h",   h );
    
    auto mesh_db = std::make_shared<AMP::Database>( "Mesh" );
    mesh_db->putScalar<int>( "dim", dim );
    mesh_db->putScalar<std::string>( "MeshName", "AMP::cube" );
    mesh_db->putScalar<std::string>( "Generator", "cube" );
    if ( dim == 1 ) {
        //mesh_db->putVector<int>( "Size", { n-3 } ); // mesh has n-2 points
        //mesh_db->putVector<double>( "Range", { h + 0.5*h, 1.0 - (h + 0.5*h) } );

        mesh_db->putVector<int>( "Size", { n-1 } ); // mesh has n points
        mesh_db->putVector<double>( "Range", { 0.5*h, 1.0 - 0.5*h } );
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

    // Check vertices are where I think they are
    #if 1
        auto iter = boxMesh->getIterator( AMP::Mesh::GeomType::Vertex );
        std::cout << "mesh spacing is h=" << h << std::endl;
        std::cout << "mesh vertices are:" << std::endl;
        for (auto node : iter) {
            std::cout << node.coord(0) << ", ";
        }
        std::cout << std::endl; 
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
    Class implementing manufactured solution 
------------------------------------------------- */
/* Abstract class representing an exact solution of a linear diffusion problem:
        u'(t) - grad \dot ( \grad u ) = s(t), u(0) = u_0
    over the spatial domain [0,1]^d.

    Robin boudary conditions are avaliable in the form of:
        a0 * u + b0 * n \dot grad u = r0 at x = 0...
        a1 * u + b1 * n \dot grad u = r1 at x = 1...
*/
class ManufacturedHeatModel {

public:

    AMP::Mesh::GeomType                              d_geomType = AMP::Mesh::GeomType::Vertex;
    std::shared_ptr<AMP::Discretization::DOFManager> d_DOFMan;
    std::shared_ptr<AMP::Mesh::BoxMesh>              d_BoxMesh;
    // The current time of the solution
    double currentTime = 0.0;

    ManufacturedHeatModel( 
        std::shared_ptr<AMP::Discretization::DOFManager> DOFMan_,
        std::shared_ptr<AMP::Mesh::BoxMesh> BoxMesh_ ) : d_DOFMan(DOFMan_), d_BoxMesh(BoxMesh_) { 

        AMP_INSIST( DOFMan_,  "Non-null DOFManager required!" );
        AMP_INSIST( BoxMesh_, "Non-null BoxMesh required!" );
    }

    // Exact solution, its gradient and corresponding source term
    double exactSolution(double x);
    double exactSolutionGradient(double x);
    double sourceTerm(double x);

    // Populate vectors with exact PDE solution and corresponing source term
    void fillWith_uexact( std::shared_ptr<AMP::LinearAlgebra::Vector> uexact  );
    void fillWith_fsource( std::shared_ptr<AMP::LinearAlgebra::Vector> fsource );

    inline double getCurrentTime() { return currentTime; };
    inline void setCurrentTime(double currentTime_) { currentTime = currentTime_; };
    
    // TODO: Make this a function where the a and b coefficients are specified, and the boundary, and the spatial location...
    inline std::vector<double> getRobinValues( std::shared_ptr<const AMP::Database> db_ ) {
        // Get Robin BC values for manufactured solution 
        // Note that the Robin BCs evaluate the solution exactly at x = 0 and x = 1, even though these are not computational vertices (they are half way between a ghost vertex and the first interior vertex, which ultimately is eliminated from the system)

        // Unpack constants from the database
        double a0 = db_->getScalar<double>( "a0" );
        double b0 = db_->getScalar<double>( "b0" );
        double a1 = db_->getScalar<double>( "a1" );
        double b1 = db_->getScalar<double>( "b1" );
        
        double u0    = exactSolution( 0.0 );
        double u1    = exactSolution( 1.0 );
        double dudx0 = exactSolutionGradient( 0.0 );
        double dudx1 = exactSolutionGradient( 1.0 );
        
        double r0 = a0 * u0 + b0 * dudx0;
        double r1 = a1 * u1 + b1 * dudx1;

        std::vector<double> r = { r0, r1 };
        return r;
    };
};


double ManufacturedHeatModel::exactSolution(double x) {
    double t = this->getCurrentTime();
    return std::sin((3.0/2.0)*M_PI*x)*std::cos(2*M_PI*t);
}

double ManufacturedHeatModel::exactSolutionGradient(double x) {
    double t = this->getCurrentTime();
    return (3.0/2.0)*M_PI*std::cos(2*M_PI*t)*std::cos((3.0/2.0)*M_PI*x);
}

double ManufacturedHeatModel::sourceTerm(double x) {
    double t = this->getCurrentTime();
    //double cxx = d_c.xx;
    double cxx = 1.0; // todo: make consistent
    return (1.0/4.0)*M_PI*(9*M_PI*cxx*std::cos(2*M_PI*t) - 8*std::sin(2*M_PI*t))*std::sin((3.0/2.0)*M_PI*x);
}

/* Populate exact solution vector.  */
void ManufacturedHeatModel::fillWith_uexact( std::shared_ptr<AMP::LinearAlgebra::Vector> uexact ) {

    auto it      = d_BoxMesh->getIterator(d_geomType); // Mesh iterator
    auto meshDim = d_BoxMesh->getDim(); // Dimension

    // Fill in exact solution and source term vectors
    for ( auto elem = it.begin(); elem != it.end(); elem++ ) {
        
        std::vector<double> u;

        if ( meshDim == 1 ) {
            double x = ( elem->coord() )[0];
            u.push_back(exactSolution(x));
        } 

        std::vector<size_t> i;
        d_DOFMan->getDOFs( elem->globalID(), i );
        uexact->setValuesByGlobalID( 1, &i[0], &u[0] );
    }
    uexact->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}


/* Populate source-term vector */
void ManufacturedHeatModel::fillWith_fsource( std::shared_ptr<AMP::LinearAlgebra::Vector> fsource ) {

    auto it      = d_BoxMesh->getIterator(d_geomType); // Mesh iterator
    auto meshDim = d_BoxMesh->getDim(); // Dimension

    // Fill in exact solution and source term vectors
    for ( auto elem = it.begin(); elem != it.end(); elem++ ) {
        
        std::vector<double> f;
        if ( meshDim == 1 ) {
            double x = ( elem->coord() )[0];
            f.push_back(sourceTerm( x ));
        } 

        std::vector<size_t> i;
        d_DOFMan->getDOFs( elem->globalID(), i );
        fsource->setValuesByGlobalID( 1, &i[0], &f[0] );
    }
    fsource->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}


/* ------------------------------------------------
    Class implementing a backward Euler operator 
------------------------------------------------- */
/* Implements the linearOperator I + gamma*L, corresponding to linear systems arising from the application of BDF methods to the linear ODEs u'(t) + L*u = s(t) */ 
class BELinearOp : public AMP::Operator::LinearOperator {

public:
    double                                         d_gamma;
    std::shared_ptr<AMP::Operator::LinearOperator> d_L;

    bool d_matrixRequired; // Does this operator have to construct a matrix?
    bool d_gammaChanged;   // Has the value of gamma been changed since the matrix was last constructed?
    
    // Constructor call's base class's constructor
    BELinearOp( std::shared_ptr<const AMP::Operator::OperatorParameters> params_ ) : 
            AMP::Operator::LinearOperator( params_ ) {

        d_gamma = -1.0; // Initialize with some dummy value
        d_matrixRequired = params_->d_db->getWithDefault<bool>( "matrixRequired", false );

    }

    /* Build the matrix A = I + gamma*L */
    void buildAndSetMatrix() {

        std::shared_ptr<AMP::LinearAlgebra::Matrix> L_mat = d_L->getMatrix();
        AMP_INSIST( L_mat, "Setting the matrix requires L to have a non-null matrix" );

        AMP_INSIST( !(this->getMatrix()), "The case of changing gamma is not implemented yet..." );

        // Should use a clone and then copy if creating a new matrix rather than this copy
        // auto A_mat = AMP::LinearAlgebra::Matrix( *L_mat ); // Copy constructor
        
        auto A_mat = L_mat;
        A_mat->scale( d_gamma ); // A <- gamma*A

        // A <- A + I
        auto DOFMan = A_mat->getRightDOFManager();
        for ( auto dof = DOFMan->beginDOF(); dof != DOFMan->endDOF(); dof++ ) {
            A_mat->addValueByGlobalID(dof, dof, 1.0);
        }
        
        // Set the matrix
        this->setMatrix( A_mat );
    }

    // Set the LinearOperator L
    void setRHSOperator( std::shared_ptr<AMP::Operator::LinearOperator> L_ ) { d_L = L_; };

    std::shared_ptr<AMP::LinearAlgebra::Vector> getRightVector() const override {
        AMP_INSIST( d_L, "RHS operator not set!" );
        return d_L->getRightVector();
    }

    // Used to register this operator in a factory
    std::string type() const override {
        return "BELinearOp";
    }

    // Set the time-step size in the operator
    void setGamma( AMP::Scalar gamma_ ) { 
        double gamma = double( gamma_ );

        d_gammaChanged = !(AMP::Utilities::approx_equal( gamma, d_gamma, 1e-10 ));

        // Build matrix if: 1. gamma has changed, and 2. one is required.
        if (d_gammaChanged) {
            d_gamma = gamma;
            if (d_matrixRequired) {
                buildAndSetMatrix();
            }
        }
    }

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u_in,
                AMP::LinearAlgebra::Vector::shared_ptr r ) override
    {
        AMP_INSIST( d_L, "RHS operator not set!" );

        // Compute r <- (I + gamma*A)*u
        d_L->apply( u_in, r );
        r->axpby (1.0, d_gamma, *u_in); // r <- 1.0*u + d_gamma * r
    }
};


/* ------------------------------------------------
    Class implementing a discrete Poisson problem 
------------------------------------------------- */
class PoissonOp : public AMP::Operator::LinearOperator {

private:
    
    PDECoefficients                     d_c;
    AMP::Mesh::GeomType                 d_geomType = AMP::Mesh::GeomType::Vertex;

    void setPDECoefficients();
    std::shared_ptr<AMP::LinearAlgebra::Matrix> getLaplacianMatrix();
    //std::shared_ptr<AMP::LinearAlgebra::Matrix> d_APoisson;
    
    double d  = 1.0; // diffusion coefficient... make consistent with c.xx
    
    // Constants in Robin BCs
    double d_a0, d_b0, d_r0;
    double d_a1, d_b1, d_r1;
    
public:

    /* Set Robin values */
    void resetRobinValues( double r0_, double r1_ ) { d_r0 = r0_; d_r1 = r1_; };

    void ApplyRobinCorrectionToMatrix( std::shared_ptr<AMP::LinearAlgebra::Matrix> A );
    void ApplyRobinCorrectionToVector( std::shared_ptr<AMP::LinearAlgebra::Vector> f );

    std::shared_ptr<AMP::Database>                   d_db;
    std::shared_ptr<AMP::Discretization::DOFManager> d_DOFMan;
    std::shared_ptr<AMP::Mesh::BoxMesh>              d_BoxMesh;

    
    // Constructor call's base class's constructor
    PoissonOp(std::shared_ptr<const AMP::Operator::OperatorParameters> params_) : 
            AMP::Operator::LinearOperator( params_ ) { 

        // Keep a pointer to my BoxMesh to save having to do this downcast repeatedly 
        d_BoxMesh = std::dynamic_pointer_cast<AMP::Mesh::BoxMesh>(this->getMesh());
        AMP_INSIST( d_BoxMesh, "Mesh must be a AMP::Mesh::BoxMesh" );

        // Set my database
        d_db = params_->d_db->getDatabase( "PDE_db" );

        // Unpack constants from the database
        d_a0 = d_db->getScalar<double>( "a0" );
        d_b0 = d_db->getScalar<double>( "b0" );
        d_r0 = d_db->getScalar<double>( "r0" );
        d_a1 = d_db->getScalar<double>( "a1" );
        d_b1 = d_db->getScalar<double>( "b1" );
        d_r1 = d_db->getScalar<double>( "r0" );

        // Set DOFManager
        this->set_DOFManager();
        AMP_INSIST(  d_DOFMan, "Requires non-null DOF" );

        // Set PDE coefficients
        this->setPDECoefficients();

        // Get the matrix
        auto A = getLaplacianMatrix( );

        // Update according to Robin data
        ApplyRobinCorrectionToMatrix( A );
        
        // Set linear operator's matrix
        this->setMatrix( A );
        //d_APoisson = A;
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


private:

    /* Build and set d_DOFMan */
    void set_DOFManager() {
        int DOFsPerElement = 1; 
        int gcw  = 1; // Ghost-cell width (stencils are at most 3-point in each direction)
        d_DOFMan = AMP::Discretization::boxMeshDOFManager::create(this->getMesh(), d_geomType, gcw, DOFsPerElement);
    }

    // Coefficients defining the PDE
    std::vector<double> getPDECoefficients1D();
    
    // FD stencils
    std::vector<double> getStencil1D();

    // Map from row index to col+data in that row
    std::map<size_t, colsDataPair> getCSRData1D();

    // Map from grid index i, or i,j, or i,j,k to a MeshElementIndex to a MeshElementId and then to the corresponding DOF
    size_t gridIndsToDOF( int i, int j = 0, int k = 0 ) {
        AMP::Mesh::BoxMesh::MeshElementIndex ind(
                        AMP::Mesh::GeomType::Vertex, 0, i, j, k );
        AMP::Mesh::MeshElementID id = d_BoxMesh->convert( ind );
        std::vector<size_t> dof;
        d_DOFMan->getDOFs(id, dof);
        return dof[0];
    };
}; 


void PoissonOp::setPDECoefficients() {
    int dim = this->getMesh()->getDim();
    if ( dim == 1 ) {
        d_c = PDECoefficients( 1, getPDECoefficients1D() );
    }
}

std::vector<double> PoissonOp::getPDECoefficients1D() {
    // PDE coefficients
    double cxx = 1.0;
    std::vector<double> c = { cxx };
    return c;
}

/* Get 3-point stencil for 1D Poisson. */
std::vector<double> PoissonOp::getStencil1D() {
    // Unpack PDE coefficients
    double cxx = d_c.xx;
    // Stencil coefficients
    double W  = -1.0*cxx; 
    double O  = +2.0*cxx; 
    double E  = -1.0*cxx;
    // Populate stencil
    std::vector<double> stencil = { W, O, E };

    // Introduce 1/h^2 scaling 
    auto h = d_db->getScalar<double>("h");
    for ( auto &s : stencil ) {
        s *= 1.0/(h*h);
    }

    return stencil;
}

/* Get CSR structure of 1D Laplacian */
std::map<size_t, colsDataPair> PoissonOp::getCSRData1D() {    

    // Get 3-point stencil
    auto stencil = getStencil1D( );

    // Get local grid index box w/ zero ghosts
    auto localBox  = getLocalNodeBox( d_BoxMesh );
    auto globalBox = getGlobalNodeBox( d_BoxMesh );

    // Create a map from the DOF to a pair a vectors
    // Map from a DOF to vector of col inds and associated data
    std::map<size_t, colsDataPair> localCSRData;

    // Iterate over local interior box
    for (auto i = localBox.first[0]; i <= localBox.last[0]; i++) {
        
        // Current row
        auto dof = gridIndsToDOF( i );

        // West boundary
        if (i == globalBox.first[0]) {
            std::vector<size_t> cols = { dof, gridIndsToDOF( i+1 )};
            std::vector<double> vals = { stencil[1], stencil[2] };
            localCSRData[dof] = { cols, vals };
            continue;
        }
        // East boundary
        if (i == globalBox.last[0]) {
            std::vector<size_t> cols = { gridIndsToDOF( i-1 ), dof };
            std::vector<double> vals = { stencil[0], stencil[1] };
            localCSRData[dof] = { cols, vals };
            continue;
        }

        // Copy of stencil
        std::vector<double> vals = stencil; 
        // Column indices, ordered consistently with the stencil
        std::vector<size_t> cols = { 
            gridIndsToDOF( i-1 ),
            dof,
            gridIndsToDOF( i+1 ) };
        localCSRData[dof] = { cols, vals };
    }  

    return localCSRData;
}


/* This correction only depends on the constants a_i and b_i; it does not depend on the constants r_i!  

Robin set up. Use cell centered discretization, with one ghost cell. The stencil for the first interior point uses the value in the ghost cell. The discretized boundary condition is solved for the ghost value as a function of the first interior DOF, which is then substituted into the stencil for the first interior DOF.

at x = 0, the BC is
a0 * E - b0 * d_E * dE/dx = r0

at x = 1, the BC is
a1 * E - b1 * d_E * dE/dx = r1.
*/
void PoissonOp::ApplyRobinCorrectionToMatrix( std::shared_ptr<AMP::LinearAlgebra::Matrix> A ) {

    auto localBox  = getLocalNodeBox( d_BoxMesh );
    auto globalBox = getGlobalNodeBox( d_BoxMesh );
    auto h = d_db->getScalar<double>("h");
    double d = d_c.xx;

    // West boundary
    if ( localBox.first[0] == globalBox.first[0] ) {

        int i = localBox.first[0];
        size_t dof = gridIndsToDOF( i );

        double alpha0 = (2*d*d_b0 + d_a0*h)/(2*d*d_b0 - d_a0*h);

        // We solve for the ghost value as:
        // Eg = alpha0*E0 + beta0
        // 1/h^2*( - Eg + 2*E0 - E1 ) = 1/h^2*( [2-alpha0]*E0 - E1 ) - 1/h^2 * beta0
        std::vector<size_t> cols = { dof, gridIndsToDOF( i+1 ) };
        std::vector<double> vals = { 1.0/(h*h) * ( 2.0 - alpha0 ), -1.0/(h*h) };

        // Set values by grid indices
        int nrows = 1;
        int ncols = 2;
        A->setValuesByGlobalID<double>( nrows, ncols, &dof, cols.data(), vals.data() );
    }

    // East boundary
    if ( localBox.last[0] == globalBox.last[0] ) {
        int i = localBox.last[0];
        size_t dof  = gridIndsToDOF( i );

        double alpha1 = (2*d*d_b1 - d_a1*h)/(2*d*d_b1 + d_a1*h);

        // We solve for the ghost value as:
        // Eg = alpha1*E1 + beta1
        // 1/h^2*( - E0 + 2*E1 - Eg ) = 1/h^2*( - E0 + [2-alpha1]*E1 ) - 1/h^2 * beta1
        std::vector<size_t> cols = { gridIndsToDOF( i-1 ), dof };
        std::vector<double> vals = { -1.0/(h*h), 1.0/(h*h) * ( 2.0 - alpha1 ) };

        // Set values by grid indices
        int nrows = 1;
        int ncols = 2;
        A->setValuesByGlobalID<double>( nrows, ncols, &dof, cols.data(), vals.data() );
    }

    A->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}

/* Applies corrections to vector f that arise from the elimination of the first set of interior vertices. This correction depends on all of the constants a_i, b_i, and r_i  */
void PoissonOp::ApplyRobinCorrectionToVector( std::shared_ptr<AMP::LinearAlgebra::Vector> f ) {

    AMP_INSIST( f, "Non-null f required!" );

    auto localBox  = getLocalNodeBox( d_BoxMesh );
    auto globalBox = getGlobalNodeBox( d_BoxMesh );
    auto h = d_db->getScalar<double>("h");
    double d = d_c.xx;

    // West boundary
    if ( localBox.first[0] == globalBox.first[0] ) {
        int i             = localBox.first[0];
        int dof           = gridIndsToDOF( i );
        double beta0      = -2*d_r0*h/(2*d*d_b0 - d_a0*h);
        double correction = 1.0/(h*h) * beta0;
        
        // Add correction to existing value of f
        double f0 = f->getValueByGlobalID( dof ) + correction;
        f->setValueByGlobalID<double>( dof, f0 );
    }

    // East boundary
    if ( localBox.last[0] == globalBox.last[0] ) {
        int i             = localBox.last[0];
        int dof           = gridIndsToDOF( i );
        double beta1      = 2*d_r1*h/(2*d*d_b1 + d_a1*h);
        double correction = 1.0/(h*h) * beta1;
        
        // Add correction to existing value of f
        double f1 = f->getValueByGlobalID( dof ) + correction;
        f->setValueByGlobalID<double>( dof, f1 );
    }

    f->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}

/* Return a constructed CSR matrix corresponding to the discretized Laplacian on the mesh */
std::shared_ptr<AMP::LinearAlgebra::Matrix> PoissonOp::getLaplacianMatrix( ) {

    auto inVec  = this->getRightVector();
    auto outVec = this->getRightVector();

    int meshDim = this->getMesh()->getDim();
    std::map<size_t, colsDataPair> localCSRData;
    if ( meshDim == 1 ) {
        localCSRData = this->getCSRData1D( );
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
            AMP::UnitTest *ut,
            const std::string &exeName ) {

    std::string input_file = "input_"  + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::logOnlyNodeZero( log_file );

    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Unpack databases from the input file
    auto PDE_db    = input_db->getDatabase( "PDE" );
    //auto solver_db = input_db->getDatabase( "LinearSolver" );
    auto ti_db     = input_db->getDatabase( "TimeIntegrator" );

    AMP_INSIST( PDE_db,    "A PDE database must be provided" );
    //AMP_INSIST( solver_db, "A LinearSolver database must be provided" );
    AMP_INSIST( ti_db, "A TimeIntegrator database must be provided" );


    /****************************************************************
    * Create a mesh                                                 *
    ****************************************************************/
    static std::shared_ptr<AMP::Mesh::BoxMesh> mesh = createBoxMesh( comm, PDE_db );

    // Print basic problem information
    AMP::pout << "--------------------------------------------------------------------------------" << std::endl;
    AMP::pout << "Solving " << static_cast<int>(mesh->getDim()) 
        << "D Poisson problem on mesh with " 
        << mesh->numGlobalElements(AMP::Mesh::GeomType::Vertex) 
        << " (==" << PDE_db->getScalar<int>( "n" )+1 << "^" << static_cast<int>(mesh->getDim()) << ")"
        << " total DOFs across " << mesh->getComm().getSize() 
        << " ranks" << std::endl;
    AMP::pout << "--------------------------------------------------------------------------------" << std::endl;


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

    // Create PoissonOp 
    auto myPoissonOp = std::make_shared<PoissonOp>( OpParameters );    
    // Create BEOp based on Poisson Op
    auto myBEOp      = std::make_shared<BELinearOp>( OpParameters );  
    myBEOp->setRHSOperator( myPoissonOp );
    

    /****************************************************************
    * Set up maunfactured solution                                  *
    ****************************************************************/
    // Create required vectors over the mesh
    auto myManufacturedHeat = std::make_shared<ManufacturedHeatModel>( myPoissonOp->d_DOFMan, myPoissonOp->d_BoxMesh );

    /****************************************************************
    * Set up relevant vectors over the mesh                         *
    ****************************************************************/
    // Create required vectors over the mesh
    auto numSolVec    = myPoissonOp->getRightVector();
    auto manSolVec    = myPoissonOp->getRightVector();
    auto errorVec     = myPoissonOp->getRightVector();
    auto BDFSourceVec = myPoissonOp->getRightVector();


    /****************************************************************
    * Set up implicit time integrator                               *
    ****************************************************************/
    // Parameters for time integrator
    double dt = PDE_db->getScalar<double>( "h" );

    // Create initial condition vector
    auto ic = myPoissonOp->getRightVector();
    myManufacturedHeat->setCurrentTime( 0.0 );
    myManufacturedHeat->fillWith_uexact( ic );
    ic->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );

    // Create vectors to hold current and new solutions
    auto sol_old = ic->clone( );
    sol_old->copyVector( ic ); 
    auto sol_new = ic->clone( );
    sol_new->copyVector( ic ); 
    
    auto tiParams           = std::make_shared<AMP::TimeIntegrator::TimeIntegratorParameters>( ti_db );
    tiParams->d_ic_vector   = ic;
    tiParams->d_operator    = myBEOp;
    // Create a source vector
    tiParams->d_pSourceTerm = BDFSourceVec;

    // Create timeIntegrator from factory
    std::shared_ptr<AMP::TimeIntegrator::TimeIntegrator> timeIntegrator = AMP::TimeIntegrator::TimeIntegratorFactory::create( tiParams );
    
    auto implicitIntegrator =
            std::dynamic_pointer_cast<AMP::TimeIntegrator::ImplicitIntegrator>( timeIntegrator );

    // required to set gamma for the user operator at each timestep
    implicitIntegrator->setTimeScalingFunction(
        std::bind( &BELinearOp::setGamma, &( *myBEOp ), std::placeholders::_1 ) );

    #if 0
    int step = 0;
    int n = PDE_db->getScalar<int>( "n" );
    std::string out_dir = "out/n" + std::to_string(n) + "/";
    std::string num_dir = out_dir + "unum";
    std::string man_dir = out_dir + "uman";    
    
    // Remove outdir and its contents if it already exists
    if (std::filesystem::is_directory(out_dir)) {
        std::filesystem::remove_all( out_dir ); 
    }
    //create outdir
    std::filesystem::create_directory( out_dir );
    
    // Write IC
    {
        double T = 0.0;
        AMP::IO::AsciiWriter vecWriter_man;
        std::string name = std::to_string( T );
        manSolVec->setName( name );
        sol_new->setName( name );
        vecWriter_man.registerVector( ic );
        vecWriter_man.writeFile( man_dir, step, T  );
        AMP::IO::AsciiWriter vecWriter_num;
        vecWriter_num.registerVector( sol_new );
        vecWriter_num.writeFile( num_dir, step, T  );
    }
    #endif


    // Integrate!
    double finalTime = timeIntegrator->getFinalTime();
    double T = 0.0;
    timeIntegrator->setInitialDt( dt );
    while ( T < finalTime ) {

        // Set the solution-independent source term; note that this approach only works for implicit multistep methods
        myManufacturedHeat->setCurrentTime( T + dt ); // Set source to new time...
        // Fill vector with PDE source term
        myManufacturedHeat->fillWith_fsource( BDFSourceVec );
        // Get exact Robin coefficients at T + dt
        auto robinCoeffs = myManufacturedHeat->getRobinValues( PDE_db );
        // Reset Robin coefficients in Poisson operator
        myPoissonOp->resetRobinValues( robinCoeffs[0], robinCoeffs[1] );
        // Apply Robin correction to source vector
        myPoissonOp->ApplyRobinCorrectionToVector( BDFSourceVec );

        // Advance the solution
        timeIntegrator->advanceSolution( dt, T == 0, sol_old, sol_new );
        if ( timeIntegrator->checkNewSolution() ) {
            timeIntegrator->updateSolution();
            sol_old->copyVector( sol_new );
        } else {
            AMP_ERROR( "Solution didn't converge" );
        }

        // Update time
        T += dt;

        /* Compare numerical solution with manufactured solution */
        myManufacturedHeat->setCurrentTime( T );
        myManufacturedHeat->fillWith_uexact( manSolVec );
        errorVec->subtract( *sol_new, *manSolVec );
        AMP::pout << "Discretization error norms:" << std::endl;
        auto enorms = getDiscreteNorms( PDE_db->getScalar<double>( "h" ), errorVec );
        AMP::pout << "||e||=(" << enorms[0] << "," << enorms[1] << "," << enorms[2] << ")" << std::endl;
        AMP::pout << "----------------------" << std::endl;

        
        // Write manufactured and numerical solution to file.
        #if 0
        step++;
        std::string name = std::to_string( T );
        manSolVec->setName( name );
        sol_new->setName( name );
        AMP::IO::AsciiWriter vecWriter_man;
        vecWriter_man.registerVector( manSolVec );
        vecWriter_man.writeFile( man_dir, step, T+dt  );
        AMP::IO::AsciiWriter vecWriter_num;
        vecWriter_num.registerVector( sol_new );
        vecWriter_num.writeFile( num_dir, step, T+dt  );
        #endif
        

        // Drop out if we've exceeded max steps
        if ( !timeIntegrator->stepsRemaining() ) {
            AMP_WARNING( "max_integrator_steps has been reached, dropping out of loop now..." );
            break;
        }
    }

    
    // // No specific solution is implemented for this problem, so this will just check that the solver converged. 
    // checkConvergence( linearSolver.get(), input_db, input_file, *ut );

}
// end of driver()



/*  The input file must contain a "PDE" database, which has 

    dim : the dimension of the problem (1, 2, or 3)
    n   : the number of mesh points (plus 1) in each grid dimension
    
    - In 1D, the PDE is -u_xx = f. Standard 3-point finite differences are used. 
    - In 2D, the PDE is -u_xx -eps*u_yy = f, but rotated an angle of theta radians counter clockwise from the positive x axis. 9-point finite differences are used.
    - In 3D, the PDE is -u_xx -epsy*u_yy - epsz*u_zz = f, but rotated according in 3D according the the so-called "extrinsic" 3-1-3 Euler angles (see https://en.wikipedia.org/wiki/Euler_angles#Definition_by_extrinsic_rotations) gamma, beta, and alpha. A new coordinate system XYZ is obtained from rotating the xyz coordinate system by: 1. first rotating gamma radians around the z axis, 2. then beta radians about the x axis, 3. then alpha radians about the z axis. 19-point finite differences are used. 
    
    - In 2D, the PDE database must also specify the constants "eps", and "theta"
    - In 2D, the PDE database must also specify the constants "epsy", "epsz", "gamma", "beta", and "alpha"

    - Mixed derivatives are discretized with central finite differences (as opposed to upwind)

    - Note that for eps = epsy = epsz = 1, the PDEs are isotropic and grid aligned (the rotation angles have no impact in this case), and moreover, in these cases the discretizations reduce to the standard 5- and 7-point stencils in 2D and 3D, respectively.

*/
int main( int argc, char **argv )
{

    AMP::AMPManager::startup( argc, argv );      

    AMP::UnitTest ut;

    // Create a global communicator
    AMP::AMP_MPI comm( AMP_COMM_WORLD );
    
    std::vector<std::string> exeNames;
    exeNames.emplace_back( "heat" );

    for ( auto &exeName : exeNames ) {
        driver( comm, &ut, exeName );
    }    
    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}


// /* Get 7-point stencil for 3D Poisson. */
// std::vector<double> PoissonOp::getStencil3D() {
    
//     auto epsy = d_db->getScalar<double>( "epsy" );
//     auto epsz = d_db->getScalar<double>( "epsz" );

//     double D  = -1.0*epsz;
//     double S  = -1.0*epsy;
//     double W  = -1.0; 
//     double O  = +2.0 + 2.0*epsy + 2.0*epsz; 
//     double E  = -1.0;
//     double N  = -1.0*epsy;
//     double U  = -1.0*epsz;

//     // Populate stencil
//     std::vector<double> stencil = { O, D, S, W, E, N, U };

//     // Introduce 1/h^2 scaling 
//     auto h = d_db->getScalar<double>( "h" );
//     for ( auto &s : stencil ) {
//         s *= 1.0/(h*h);
//     }

//     return stencil;
// }


// /* Get 9-point stencil for upwind FD disretization of rotated anisotropic 2D Poisson. */
// std::vector<double> PoissonOp::getStencil2D() {
    
//     double eps   = this->d_db->getScalar<double>( "eps" );
//     double theta = this->d_db->getScalar<double>( "theta" );

//     AMP_INSIST( theta >= 0.0 && theta <= M_PI/2.0, "Upwind discretization only valid for theta in [0,pi/2]" );

//     double c     = cos(theta);
//     double s     = sin(theta);

//     // PDE coefficients
//     double alpha = c*c + eps * s*s;
//     double beta  = eps * c*c + s*s;
//     double gamma = 2*(1 - eps) * c * s;
//     gamma *= 0.5; // gamma only ever appears multiplied by 1/2, except at O.
//     double NW = 0.0;
//     double N  =           -   beta +   gamma; 
//     double NE =                    -   gamma; 
//     double W  = -   alpha          +   gamma; 
//     double O  = + 2*alpha + 2*beta - 2*gamma; 
//     double E  = -   alpha          +   gamma;
//     double SW =                    -   gamma;
//     double S  =           -   beta +   gamma;
//     double SE = 0.0; 

//     // Populate stencil
//     std::vector<double> stencil = { O, SW, S, SE, W, E, NW, N, NE };

//     // Introduce 1/h^2 scaling 
//     auto h = d_db->getScalar<double>( "h" );
//     for ( auto &s : stencil ) {
//         s *= 1.0/(h*h);
//     }

//     return stencil;
// }
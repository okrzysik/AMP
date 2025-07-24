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

/* Numerical implementation of the heat equation with finite difference spatial discretization and Robin boundary conditions.

    Some of the code is generated in the associated -helper.py file.
*/


/* Representation of constant coefficients in front of derivatives in second-order differential operator
In 1D: cxx*u_xx
In 2D: cxx*u_xx + cyy*u_yy 
*/
struct PDECoefficients {
    double xx = std::nan("");
    double yy = std::nan("");

    PDECoefficients() {};
    PDECoefficients( int dim, std::vector<double> c ) {
        if ( dim == 1 ) {
            xx = c[0];
        } else if ( dim == 2 ) {
            xx = c[0], yy = c[1];
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
Since there are n cells in [0, 1], the width of each must be 1/n
The computational vertices are the mid points of each cell, and we have n active vertices:
        0.5*h, h + 0.5*h, ..., 1 - (h + 0.5*), 1 - 0.5*h. 

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
        mesh_db->putVector<int>( "Size", { n-1 } ); // mesh has n points
        mesh_db->putVector<double>( "Range", { 0.5*h, 1.0 - 0.5*h } );
    } else if ( dim == 2 ) {
        mesh_db->putVector<int>( "Size", { n-1, n-1 } ); // mesh has n x n points
        mesh_db->putVector<double>( "Range", { 0.5*h, 1.0 - 0.5*h, 0.5*h, 1.0 - 0.5*h } );
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
    #if 0
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
    Class implementing a backward Euler operator 
------------------------------------------------- */
/* Implements the linearOperator I + gamma*L, corresponding to linear systems arising from the application of BDF methods to the linear ODEs u'(t) + L*u = s(t) */ 
class BELinearOp : public AMP::Operator::LinearOperator {

public:
    double                                         d_gamma;
    std::shared_ptr<AMP::Operator::LinearOperator> d_L;

    bool d_matrixRequired; // Does this operator have to construct a matrix?
    bool d_gammaChanged;   // Has the value of gamma been changed since the matrix was last constructed?
    int d_print_info_level;
    
    // Constructor call's base class's constructor
    BELinearOp( std::shared_ptr<const AMP::Operator::OperatorParameters> params_ ) : 
            AMP::Operator::LinearOperator( params_ ) {

        d_gamma = -1.0; // Initialize with dummy value
        d_matrixRequired   = params_->d_db->getScalar<bool>( "matrixRequired" );
        d_print_info_level = params_->d_db->getWithDefault<int>( "print_info_level", 0 );
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
        if ( d_print_info_level > 0 )
            AMP::pout << "BELinearOp::setGamma() ";

        double gamma = double( gamma_ );
        d_gammaChanged = !(AMP::Utilities::approx_equal( gamma, d_gamma, 1e-10 ));

        // Build matrix if: 1. gamma has changed, and 2. one is required.
        if (d_gammaChanged) {
            if ( d_print_info_level > 0 )
                AMP::pout << "gamma has changed: gamma_old=" << d_gamma << ", gamma_new=" << gamma << std::endl;
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

        // AMP::Operator::LinearOperator::apply( u_in, r );
    }

private:
    /* Build the matrix A = I + gamma*L. Note that we create A from a copy of the matrix L; this is the simplest way to handle if gamma has changed... */
    void buildAndSetMatrix() {

        if ( d_print_info_level > 0 )
            AMP::pout << "BELinearOp::buildAndSetMatrix()" << std::endl;

        std::shared_ptr<AMP::LinearAlgebra::Matrix> L_mat = d_L->getMatrix();
        AMP_INSIST( L_mat, "Setting the matrix requires L to have a non-null matrix" );

        //AMP_INSIST( !(this->getMatrix()), "The case of changing gamma is not implemented yet..." );
        //auto A_mat = L_mat;
        
        // A <- L
        auto A_mat = L_mat->clone();
        A_mat->copy( L_mat );

        // A <- gamma*A
        A_mat->scale( d_gamma ); 

        // A <- A + I
        auto DOFMan = A_mat->getRightDOFManager(  );
        for ( auto dof = DOFMan->beginDOF(); dof != DOFMan->endDOF(); dof++ ) {
            A_mat->addValueByGlobalID(dof, dof, 1.0);
        }
        A_mat->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
        
        // Set the matrix
        this->setMatrix( A_mat );

        AMP_INSIST( this->getMatrix(), "Why is matrix null?" );
    }
};


/* ------------------------------------------------
    Class implementing a discrete Poisson problem 
------------------------------------------------- */
class PoissonOp : public AMP::Operator::LinearOperator {

public:
    
    AMP::Mesh::GeomType                              d_geomType = AMP::Mesh::GeomType::Vertex;
    PDECoefficients                                  d_c;
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
        d_a1 = d_db->getScalar<double>( "a1" );
        d_b1 = d_db->getScalar<double>( "b1" );
        d_r1 = d_db->getScalar<double>( "r1" );
        d_a2 = d_db->getScalar<double>( "a2" );
        d_b2 = d_db->getScalar<double>( "b2" );
        d_r2 = d_db->getScalar<double>( "r2" );
        //
        if ( d_BoxMesh->getDim() >= 2 ) {
            d_a3 = d_db->getScalar<double>( "a3" );
            d_b3 = d_db->getScalar<double>( "b3" );
            d_r3 = d_db->getScalar<double>( "r3" );
            d_a4 = d_db->getScalar<double>( "a4" );
            d_b4 = d_db->getScalar<double>( "b4" );
            d_r4 = d_db->getScalar<double>( "r4" );
        }

        // Specify Robin return function as the default
        std::function<double( int, double, double, AMP::Mesh::MeshElement & )> wrapper = [&]( int boundary, double, double, AMP::Mesh::MeshElement & ) { return robinFunctionDefault( boundary ); };
        this->setRobinFunction( wrapper );
        
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

    // The user can specify any Robin return function with this signature; if they do then this will overwrite the default.
    void setRobinFunction( std::function<double(int, double, double, AMP::Mesh::MeshElement &)> fn_ ) { d_robinFunction = fn_; };

    // Populate a vector with the given function
    void fillWithFunction( std::shared_ptr<AMP::LinearAlgebra::Vector> u, std::function<double(AMP::Mesh::MeshElement &)> fun );

    void ApplyRobinCorrectionToVector( std::shared_ptr<AMP::LinearAlgebra::Vector> f );

private:

    // Constants in Robin BCs
    double d_a1, d_b1, d_r1; // West boundary,  x = 0
    double d_a2, d_b2, d_r2; // East boundary,  x = 1
    double d_a3, d_b3, d_r3; // South boundary, y = 0
    double d_a4, d_b4, d_r4; // North boundary, y = 1
      
    std::shared_ptr<AMP::LinearAlgebra::Matrix> getLaplacianMatrix();
  
    /* Build and set d_DOFMan */
    void set_DOFManager() {
        int DOFsPerElement = 1; 
        int gcw  = 1; // Ghost-cell width (stencils are at most 3-point in each direction)
        d_DOFMan = AMP::Discretization::boxMeshDOFManager::create(this->getMesh(), d_geomType, gcw, DOFsPerElement);
    }


    void setPDECoefficients();

    // Coefficients defining the PDE
    std::vector<double> getPDECoefficients1D();
    std::vector<double> getPDECoefficients2D();
    
    // FD stencils
    std::vector<double> getStencil1D();
    std::vector<double> getStencil2D();

    // Map from row index to col+data in that row
    std::map<size_t, colsDataPair> getCSRData1D();
    std::map<size_t, colsDataPair> getCSRData2D();

    // Map from grid index i, or i,j, or i,j,k to a MeshElementIndex to a MeshElementId and then to the corresponding DOF
    size_t gridIndsToDOF( int i, int j = 0, int k = 0 ) {
        AMP::Mesh::BoxMesh::MeshElementIndex ind(
                        AMP::Mesh::GeomType::Vertex, 0, i, j, k );
        AMP::Mesh::MeshElementID id = d_BoxMesh->convert( ind );
        std::vector<size_t> dof;
        d_DOFMan->getDOFs(id, dof);
        return dof[0];
    };

    // Relating to Robin BCs

    // Prototype of function returning value of Robin BC on given boundary at given node. The user can specify any function with this signature
    std::function<double( int boundary, double a, double b, AMP::Mesh::MeshElement & node )> d_robinFunction;

    // Default function for returning Robin values; this can be overridden by the user
    double robinFunctionDefault( int boundary ) {
        if ( boundary == 1 ) {
            return d_r1;
        } else if ( boundary == 2 ) {
            return d_r2;
        } else if ( boundary == 3 ) {
            return d_r3;
        } else if ( boundary == 4 ) {
            return d_r4;
        } else { 
            AMP_ERROR( "Invalid boundary" );
        }
    }

    void ApplyRobinCorrectionToMatrix( std::shared_ptr<AMP::LinearAlgebra::Matrix> A );
    double RobinCorrectionMatrixCoefficient( int boundary );
    double RobinCorrectionVectorCoefficient( int boundary, AMP::Mesh::MeshElement & node );
}; 


void PoissonOp::setPDECoefficients() {
    int dim = this->getMesh()->getDim();
    if ( dim == 1 ) {
        d_c = PDECoefficients( 1, getPDECoefficients1D() );
    } else if ( dim == 2 ) {
        d_c = PDECoefficients( 2, getPDECoefficients2D() );
    }
}

std::vector<double> PoissonOp::getPDECoefficients1D() {
    // PDE coefficients
    double cxx = 1.0;
    std::vector<double> c = { cxx };
    return c;
}

std::vector<double> PoissonOp::getPDECoefficients2D() {
    // PDE coefficients
    double cxx = 1.0;
    double cyy = 1.0;
    //double cyy = 0.0;
    std::vector<double> c = { cxx, cyy };
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


/* Get 9-point stencil for 2D Poisson that discretizes the operator
    -cxx*u_xx
    -cyy*u_yy

Standard 3-point differences are used for the _ii terms.
*/
std::vector<double> PoissonOp::getStencil2D() {
    
    // Unpack PDE coefficients
    double cxx = d_c.xx;
    double cyy = d_c.yy;

    double O = 0.0;
    // -cxx*u_xx
    double W  = -1.0*cxx; 
           O += +2.0*cxx;  
    double E  = -1.0*cxx;
    // -cyy =*u_yy
    double S  = -1.0*cyy;
           O += +2.0*cyy;
    double N  = -1.0*cyy;

    // Populate stencil
    std::vector<double> stencil = { S, W, O, E, N };

    // Introduce 1/h^2 scaling 
    auto h = d_db->getScalar<double>( "h" );
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

        // In boundary rows, just set interior stencil connections
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

/* Get CSR structure of rotated anisotropic 2D Laplacian */
std::map<size_t, colsDataPair> PoissonOp::getCSRData2D() {    

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
            size_t dof = gridIndsToDOF( i, j );

            // At a boundary DOF Set connections to interior DOFs in boundary rows
            if (j == globalBox.first[1] || j == globalBox.last[1] || i == globalBox.first[0] || i == globalBox.last[0]) {

                size_t DOF_O = dof; // O
                // Unpack stencil connections for readability
                double val_S = stencil[0];
                double val_W = stencil[1];
                double val_O = stencil[2];
                double val_E = stencil[3];
                double val_N = stencil[4];
                
                // WEST boundary. No connection to WEST
                if ( i == globalBox.first[0] ) {
                    size_t DOF_E = gridIndsToDOF( i+1, j   ); // E
                    
                    // WEST boundary, SW corner. No connection to SOUTH or WEST
                    if ( j == globalBox.first[1] ) {
                        size_t DOF_N = gridIndsToDOF( i ,  j+1 ); // N
                        std::vector<size_t> cols = { DOF_O, DOF_E, DOF_N };
                        std::vector<double> vals = { val_O, val_E, val_N };
                        localCSRData[dof] = { cols, vals };

                    // WEST boundary, NW corner. No connection to NORTH or WEST
                    } else if ( j == globalBox.last[1] ) {
                        size_t DOF_S = gridIndsToDOF( i ,  j-1 ); // S
                        std::vector<size_t> cols = { DOF_S, DOF_O, DOF_E };
                        std::vector<double> vals = { val_S, val_O, val_E };
                        localCSRData[dof] = { cols, vals };
                    
                    // WEST boundary, non-corner DOF. No connection to WEST
                    } else {
                        size_t DOF_S = gridIndsToDOF( i ,  j-1 ); // S
                        size_t DOF_N = gridIndsToDOF( i ,  j+1 ); // N
                        std::vector<size_t> cols = { DOF_S, DOF_O, DOF_E, DOF_N };
                        std::vector<double> vals = { val_S, val_O, val_E, val_N };
                        localCSRData[dof] = { cols, vals };
                    }

                // EAST boundary. No connection to EAST
                } else if ( i == globalBox.last[0] ) {
                    size_t DOF_W = gridIndsToDOF( i-1, j   ); // W
                    
                    // EAST boundary, SE corner. No connection to SOUTH or EAST
                    if ( j == globalBox.first[1] ) {
                        size_t DOF_N = gridIndsToDOF( i ,  j+1 ); // N
                        std::vector<size_t> cols = { DOF_W, DOF_O, DOF_N };
                        std::vector<double> vals = { val_W, val_O, val_N };
                        localCSRData[dof] = { cols, vals };

                    // EAST boundary, NE corner. No connection to NORTH or EAST
                    } else if ( j == globalBox.last[1] ) {
                        size_t DOF_S = gridIndsToDOF( i ,  j-1 ); // S
                        std::vector<size_t> cols = { DOF_S, DOF_W, DOF_O };
                        std::vector<double> vals = { val_S, val_W, val_O };
                        localCSRData[dof] = { cols, vals };

                    // EAST boundary, non-corner DOF. No connection to EAST 
                    } else {
                        size_t DOF_S = gridIndsToDOF( i ,  j-1 ); // S
                        size_t DOF_N = gridIndsToDOF( i ,  j+1 ); // N
                        std::vector<size_t> cols = { DOF_S, DOF_W, DOF_O, DOF_N };
                        std::vector<double> vals = { val_S, val_W, val_O, val_N };
                        localCSRData[dof] = { cols, vals };
                    }

                // SOUTH boundary, non-corner DOF. No connection to SOUTH
                } else if ( j == globalBox.first[1] ) {
                    size_t DOF_W = gridIndsToDOF( i-1, j   ); // W
                    size_t DOF_O = dof;                       // O
                    size_t DOF_E = gridIndsToDOF( i+1, j   ); // E
                    size_t DOF_N = gridIndsToDOF( i ,  j+1 ); // N

                    std::vector<size_t> cols = { DOF_W, DOF_O, DOF_E, DOF_N };
                    std::vector<double> vals = { val_W, val_O, val_E, val_N };
                    localCSRData[dof] = { cols, vals };
                    
                // NORTH boundary, non-corner DOF. No connection to NORTH
                } else if ( j == globalBox.last[1] ) {
                    size_t DOF_S = gridIndsToDOF( i ,  j-1 ); // S
                    size_t DOF_W = gridIndsToDOF( i-1, j   ); // W
                    size_t DOF_O = dof;                       // O
                    size_t DOF_E = gridIndsToDOF( i+1, j   ); // E

                    std::vector<size_t> cols = { DOF_S, DOF_W, DOF_O, DOF_E };
                    std::vector<double> vals = { val_S, val_W, val_O, val_E };
                    localCSRData[dof] = { cols, vals };
                }
            
                // Dealt with boundary DOF, move to next one    
                continue;
            }
            
            // Copy of stencil
            std::vector<double> vals = stencil; 
            // Column indices, ordered consistently with the stencil
            std::vector<size_t> cols = { 
                gridIndsToDOF( i ,  j-1 ), // S
                gridIndsToDOF( i-1, j   ), // W
                dof,                       // O
                gridIndsToDOF( i+1, j   ), // E
                gridIndsToDOF( i ,  j+1 ), // N
            };

            localCSRData[dof] = { cols, vals };
        }
    }   

    return localCSRData;
}


/* Equation for u_{ij} is: 1/h^2 * (-u_{ghost} + c*u_{ij} + ... ) = f_{ij} for some constant c

We write u_{ghost} = alpha * u_{ij} + beta, so the equation becomes
    1/h^2 * ( [c-alpha]*u_{ij} + ... ) = f_{ij} + beta/h^2

Thus, the matrix is corrected by adding -alpha/h^2 to the diagonal in row ij, and the RHS vector is corrected by adding +beta/h^2 to row ij.

Note that on a given boundary alpha is constant across that boundary
Note that interior corner points have two ghost points, such that they are doubly corrected 
*/
void PoissonOp::ApplyRobinCorrectionToMatrix( std::shared_ptr<AMP::LinearAlgebra::Matrix> A ) {

    AMP_INSIST( A, "Non-null A required!" );

    // Select boundaries we need to iterate over
    std::vector<int> boundary_ids = { 1, 2 };
    if ( this->getMesh()->getDim() >= 2 ) {
        boundary_ids.push_back( 3 );
        boundary_ids.push_back( 4 );
    }

    auto h = d_db->getScalar<double>("h");

    // Iterate across all boundaries
    for ( auto boundary_id : boundary_ids ) {

        // Get on-process iterator over current boundary
        auto it = d_BoxMesh->getBoundaryIDIterator( d_geomType, boundary_id );
        
        // Get correction on current boundary (it's constant for all DOFs on this boundary)
        double alpha      = RobinCorrectionMatrixCoefficient( boundary_id );
        double correction = -alpha / (h*h);

        // Add correction to all nodes on current boundary; if there are none then "it" is empty
        for ( auto node = it.begin(); node != it.end(); node++ ) {
            std::vector<size_t> dof;
            d_DOFMan->getDOFs( node->globalID(), dof);
            A->addValueByGlobalID( dof[0], dof[0], correction );
        }
    }
    A->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}

/* Applies corrections to vector f that arise from the elimination of the first set of interior vertices. This correction depends on all of the constants a_i, b_i, and r_i  */
void PoissonOp::ApplyRobinCorrectionToVector( std::shared_ptr<AMP::LinearAlgebra::Vector> f ) {

    AMP_INSIST( f, "Non-null f required!" );

    // Select boundaries we need to iterate over
    std::vector<int> boundary_ids = { 1, 2 };
    if ( this->getMesh()->getDim() >= 2 ) {
        boundary_ids.push_back( 3 );
        boundary_ids.push_back( 4 );
    }

    auto h = d_db->getScalar<double>("h");

    // Iterate across all boundaries
    for ( auto boundary_id : boundary_ids ) {

        // Get on-process iterator over current boundary
        auto it = d_BoxMesh->getBoundaryIDIterator( d_geomType, boundary_id );
        
        // Add correction to all nodes on current boundary; if there are none then "it" is empty
        for ( auto node = it.begin(); node != it.end(); node++ ) {
            // Get correction at current point on current boundary
            double beta       = RobinCorrectionVectorCoefficient( boundary_id, *node );
            double correction = beta / (h*h);

            std::vector<size_t> dof;
            d_DOFMan->getDOFs( node->globalID(), dof);
            f->addValuesByGlobalID<double>( 1, &dof[0], &correction );
        }
    }
    f->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}

/* We write u_ghost = alpha_{} * u_interior + beta_{} 
    This function returns the coefficient alpha on the given boundary.
    Note that up to the constants in the BC changing, the formula is the same on all boundaries.
    The formulae are also dimension-agnostic, in that they hold for 1D, 2D, 3D, etc.. */
double PoissonOp::RobinCorrectionMatrixCoefficient( int boundary ) {
    
    // Get constants a, b, c
    double a, b, c;
    if ( boundary == 1 ) {
        a = d_a1, b = d_b1, c = d_c.xx;
    } else if ( boundary == 2 ) {
        a = d_a2, b = d_b2, c = d_c.xx;
    } else if ( boundary == 3 ) {
        a = d_a3, b = d_b3, c = d_c.yy;
    } else if ( boundary == 4 ) {
        a = d_a4, b = d_b4, c = d_c.yy;
    } else {
        AMP_ERROR( "Invalid boundary" );
    }

    double h     = d_db->getScalar<double>( "h" );
    double alpha = (2*c*b - a*h)/(2*c*b + a*h);
    return alpha;
}

/* We write u_ghost = alpha_{} * u_interior + beta_{} 
    This function returns the coefficient beta on the given boundary. 
    Note that up to the constants in the BC changing, the formula is the same on all boundaries.
    The formulae are also dimension-agnostic, in that they hold for 1D, 2D, 3D, etc.. */
double PoissonOp::RobinCorrectionVectorCoefficient( int boundary, AMP::Mesh::MeshElement &node ) {
    
    // Get constants a, b, c
    double a, b, c;
    if ( boundary == 1 ) {
        a = d_a1, b = d_b1, c = d_c.xx;
    } else if ( boundary == 2 ) {
        a = d_a2, b = d_b2, c = d_c.xx;
    } else if ( boundary == 3 ) {
        a = d_a3, b = d_b3, c = d_c.yy;
    } else if ( boundary == 4 ) {
        a = d_a4, b = d_b4, c = d_c.yy;
    } else {
        AMP_ERROR( "Invalid boundary" );
    }

    double h    = d_db->getScalar<double>( "h" );
    double r    = d_robinFunction( boundary, a, b, node );
    double beta = 2*r*h/(2*c*b + a*h);
    return beta;
}

/* Return a constructed CSR matrix corresponding to the discretized Laplacian on the mesh */
std::shared_ptr<AMP::LinearAlgebra::Matrix> PoissonOp::getLaplacianMatrix( ) {

    auto inVec  = this->getRightVector();
    auto outVec = this->getRightVector();

    int meshDim = this->getMesh()->getDim();
    std::map<size_t, colsDataPair> localCSRData;
    if ( meshDim == 1 ) {
        localCSRData = this->getCSRData1D( );
    } else if ( meshDim == 2 ) {
        localCSRData = this->getCSRData2D( );
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


/* Populate vector with function that takes a reference to a MeshElement and returns a double.  */
void PoissonOp::fillWithFunction( std::shared_ptr<AMP::LinearAlgebra::Vector> vec, std::function<double(AMP::Mesh::MeshElement &)> fun ) {

    auto it = d_BoxMesh->getIterator(d_geomType); // Mesh iterator

    // Fill in exact solution vector
    for ( auto elem = it.begin(); elem != it.end(); elem++ ) {
        std::vector<double> u;
        u.push_back( fun( *elem ) );
        std::vector<size_t> i;
        d_DOFMan->getDOFs( elem->globalID(), i );
        vec->setValuesByGlobalID( 1, &i[0], &u[0] );
    }
    vec->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}



/* ------------------------------------------------
    Class implementing heat equation 
------------------------------------------------- */
/* Abstract base class representing a linear heat problem:
        u'(t) - grad \dot ( D * \grad u ) = s(t), u(0) = u_0
    over the spatial domain [0,1]^d, for d = 1 or d = 2.
*/
class HeatModel {

public:

    // The current time of the solution
    double d_currentTime = 0.0;

    bool d_exactSolutionAvailable = false;

    // PoissonOperator on which this is built
    std::shared_ptr<PoissonOp> d_PoissonOp;

    HeatModel( std::shared_ptr<PoissonOp> PoissonOp_) : 
        d_PoissonOp( PoissonOp_ ){ 

        AMP_INSIST( d_PoissonOp,  "Non-null PoissonOp required!" );
    }

    inline double getCurrentTime() { return d_currentTime; };
    inline void setCurrentTime(double currentTime_) { d_currentTime = currentTime_; };

    /* Pure virtual functions */
    virtual double sourceTerm( AMP::Mesh::MeshElement &node ) = 0;
    virtual double initialCondition( AMP::Mesh::MeshElement &node ) = 0;

    /* Virtual functions */
    virtual double exactSolution( AMP::Mesh::MeshElement & ) {
        AMP_ERROR( "Base class cannot provide an implementation of this function" );
    }
};


/* ------------------------------------------------
    Class implementing manufactured heat equation 
------------------------------------------------- */
/* Abstract class representing a linear heat problem:
        u'(t) - grad \dot ( D * \grad u ) = s(t), u(0) = u_0
    over the spatial domain [0,1]^d, for d = 1 or d = 2.

    With a zero source s(t) = 0 and some specific initial condition
*/
class BasicHeatModel : public HeatModel {

public:

    // Call base class' constructor
    BasicHeatModel( std::shared_ptr<PoissonOp> PoissonOp_) : HeatModel( PoissonOp_ ) { }
    
    // Implementation of pure virtual function
    // Dimesionless wrapper around the exact source term functions
    double sourceTerm( AMP::Mesh::MeshElement & ) override {
        return 0.0;
    }

    // Implementation of pure virtual function
    double initialCondition( AMP::Mesh::MeshElement &node ) override {
        return initialCondition_( node );
    }

private:

    // Dimesionless wrapper around the exact solution functions
    double initialCondition_( AMP::Mesh::MeshElement &node ) {
        int meshDim = d_PoissonOp->getMesh()->getDim();
        if ( meshDim == 1 ) {
            double x = ( node.coord() )[0];
            return initialCondition_( x );
        } else if ( meshDim == 2 ) {
            double x = ( node.coord() )[0];
            double y = ( node.coord() )[1];
            return initialCondition_( x, y );
        } else {
            AMP_ERROR( "Invalid dimension" );
        }
    }
    double initialCondition_(double x) {
        double t = 0.0;
        return std::sin((3.0/2.0)*M_PI*x)*std::cos(2*M_PI*t);
    }
    double initialCondition_(double x, double y) {
        double t = 0.0;
        return std::sin((3.0/2.0)*M_PI*x)*std::cos(2*M_PI*t)*std::cos(3*M_PI*y);
    }
};


/* ------------------------------------------------
    Class implementing manufactured heat equation 
------------------------------------------------- */
/* Abstract class representing a linear heat problem:
        u'(t) - grad \dot ( D * \grad u ) = s(t), u(0) = u_0
    over the spatial domain [0,1]^d, for d = 1 or d = 2.

    The diffusion tensor is D = cxx in 1D, and D = diag( cxx, cyy ) in 2D.

    This manufactured problem is constructed using a PoissonOp, which is a discretization of the term - grad \dot ( D * \grad u ).

    1D. Robin boundary conditions are avaliable in the form of:
        a1 * u + b1 * +cxx * du/dx = r1 at x = 0...
        a2 * u + b2 * -cyy * du/dx = r2 at x = 1...

    2D. Robin boundary conditions are avaliable in the form of:
        a1 * u + b1 * -cxx * du/dx = r1 at x = 0...
        a2 * u + b2 * +cxx * du/dx = r2 at x = 1...
        a3 * u + b3 * -cyy * du/dy = r3 at y = 0...
        a4 * u + b4 * +cyy * du/dy = r4 at y = 1...
*/
class ManufacturedHeatModel : public HeatModel {

public:

    // Call base class' constructor
    ManufacturedHeatModel( std::shared_ptr<PoissonOp> PoissonOp_) : HeatModel( PoissonOp_ ) { 
        // Set flag indicating this class does provide an implementation of exactSolution
        d_exactSolutionAvailable = true;
    }
    
    // Implementation of pure virtual function
    // Dimesionless wrapper around the exact source term functions
    double sourceTerm( AMP::Mesh::MeshElement &node ) override {
        int meshDim = d_PoissonOp->getMesh()->getDim();
        if ( meshDim == 1 ) {
            double x = ( node.coord() )[0];
            return sourceTerm_( x );
        } else if ( meshDim == 2 ) {
            double x = ( node.coord() )[0];
            double y = ( node.coord() )[1];
            return sourceTerm_( x, y );
        } else {
            AMP_ERROR( "Invalid dimension" );
        }
    }

    // Implementation of pure virtual function
    double initialCondition( AMP::Mesh::MeshElement &node ) override {
        double currentTime_ = this->getCurrentTime();
        this->setCurrentTime( 0.0 );
        double ic = exactSolution( node );
        this->setCurrentTime( currentTime_ );
        return ic;
    }

    // Dimesionless wrapper around the exact solution functions
    double exactSolution( AMP::Mesh::MeshElement &node ) override {
        int meshDim = d_PoissonOp->getMesh()->getDim();
        if ( meshDim == 1 ) {
            double x = ( node.coord() )[0];
            return exactSolution_( x );
        } else if ( meshDim == 2 ) {
            double x = ( node.coord() )[0];
            double y = ( node.coord() )[1];
            return exactSolution_( x, y );
        } else {
            AMP_ERROR( "Invalid dimension" );
        }
    }

    // Dimesionless wrapper around the Robin functions
    double getRobinValue( int boundary, double a, double b, AMP::Mesh::MeshElement & node ) {
        int meshDim = d_PoissonOp->getMesh()->getDim();
        if ( meshDim == 1 ) {
            return getRobinValue1D( boundary, a, b );
        } else if ( meshDim == 2 ) {
            double x = ( node.coord() )[0];
            double y = ( node.coord() )[1];
            return getRobinValue2D( boundary, a, b, x, y );
        } else {
            AMP_ERROR( "Invalid dimension" );
        }
    }
    
private:
    // Exact solution, its gradient and corresponding source term
    double exactSolution_(double x);
    double exactSolutionGradient(double x);
    double sourceTerm_(double x);
    //
    double exactSolution_(double x, double y);
    double exactSolutionGradient(double x, double y, std::string component);
    double sourceTerm_(double x, double y);

    /* 1D. The boundaries of 0 and 1 are hard coded here. */
    double getRobinValue1D( int boundary, double a, double b ) {
        double cxx = d_PoissonOp->d_c.xx;
        double x = -1.0;
        double normal_sign = 0.0;
        if ( boundary == 1 ) { // West
            x = 0.0;
            // Normal vector is -hat{x}
            normal_sign      = -1.0 * cxx;
        } else if ( boundary == 2 ) { // East
            x = 1.0;
            // Normal vector is +hat{x}
            normal_sign      = +1.0 * cxx;
        }
        double u    = exactSolution_( x );
        double dudx = exactSolutionGradient( x );
        return a * u + normal_sign * b * dudx;
    }

    /* 2D. point is the value of the variable not on the boundary. The boundaries of 0 and 1 are hard coded here, so that the value (either x or y) that is parsed is ignored depending on the value of "boundary" (this just makes calling this function easier) */
    double getRobinValue2D( int boundary, double a, double b, double x, double y ) {
        double cxx = d_PoissonOp->d_c.xx;
        double cyy = d_PoissonOp->d_c.yy;
        std::string normal_direction = "";
        double normal_sign = 0.0;
        if ( boundary == 1 ) { // West
            x = 0.0;
            // Normal vector is -hat{x}
            normal_direction = "x";
            normal_sign      = -1.0 * cxx;
        } else if ( boundary == 2 ) { // East
            x = 1.0;
            // Normal vector is +hat{x}
            normal_direction = "x";
            normal_sign      = +1.0 * cxx;
        } else if ( boundary == 3 ) { // South
            y = 0.0;
            // Normal vector is -hat{y}
            normal_direction = "y";
            normal_sign      = -1.0 * cyy;
        } else if ( boundary == 4 ) { // North
            y = 1.0;
            // Normal vector is +hat{y}
            normal_direction = "y";
            normal_sign      = +1.0 * cyy;
        }
        double u    = exactSolution_( x, y );
        double dudn = exactSolutionGradient( x, y, normal_direction );
        return a * u + normal_sign * b * dudn;
    }
}; 

double ManufacturedHeatModel::exactSolution_(double x) {
    double t = this->getCurrentTime();
    return std::sin((3.0/2.0)*M_PI*x)*std::cos(2*M_PI*t);
}

double ManufacturedHeatModel::exactSolutionGradient(double x) {
    double t = this->getCurrentTime();
    return (3.0/2.0)*M_PI*std::cos(2*M_PI*t)*std::cos((3.0/2.0)*M_PI*x);
}

double ManufacturedHeatModel::sourceTerm_(double x) {
    double t = this->getCurrentTime();
    double cxx = d_PoissonOp->d_c.xx;
    return (1.0/4.0)*M_PI*(9*M_PI*cxx*std::cos(2*M_PI*t) - 8*std::sin(2*M_PI*t))*std::sin((3.0/2.0)*M_PI*x);
}

double ManufacturedHeatModel::exactSolution_(double x, double y) {
    double t = this->getCurrentTime();
    return std::sin((3.0/2.0)*M_PI*x)*std::cos(2*M_PI*t)*std::cos(3*M_PI*y);
}

double ManufacturedHeatModel::exactSolutionGradient(double x, double y, std::string component) {
    double t = this->getCurrentTime();
    double grad = 0.0;
    if ( component == "x" ) {
        grad = (3.0/2.0)*M_PI*std::cos(2*M_PI*t)*std::cos((3.0/2.0)*M_PI*x)*std::cos(3*M_PI*y);
    } else if ( component == "y" ) {
        grad = -3*M_PI*std::sin((3.0/2.0)*M_PI*x)*std::sin(3*M_PI*y)*std::cos(2*M_PI*t);
    } else {
        AMP_ERROR( "Component not recognised" );
    }
    return grad;
}

double ManufacturedHeatModel::sourceTerm_(double x, double y) {
    double t = this->getCurrentTime();
    double cxx = d_PoissonOp->d_c.xx;
    double cyy = d_PoissonOp->d_c.yy;
    return (1.0/4.0)*M_PI*(9*M_PI*cxx*std::cos(2*M_PI*t) + 36*M_PI*cyy*std::cos(2*M_PI*t) - 8*std::sin(2*M_PI*t))*std::sin((3.0/2.0)*M_PI*x)*std::cos(3*M_PI*y);
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
    OpDB->putScalar<int>( "print_info_level", 1 );
    OpDB->putScalar<std::string>( "name", "PoissonOp" );  
    auto OpParameters = std::make_shared<AMP::Operator::OperatorParameters>( OpDB );
    OpParameters->d_name = "PoissonOp";
    OpParameters->d_Mesh = mesh;

    // Create PoissonOp 
    auto myPoissonOp = std::make_shared<PoissonOp>( OpParameters );    
    // Create BEOp based on Poisson Op
    OpDB->putScalar<bool>( "matrixRequired", true );
    auto myBEOp      = std::make_shared<BELinearOp>( OpParameters );  
    myBEOp->setRHSOperator( myPoissonOp );
    

    /****************************************************************
    * Create heat equation model                                    *
    ****************************************************************/
    std::shared_ptr<HeatModel> myHeatModel;
    auto useManufacturedModel = PDE_db->getWithDefault<bool>( "manufacturedModel", false );
    auto useBasicModel = !useManufacturedModel;

    if ( useBasicModel ) {
        auto myHeatModel_ = std::make_shared<BasicHeatModel>( myPoissonOp );
        myHeatModel        = myHeatModel_;

    } else {
        AMP::pout << "Manufactured heat model is being used" << std::endl;
        auto myHeatModel_ = std::make_shared<ManufacturedHeatModel>( myPoissonOp );
        
        // Point the Robin BC values in the PoissonOp to those given by the manufactured problem
        myPoissonOp->setRobinFunction( std::bind( &ManufacturedHeatModel::getRobinValue, &( *myHeatModel_ ), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4 ) );
        myHeatModel = myHeatModel_;        
    }

    // Create hassle-free wrappers around ic, source term and exact solution
    auto icFun        = std::bind( &HeatModel::initialCondition, &( *myHeatModel ), std::placeholders::_1 );
    auto PDESourceFun = std::bind( &HeatModel::sourceTerm, &( *myHeatModel ), std::placeholders::_1 );
    auto uexactFun    = std::bind( &HeatModel::exactSolution, &( *myHeatModel ), std::placeholders::_1 );


    /****************************************************************
    * Set up relevant vectors                                       *
    ****************************************************************/
    // Create required vectors over the mesh
    auto numSolVec    = myPoissonOp->getRightVector();
    auto manSolVec    = myPoissonOp->getRightVector();
    auto errorVec     = myPoissonOp->getRightVector();
    auto BDFSourceVec = myPoissonOp->getRightVector();

    // Create initial condition vector
    auto ic = myPoissonOp->getRightVector();
    myPoissonOp->fillWithFunction( ic, icFun );
    ic->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );

    // Create vectors to hold current and new solution (when integrating)
    auto sol_old = ic->clone( );
    sol_old->copyVector( ic ); 
    auto sol_new = ic->clone( );
    sol_new->copyVector( ic ); 

    /****************************************************************
    * Set up implicit time integrator                               *
    ****************************************************************/
    // Ensure BDF integrator is being used
    auto bdf_ti = { "Backward Euler", "BDF1", "BDF2", "BDF3", "BDF4", "BDF5", "BDF6" };
    auto ti_name   = ti_db->getScalar<std::string>( "name" );
    auto is_bdf = ( std::find( bdf_ti.begin(), bdf_ti.end(), ti_name ) != bdf_ti.end() );
    AMP_INSIST( is_bdf, "Implementation assumes BDF integrator" );

    // Parameters for time integrator
    double dt = 0.5*PDE_db->getScalar<double>( "h" );
    auto tiParams           = std::make_shared<AMP::TimeIntegrator::TimeIntegratorParameters>( ti_db );
    tiParams->d_ic_vector   = ic;
    tiParams->d_operator    = myBEOp;
    tiParams->d_pSourceTerm = BDFSourceVec; // Point source vector to our source vector

    // Create timeIntegrator from factory
    std::shared_ptr<AMP::TimeIntegrator::TimeIntegrator> timeIntegrator = AMP::TimeIntegrator::TimeIntegratorFactory::create( tiParams );
    
    // Cast to implicit integrator
    auto implicitIntegrator =
            std::dynamic_pointer_cast<AMP::TimeIntegrator::ImplicitIntegrator>( timeIntegrator );

    // Tell implicitIntegrator how to tell our operator what the time step is
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
        myHeatModel->setCurrentTime( T + dt ); // Set manufactured solution to new time---this ensures the source term and Robin values are sampled at the new time.
        // Fill vector with PDE source term
        myPoissonOp->fillWithFunction( BDFSourceVec, PDESourceFun );
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
        if ( myHeatModel->d_exactSolutionAvailable ) {
            myHeatModel->setCurrentTime( T );
            myPoissonOp->fillWithFunction( manSolVec, uexactFun );
            errorVec->subtract( *sol_new, *manSolVec );
            AMP::pout << "----------------------------------------" << std::endl;
            AMP::pout << "Manufactured discretization error norms:" << std::endl;
            auto enorms = getDiscreteNorms( PDE_db->getScalar<double>( "h" ), errorVec );
            AMP::pout << "||e||=(" << enorms[0] << "," << enorms[1] << "," << enorms[2] << ")" << std::endl;
            AMP::pout << "----------------------------------------" << std::endl;

            // Write manufactured and numerical solution to file.
            #if 0
            step++;
            std::string name = std::to_string( T );
            manSolVec->setName( name );
            sol_new->setName( name );
            AMP::IO::AsciiWriter vecWriter_man;
            vecWriter_man.registerVector( manSolVec );
            vecWriter_man.writeFile( man_dir, step, T  );
            AMP::IO::AsciiWriter vecWriter_num;
            vecWriter_num.registerVector( sol_new );
            vecWriter_num.writeFile( num_dir, step, T  );
            #endif
        }
        

        // Drop out if we've exceeded max steps
        if ( !timeIntegrator->stepsRemaining() ) {
            AMP_WARNING( "max_integrator_steps has been reached, dropping out of loop now..." );
            break;
        }
    }

    
    // No specific solution is implemented for this problem, so this will just check that the solver converged. 
    checkConvergence( implicitIntegrator->getSolver().get(), input_db, input_file, *ut );

}
// end of driver()



/* Numerical implementation of the heat equation with finite difference spatial discretization and Robin boundary conditions.

    Some of the code is generated in the associated -helper.py file.

    In 1D, the PDE is u_t - cxx*u_xx = f
    In 2D, the PDE is u_t - (cxx*u_xx - cyy*u_yy) = f

    The input file must contain the following databases:
    1. PDE
    2. TimeIntergrator
        2.1 LinearSolver

    The PDE database must contains the following keys:
        dim - PDE dimension, 1 or 2.
        n   - number of grid points in each spatial dimension
        aj, bj, rj for j = 1,2 in 1D, and j = 1,2,3,4 in 2D. On the jth spatial boundary a Robin condition is imposed of the form aj*u + bj*du/dn = rj, with constants specified in the input file.

    Optional keys in the PDE database are:
        manufacturedModel - Boolean specifying whether to use a manufactured source and Robin BC values with corresponding exact solution implemented. If this is boolean is true then the Robin constants rj are ignored, and the numerical solution is compared against the manufactured solution at each time point. If this key is flase/not provided then a zero source term is used with some initial condition.


    The TimeIntegrator must be of BDF type.
*/
int main( int argc, char **argv )
{

    AMP::AMPManager::startup( argc, argv );      

    AMP::UnitTest ut;

    // Create a global communicator
    AMP::AMP_MPI comm( AMP_COMM_WORLD );
    
    std::vector<std::string> exeNames;
    exeNames.emplace_back( "testLinearSolvers-HeatRobinFD-CG" );

    for ( auto &exeName : exeNames ) {
        driver( comm, &ut, exeName );
    }    
    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}


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
        mesh_db->putVector<int>( "Size", { n-3 } ); // mesh has n-2 points
        mesh_db->putVector<double>( "Range", { h + 0.5*h, 1.0 - (h + 0.5*h) } );
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
    Class implementing a discrete Poisson problem 
------------------------------------------------- */
class PoissonOp : public AMP::Operator::LinearOperator {

private:
    
    std::shared_ptr<AMP::Mesh::BoxMesh> d_BoxMesh;
    PDECoefficients                     d_c;
    AMP::Mesh::GeomType                 d_geomType = AMP::Mesh::GeomType::Vertex;

    void setPDECoefficients();
    std::shared_ptr<AMP::LinearAlgebra::Matrix> getLaplacianMatrix();
    std::shared_ptr<AMP::LinearAlgebra::Matrix> d_APoisson;
    
    double d  = 1.0; // diffusion coefficient... make consistent with c.xx
    
    double a0 = 0.25, b0 = -1.0/6.0;
    double a1 = 0.25, b1 = +1.0/6.0;

    // Dirichlet BCs
    // double a0 = 1.0, b0 = 0.0;
    // double a1 = 1.0, b1 = 0.0;

    mutable double r0 = std::nan(""), r1 = std::nan("");

    void ApplyRobinCorrectionToMatrix( std::shared_ptr<AMP::LinearAlgebra::Matrix> A );


public:

    
    void ApplyRobinCorrectionToVector( std::shared_ptr<AMP::LinearAlgebra::Vector> f );

    mutable double currentTime = 0.0;
    std::shared_ptr<AMP::Database>                   d_db;
    std::shared_ptr<AMP::Discretization::DOFManager> d_DOFMan;
    double d_gamma;
    
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

        // Set PDE coefficients
        this->setPDECoefficients();

        // Get the matrix
        auto A = getLaplacianMatrix( );

        // Update according to Robin data
        ApplyRobinCorrectionToMatrix( A );
        
        // Set linear operator's matrix
        //this->setMatrix( A );
        d_APoisson = A;
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
    void fillWith_uexact( std::shared_ptr<AMP::LinearAlgebra::Vector> uexact  );
    void fillWith_fsource( std::shared_ptr<AMP::LinearAlgebra::Vector> fsource );
   

private:


    /* Build and set d_DOFMan */
    void set_DOFManager() {
        int DOFsPerElement = 1; 
        int gcw  = 1; // Ghost-cell width (stencils are at most 3-point in each direction)
        d_DOFMan = AMP::Discretization::boxMeshDOFManager::create(this->getMesh(), d_geomType, gcw, DOFsPerElement);
    }

    // Coefficients defining the PDE
    std::vector<double> getPDECoefficients1D();

    // Exact solution and corresponding source term
    double exactSolution(double x);
    double exactSolutionGradient(double x);
    double sourceTerm(double x);
    
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

public:
    inline double getCurrentTime() { return currentTime; };
    inline void setCurrentTime(double currentTime_) { currentTime = currentTime_; 
    
        // Set Robin BC values for manufactured solution 
        // Note that the Robin BCs evaluate the solution exactly at x = 0 and x = 1, even though these are not computational vertices (they are half way between a ghost vertex and the first interior vertex, which ultimately is eliminated from the system)
        double E0    = exactSolution( 0.0 );
        double E1    = exactSolution( 1.0 );
        double dEdx0 = exactSolutionGradient( 0.0 );
        double dEdx1 = exactSolutionGradient( 1.0 );
        r0 = a0 * E0 + b0 * dEdx0;
        r1 = a1 * E1 + b1 * dEdx1;
    };

    // Set the time-step size in the operator
    void setGamma( AMP::Scalar gamma )
    {
        d_gamma = double( gamma );
    }

    // only required if we are doing multi-physics and scaling of components needed
    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u_in,
                AMP::LinearAlgebra::Vector::shared_ptr r ) override
    {
        AMP::LinearAlgebra::Vector::const_shared_ptr u;

        if ( d_pSolutionScaling ) {

            AMP_ASSERT( d_pFunctionScaling );

            if ( !d_pScratchSolVector ) {
                d_pScratchSolVector = u_in->clone();
            }

            d_pScratchSolVector->multiply( *u_in, *d_pSolutionScaling );
            d_pScratchSolVector->makeConsistent();
            u = d_pScratchSolVector;

        } else {
            u = u_in;
        }

        // Compute r <- (I + gamma*A)*u
        d_APoisson->mult( u, r ); // r <- A*u
        r->axpby (1.0, d_gamma, *u); // r <- 1.0*u + d_gamma * r
 	
        if ( d_pFunctionScaling ) {
            r->divide( *r, *d_pFunctionScaling );
        }
    }

    // only required if we are doing multi-physics and scaling of components needed
    void setComponentScalings( std::shared_ptr<AMP::LinearAlgebra::Vector> s,
                               std::shared_ptr<AMP::LinearAlgebra::Vector> f )
    {
        d_pSolutionScaling = s;
        d_pFunctionScaling = f;
    }

    // only required if we are doing multi-physics and scaling of components needed
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_pSolutionScaling;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_pFunctionScaling;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_pScratchSolVector;

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


double PoissonOp::exactSolution(double x) {
    double t = this->getCurrentTime();
    return std::sin((3.0/2.0)*M_PI*x)*std::cos(2*M_PI*t);
}

double PoissonOp::exactSolutionGradient(double x) {
    double t = this->getCurrentTime();
    return (3.0/2.0)*M_PI*std::cos(2*M_PI*t)*std::cos((3.0/2.0)*M_PI*x);
}

double PoissonOp::sourceTerm(double x) {
    double t = this->getCurrentTime();
    double cxx = d_c.xx;
    return (1.0/4.0)*M_PI*(9*M_PI*cxx*std::cos(2*M_PI*t) - 8*std::sin(2*M_PI*t))*std::sin((3.0/2.0)*M_PI*x);
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



/* Populate exact solution vector.  */
void PoissonOp::fillWith_uexact( std::shared_ptr<AMP::LinearAlgebra::Vector> uexact ) {

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
void PoissonOp::fillWith_fsource( std::shared_ptr<AMP::LinearAlgebra::Vector> fsource ) {

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

/* 
Robin set up. Use cell centered discretization. One ghost cell is used. The boundary condition explicitly uses the ghost cell and the first interior cell. The boundary condition is used as the equation for the first interior cell, and this equation is eliminated.

There are n-1 active computational vertices, 0, 1, ..., n-2. There are two ghost vertices, and the first interior vertices are eliminated.

at x = 0, the BC is
a0 * E - b0 * d_E * dE/dx = r0

We solve this as:
    E0 = alpha1*E1 + alpha2*E2 + beta0
    1/h^2*( - E0 + 2*E1 - E2 ) = 1/h^2*( [2-alpha1]*E1 - [1+alpha2]*E2 ) - 1/h^2 * beta0

at x = 1, the BC is
a1 * E - b1 * d_E * dE/dx = r1.

Suppose the last interior point is E3, then, we solve this as:
    beta1 + alpha1*E1 + alpha2*E2 = E3
The stencil at the last active interior vertex is:
    1/h^2*( - E1 + 2*E2 - E3 ) = 1/h^2*( - [1+alpha1]*E1 + [2-alpha2]*E2 ) - 1/h^2 * beta1

*/
void PoissonOp::ApplyRobinCorrectionToMatrix( std::shared_ptr<AMP::LinearAlgebra::Matrix> A ) {

    auto localBox  = getLocalNodeBox( d_BoxMesh );
    auto globalBox = getGlobalNodeBox( d_BoxMesh );
    auto h = d_db->getScalar<double>("h");

    // West boundary
    if ( localBox.first[0] == globalBox.first[0] ) {

        int i = localBox.first[0];
        size_t dof = gridIndsToDOF( i );

        double alpha1 = (3*a0*h - 6*b0*d)/(4*a0*h - 4*b0*d);
        double alpha2 = (-a0*h + 2*b0*d)/(4*a0*h - 4*b0*d);

        //     We solve this as:
        // E0 = alpha1*E1 + alpha2*E2 + beta0
        // 1/h^2*( - E0 + 2*E1 - E2 ) = 1/h^2*( [2-alpha1]*E1 - [1+alpha2]*E2 ) - 1/h^2 * beta0
        std::vector<size_t> cols = { dof, gridIndsToDOF( i+1 ) };
        std::vector<double> vals = { 1.0/(h*h) * ( 2.0 - alpha1 ), -1.0/(h*h) * ( 1.0 + alpha2 ) };

        // Set values by grid indices
        int nrows = 1;
        int ncols = 2;
        A->setValuesByGlobalID<double>( nrows, ncols, &dof, cols.data(), vals.data() );
        A->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
    }

    // East boundary
    if ( localBox.last[0] == globalBox.last[0] ) {
        int i = localBox.last[0];
        size_t dof  = gridIndsToDOF( i );

        double alpha1 = (-a1*h - 2*b1*d)/(4*a1*h + 4*b1*d);
        double alpha2 = (3*a1*h + 6*b1*d)/(4*a1*h + 4*b1*d);

        // beta1 + alpha1*E1 + alpha2*E2 = E3
        // The stencil at the last active interior vertex is:
        //     1/h^2*( - E1 + 2*E2 - E3 ) = 1/h^2*( - [1+alpha1]*E1 + [2-alpha2]*E2 ) - 1/h^2 * beta1
        std::vector<size_t> cols = { gridIndsToDOF( i-1 ), dof };
        std::vector<double> vals = { -1.0/(h*h) * ( 1.0 + alpha1 ), 1.0/(h*h) * ( 2.0 - alpha2 ) };

        // Set values by grid indices
        int nrows = 1;
        int ncols = 2;
        
        A->setValuesByGlobalID<double>( nrows, ncols, &dof, cols.data(), vals.data() );
        A->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
    }
}

/* Correctly set values in boundary points of vector */
void PoissonOp::ApplyRobinCorrectionToVector( std::shared_ptr<AMP::LinearAlgebra::Vector> f ) {

    auto localBox  = getLocalNodeBox( d_BoxMesh );
    auto globalBox = getGlobalNodeBox( d_BoxMesh );
    auto h = d_db->getScalar<double>("h");

    
    // West boundary
    if ( localBox.first[0] == globalBox.first[0] ) {

        int i     = localBox.first[0];
        int dof   = gridIndsToDOF( i );
        AMP::Mesh::BoxMesh::MeshElementIndex ind(
                        AMP::Mesh::GeomType::Vertex, 0, i);
        auto elem = d_BoxMesh->getElement( ind );

        double x0 = ( elem.coord() )[0];
        double f0 = sourceTerm( x0 );

        double beta0  = h*r0/(2*a0*h - 2*b0*d);
        double val    = 1.0/(h*h) * beta0 + f0;
        
        // Set values
        f->setValueByGlobalID<double>( dof, val );
        f->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
    }

    // East boundary
    if ( localBox.last[0] == globalBox.last[0] ) {
        
        int i     = localBox.last[0];
        int dof   = gridIndsToDOF( i );
        AMP::Mesh::BoxMesh::MeshElementIndex ind(
                        AMP::Mesh::GeomType::Vertex, 0, i);
        auto elem = d_BoxMesh->getElement( ind );

        double x1 = ( elem.coord() )[0];
        double f1 = sourceTerm( x1 );

        double beta1  = h*r1/(2*a1*h + 2*b1*d);
        double val    = 1.0/(h*h) * beta1 + f1;
        
        // Set values by grid indices
        f->setValueByGlobalID<double>( dof, val );
        f->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
    }
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


/* Linear system we'll have to solve is [I - gamma*dt]*u = f. For gamma*dt < 1 the matrix is SPD. */
void updateDatabaseIfImplicit( std::shared_ptr<AMP::Database> db, std::string implicitSolverName )
{
    AMP_ASSERT( db );
    auto imp_ti      = { "Backward Euler", "BDF1", "BDF2", "BDF3", "BDF4", "BDF5", "BDF6" };
    auto name        = db->getScalar<std::string>( "name" );
    auto is_implicit = ( std::find( imp_ti.begin(), imp_ti.end(), name ) != imp_ti.end() );
    db->putScalar<bool>( "is_implicit", is_implicit );
    if ( is_implicit ) {
        db->putScalar<bool>( "user_managed_time_operator", true );
        db->putScalar<std::string>( "implicit_integrator", name );
        db->putScalar<std::string>( "solver_name", "Solver" );
        db->putScalar<std::string>( "timestep_selection_strategy", "constant" );
        db->putScalar<bool>( "use_predictor", false );
        // add the next line to turn off component scaling for multi-physics
        //        db->putScalar<bool>( "auto_component_scaling", false );
        int print_info_level = 1;
        // if (implicitSolverName == "CGSolver") 
        //     print_info_level = 0;
        auto solver_db = AMP::Database::create( "name",
                                                implicitSolverName,
                                                "print_info_level",
                                                print_info_level,
                                                "max_iterations",
                                                256,
                                                "absolute_tolerance",
                                                1.0e-8,
                                                "relative_tolerance",
                                                1.0e-8 );
        db->putDatabase( "Solver", std::move( solver_db ) );
    }
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

    AMP_INSIST( PDE_db,    "A PDE database must be provided" );
    //AMP_INSIST( solver_db, "A LinearSolver database must be provided" );

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
    auto umanVec    = AOp->getRightVector();
    auto fsourceVec = AOp->getRightVector();



    // // Create a MyOperator
    // const auto opDB = std::make_shared<AMP::Database>( "OperatorDB" );
    // opDB->putScalar<std::string>( "name", "MyOperator" );  
    // auto opParams = std::make_shared<AMP::Operator::OperatorParameters>( opDB );
    // auto myOp = std::make_shared<MyOperator>( opParams );
    // myOp->setDOFManager( comm );

    // Parameters for time integraor
    double finalTime       = 0.5;
    double dt              = 0.001;
    int maxIntegratorSteps = 300;
    //auto integratorName = "BDF2";
    auto integratorName = "BDF3";
    //auto integratorName = "ExplicitEuler";
    //auto implicitSolverName ="CGSolver";
    auto implicitSolverName ="GMRESSolver";

    // Create initial condition vector
    auto ic = AOp->getRightVector();
    AOp->setCurrentTime( 0.0 );
    AOp->fillWith_uexact( ic );
    ic->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );

    // Create vectors to hold current and new solutions
    auto sol_old = ic->clone( );
    sol_old->copyVector( ic ); 
    //sol_old->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
    auto sol_new = ic->clone( );
    sol_new->copyVector( ic ); 
    //sol_new->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );

    // Create Database to construct TimeIntegratorParameters from
    std::shared_ptr<AMP::Database> db = AMP::Database::create( "name",
                                                               integratorName,
                                                               "initial_time",
                                                               0.0,
                                                               "final_time",
                                                               finalTime,
                                                               "max_integrator_steps",
                                                               maxIntegratorSteps,
                                                               "print_info_level",
                                                               2 );
    auto tiParams           = std::make_shared<AMP::TimeIntegrator::TimeIntegratorParameters>( db );
    tiParams->d_ic_vector   = ic;
    tiParams->d_operator    = AOp;
    // Create a source vector
    tiParams->d_pSourceTerm = AOp->getRightVector();
    
    updateDatabaseIfImplicit( db, implicitSolverName );

    // Create timeIntegrator from factory
    std::shared_ptr<AMP::TimeIntegrator::TimeIntegrator> timeIntegrator = AMP::TimeIntegrator::TimeIntegratorFactory::create( tiParams );
    
    auto implicitIntegrator =
            std::dynamic_pointer_cast<AMP::TimeIntegrator::ImplicitIntegrator>( timeIntegrator );

    // required to set gamma for the user operator at each timestep
    implicitIntegrator->setTimeScalingFunction(
        std::bind( &PoissonOp::setGamma, &( *AOp ), std::placeholders::_1 ) );

    // required only if a user operator is being used which is multi-physics
    implicitIntegrator->setComponentScalingFunction(
        std::bind( &PoissonOp::setComponentScalings,
                    &( *AOp ),
                    std::placeholders::_1,
                    std::placeholders::_2 ) );

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
        umanVec->setName( name );
        sol_new->setName( name );
        vecWriter_man.registerVector( ic );
        vecWriter_man.writeFile( man_dir, step, T  );
        AMP::IO::AsciiWriter vecWriter_num;
        vecWriter_num.registerVector( sol_new );
        vecWriter_num.writeFile( num_dir, step, T  );
    }


    // Integrate!
    double T = 0.0;
    timeIntegrator->setInitialDt( dt );
    while ( T < finalTime ) {

        // Set the solution-independent source term; note that this approach can only work for multistep methods, since then the new solution 
        if ( db->getScalar<bool>( "is_implicit" ) ) {
            AOp->setCurrentTime( T + dt ); // Set source to new time...
        } else {
            AOp->setCurrentTime( T ); // Set source to current time...
        }
        AOp->fillWith_fsource( tiParams->d_pSourceTerm );
        AOp->ApplyRobinCorrectionToVector( tiParams->d_pSourceTerm );

        // Advance the solution
        timeIntegrator->advanceSolution( dt, T == 0, sol_old, sol_new );
        if ( timeIntegrator->checkNewSolution() ) {
            timeIntegrator->updateSolution();
            sol_new->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
            sol_old->copyVector( sol_new );
            
        } else {
            AMP_ERROR( "Solution didn't converge" );
        }

        // Update time
        T += dt;

        // Print exact vs approximate solution 
        AOp->setCurrentTime( T );
        AOp->fillWith_uexact( umanVec );
        
        step++;
        std::string name = std::to_string( T );
        umanVec->setName( name );
        sol_new->setName( name );
        AMP::IO::AsciiWriter vecWriter_man;
        vecWriter_man.registerVector( umanVec );
        vecWriter_man.writeFile( man_dir, step, T+dt  );
        AMP::IO::AsciiWriter vecWriter_num;
        vecWriter_num.registerVector( sol_new );
        vecWriter_num.writeFile( num_dir, step, T+dt  );
        
        
        umanVec->subtract( *sol_new, *umanVec );

        AMP::pout << "Discretization error norms:" << std::endl;
        auto enorms = getDiscreteNorms( PDE_db->getScalar<double>( "h" ), umanVec );
        AMP::pout << "||e||=(" << enorms[0] << "," << enorms[1] << "," << enorms[2] << ")" << std::endl;
        AMP::pout << "----------------------" << std::endl;

        //auto e = umanVec->maxNorm();
        //AMP::pout << "t=" << T << ": error=" << e << std::endl;

        

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
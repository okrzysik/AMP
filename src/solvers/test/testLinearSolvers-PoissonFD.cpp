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

#include <iostream>


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

// CSR structure for identity row scaled by constant c
colsDataPair localCSR_scaledIdentity(size_t dof, double c) { 
    std::vector<double> vals = { c };  // data
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


/* Create a d-dimensional BoxMesh over [0,1]^d with n+1 points in each direction. Also store the mesh spacing h in the provided database. */
std::shared_ptr<AMP::Mesh::BoxMesh> createBoxMesh( AMP::AMP_MPI comm, std::shared_ptr<AMP::Database> PDE_db )
{
    auto n   = PDE_db->getScalar<int>( "n" );
    auto dim = PDE_db->getScalar<int>( "dim" );
    // Store the mesh size: There are n+1 points on the mesh, including the boundaries
    double h = 1.0 / n;
    PDE_db->putScalar( "h",   h );

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
    
    std::shared_ptr<AMP::Mesh::BoxMesh> d_BoxMesh;
    AMP::Mesh::GeomType                 d_geomType = AMP::Mesh::GeomType::Vertex;

    void setPDECoefficients();
    std::shared_ptr<AMP::LinearAlgebra::Matrix> getLaplacianMatrix();
    
public:

    PDECoefficients                                  d_c;
    std::shared_ptr<AMP::Database>                   d_db;
    std::shared_ptr<AMP::Discretization::DOFManager> d_DOFMan;
    

    // Constructor call's base class's constructor
    PoissonOp(std::shared_ptr<const AMP::Operator::OperatorParameters> params_) : 
            AMP::Operator::LinearOperator( params_ ) { 

        // Keep a pointer to my BoxMesh to save having to do this downcast repeatedly 
        d_BoxMesh = std::dynamic_pointer_cast<AMP::Mesh::BoxMesh>(this->getMesh());
        AMP_INSIST( d_BoxMesh, "Mesh must be a AMP::Mesh::BoxMesh" );

        // Set my database
        d_db = params_->d_db->getDatabase( "PDE_db" );

        // Unpack constants from the database
        d_u1 = d_db->getScalar<double>( "u1" );
        d_u2 = d_db->getScalar<double>( "u2" );
        //
        if ( d_BoxMesh->getDim() >= 2 ) {
            d_u3 = d_db->getScalar<double>( "u3" );
            d_u4 = d_db->getScalar<double>( "u4" );
        }
        if ( d_BoxMesh->getDim() >= 3 ) {
            d_u5 = d_db->getScalar<double>( "u5" );
            d_u6 = d_db->getScalar<double>( "u6" );
        }

        // Specify Dirichlet return function as the default
        std::function<double( int, AMP::Mesh::MeshElement & )> wrapper = [&]( int boundary, AMP::Mesh::MeshElement & ) { return dirichletFunctionDefault( boundary ); };
        this->setDirichletFunction( wrapper );

        // Set DOFManager
        this->set_DOFManager();
        AMP_INSIST(  d_DOFMan, "Requires non-null DOF" );

        // Set PDE coefficients
        this->setPDECoefficients();

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

    
    // Populate a vector with the given function
    void fillWithFunction( std::shared_ptr<AMP::LinearAlgebra::Vector> u, std::function<double(AMP::Mesh::MeshElement &)> fun );

    // The user can specify any Dirichlet return function with this signature; if they do then this will overwrite the default.
    void setDirichletFunction( std::function<double(int, AMP::Mesh::MeshElement &)> fn_ ) { d_dirichletFunction = fn_; };

    void applyDirichletCorrectionToVector( std::shared_ptr<AMP::LinearAlgebra::Vector> f );

private:

    // Constants in Dirichlet BCs
    double d_u1, d_u2, d_u3, d_u4, d_u5, d_u6; 

      // Prototype of function returning value of Dirichlet BC on given boundary at given node. The user can specify any function with this signature
    std::function<double( int boundary, AMP::Mesh::MeshElement & node )> d_dirichletFunction;

    /* Build and set d_DOFMan */
    void set_DOFManager() {
        int DOFsPerElement = 1; 
        int gcw  = 1; // Ghost-cell width (stencils are at most 3-point in each direction)
        d_DOFMan = AMP::Discretization::boxMeshDOFManager::create(this->getMesh(), d_geomType, gcw, DOFsPerElement);
    }

    // Coefficients defining the PDE
    std::vector<double> getPDECoefficients1D();
    std::vector<double> getPDECoefficients2D();
    std::vector<double> getPDECoefficients3D();
    
    // FD stencils
    std::vector<double> getStencil1D();
    std::vector<double> getStencil2D();
    std::vector<double> getStencil3D();

    // Map from row index to col+data in that row
    std::map<size_t, colsDataPair> getCSRData1D();
    std::map<size_t, colsDataPair> getCSRData2D();
    std::map<size_t, colsDataPair> getCSRData3D();

    // Map from grid index i, or i,j, or i,j,k to a MeshElementIndex to a MeshElementId and then to the corresponding DOF
    size_t gridIndsToDOF( int i, int j = 0, int k = 0 ) {
        AMP::Mesh::BoxMesh::MeshElementIndex ind(
                        AMP::Mesh::GeomType::Vertex, 0, i, j, k );
        AMP::Mesh::MeshElementID id = d_BoxMesh->convert( ind );
        std::vector<size_t> dof;
        d_DOFMan->getDOFs(id, dof);
        return dof[0];
    };

    // Default function for returning Dirichlet values; this can be overridden by the user
    double dirichletFunctionDefault( int boundary ) {
        if ( boundary == 1 ) {
            return d_u1;
        } else if ( boundary == 2 ) {
            return d_u2;
        } else if ( boundary == 3 ) {
            return d_u3;
        } else if ( boundary == 4 ) {
            return d_u4;
        } else if ( boundary == 5 ) {
            return d_u5;
        } else if ( boundary == 6 ) {
            return d_u6;
        } else { 
            AMP_ERROR( "Invalid boundary" );
        }
    }
}; 

/* Applies corrections to vector f for non-eliminated boundary points; this amounts to putting the given Dirichlet boundary values in the corresponding rows of f (and scaling them by whatever is in the diagonal entry of A).
Note that corner points are set multiple times since they're owned by multiple boundaries.
*/
void PoissonOp::applyDirichletCorrectionToVector( std::shared_ptr<AMP::LinearAlgebra::Vector> f ) {

    AMP_INSIST( f, "Non-null f required!" );

    // Select boundaries we need to iterate over
    std::vector<int> boundary_ids = { 1, 2 };
    if ( this->getMesh()->getDim() >= 2 ) {
        boundary_ids.push_back( 3 );
        boundary_ids.push_back( 4 );
    }
    if ( this->getMesh()->getDim() >= 3 ) {
        boundary_ids.push_back( 5 );
        boundary_ids.push_back( 6 );
    }

    // Placeholder for function value on boundary
    double boundaryValue;

    // Iterate across all boundaries
    for ( auto boundary_id : boundary_ids ) {

        // Get on-process iterator over current boundary
        auto it = d_BoxMesh->getBoundaryIDIterator( d_geomType, boundary_id );
        
        // Add correction to all nodes on current boundary; if there are none then "it" is empty
        for ( auto node = it.begin(); node != it.end(); node++ ) {
            // Get DOF
            std::vector<size_t> dof;
            d_DOFMan->getDOFs( node->globalID(), dof);
            
            // Get function at current point on current boundary
            boundaryValue = d_dirichletFunction( boundary_id, *node );
            // Scale it by diagonal entry of A
            boundaryValue *= this->getMatrix()->getValueByGlobalID( dof[0], dof[0] );

            f->setValueByGlobalID<double>( dof[0], boundaryValue );
        }
    }
}


/* Populate vector with function that takes a reference to a MeshElement and returns a double.  */
void PoissonOp::fillWithFunction( std::shared_ptr<AMP::LinearAlgebra::Vector> vec, std::function<double(AMP::Mesh::MeshElement &)> fun ) {

    double u; // Placeholder for funcation evaluation
    
    // Fill in exact solution vector
    auto it = d_BoxMesh->getIterator(d_geomType); // Mesh iterator
    for ( auto elem = it.begin(); elem != it.end(); elem++ ) {
        u = fun( *elem );
        std::vector<size_t> i;
        d_DOFMan->getDOFs( elem->globalID(), i );
        vec->setValueByGlobalID( i[0], u );
    }
}


void PoissonOp::setPDECoefficients() {
    int dim = this->getMesh()->getDim();
    if ( dim == 1 ) {
        d_c = PDECoefficients( 1, getPDECoefficients1D() );
    } else if ( dim == 2 ) {
        d_c = PDECoefficients( 2, getPDECoefficients2D() );
    } else if ( dim == 3 ) {
        d_c = PDECoefficients( 3, getPDECoefficients3D() );
    }
}


std::vector<double> PoissonOp::getPDECoefficients1D() {
    // PDE coefficients
    double cxx = 1.0;
    std::vector<double> c = { cxx };
    return c;
}

std::vector<double> PoissonOp::getPDECoefficients2D() {
    auto eps    = d_db->getScalar<double>("eps");
    auto theta  = d_db->getScalar<double>("theta");
    
    // Trig functions of angle
    double cth = cos(theta);
    double sth = sin(theta);
    // Diffusion tensor coefficients
    double d11 = std::pow(cth, 2) + eps*std::pow(sth, 2);
    double d22 = std::pow(cth, 2)*eps + std::pow(sth, 2);
    double d12 = cth*sth*(1 - eps);
    // PDE coefficients
    double cxx = d11;
    double cyy = d22;
    double cxy = 2.0*d12;
    std::vector<double> c = { cxx, cyy, cxy };
    return c;
}


std::vector<double> PoissonOp::getPDECoefficients3D() {
    auto epsy  = d_db->getScalar<double>("epsy");
    auto epsz  = d_db->getScalar<double>("epsz");
    auto alpha = d_db->getScalar<double>("alpha");
    auto beta  = d_db->getScalar<double>("beta");
    auto gamma = d_db->getScalar<double>("gamma");

    // Trig functions of angles
    double ca = cos(alpha);
    double sa = sin(alpha);
    double cb = cos(beta);
    double sb = sin(beta);
    double cg = cos(gamma);
    double sg = sin(gamma);
    // Diffusion tensor coefficients
    double d11 = epsy*std::pow(ca*sg + cb*cg*sa, 2) + epsz*std::pow(sa, 2)*std::pow(sb, 2) + std::pow(ca*cg - cb*sa*sg, 2);
    double d22 = std::pow(ca, 2)*epsz*std::pow(sb, 2) + epsy*std::pow(ca*cb*cg - sa*sg, 2) + std::pow(ca*cb*sg + cg*sa, 2);
    double d33 = std::pow(cb, 2)*epsz + std::pow(cg, 2)*epsy*std::pow(sb, 2) + std::pow(sb, 2)*std::pow(sg, 2);
    double d12 = -ca*epsz*sa*std::pow(sb, 2) - epsy*(ca*sg + cb*cg*sa)*(ca*cb*cg - sa*sg) + (ca*cg - cb*sa*sg)*(ca*cb*sg + cg*sa);
    double d13 = sb*(cb*epsz*sa - cg*epsy*(ca*sg + cb*cg*sa) + sg*(ca*cg - cb*sa*sg));
    double d23 = sb*(-ca*cb*epsz + cg*epsy*(ca*cb*cg - sa*sg) + sg*(ca*cb*sg + cg*sa));
    // PDE coefficients
    double cxx = d11;
    double cyy = d22;
    double czz = d33;
    double cxy = 2.0*d12;
    double cxz = 2.0*d13;
    double cyz = 2.0*d23;
    std::vector<double> c = { cxx, cyy, czz, cxy, cxz, cyz };
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
    -cxy*u_xy 

Standard 3-point differences are used for the _ii terms, and 4-point central differences are used for the mixed _ij term.
*/
std::vector<double> PoissonOp::getStencil2D() {
    
    // Unpack PDE coefficients
    double cxx = d_c.xx;
    double cyy = d_c.yy;
    double cxy = d_c.xy;

    // The non-mixed terms
    double O = 0.0;
    // -cxx*u_xx
    double W  = -1.0*cxx; 
           O += +2.0*cxx;  
    double E  = -1.0*cxx;
    // -cyy =*u_yy
    double S  = -1.0*cyy;
           O += +2.0*cyy;
    double N  = -1.0*cyy;

    // Mixed term
    // -cxy*u_xy
    double SW = -0.25*cxy;
    double SE = +0.25*cxy;
    double NW = +0.25*cxy;
    double NE = -0.25*cxy;

    // Populate stencil
    std::vector<double> stencil = { SW, S, SE, W, O, E, NW, N, NE };

    // Introduce 1/h^2 scaling 
    auto h = d_db->getScalar<double>( "h" );
    for ( auto &s : stencil ) {
        s *= 1.0/(h*h);
    }

    return stencil;
}



/* Get 19-point stencil for 3D Poisson that discretizes the operator
    -cxx*u_xx
    -cyy*u_yy
    -czz*u_zz
    -cxy*u_xy
    -cxz*u_xz
    -czz*u_yz

Standard 3-point differences are used for the _ii terms, and 4-point central differences are used for the mixed _ij terms.

The stencil is 19 point rather than 27 point because none of the corers of the 3D cube [-1,0,+1]^3 are included in it. That is, the discretization is defined by a combination of differences that live in the xy, xz, and yz planes. 
*/
std::vector<double> PoissonOp::getStencil3D() {
    
    // Unpack PDE coefficients
    double cxx = d_c.xx;
    double cyy = d_c.yy;
    double czz = d_c.zz;
    double cxy = d_c.xy;
    double cxz = d_c.xz;
    double cyz = d_c.yz; 

    /*
    z = -1 plane (D == down)
    DNE DN DNE
    DW  DO DE
    DSW DS DSE

    z = 0 plane (M == middle)
    MNE MN MNE
    MW  MO ME
    MSW MS MSE

    z = +1 plane (U == up)
    UNE UN UNE
    UW  UO UE
    USW US USE
    */

    // The non-mixed terms
    double MO = 0.0;
    // -cxx*u_xx
    double MW  = -1.0*cxx; 
           MO += +2.0*cxx;  
    double ME  = -1.0*cxx;
    // -czz =*u_yy
    double MS  = -1.0*cyy;
           MO += +2.0*cyy;
    double MN  = -1.0*cyy;
    // -czz*u_zz
    double DO  = -1.0*czz;
           MO += +2.0*czz;
    double UO  = -1.0*czz;

    // Mixed terms
    // -cxy*u_xy
    double MSW = -0.25*cxy;
    double MSE = +0.25*cxy;
    double MNW = +0.25*cxy;
    double MNE = -0.25*cxy;
    // -cxz*u_xz
    double DW  = -0.25*cxz;
    double DE  = +0.25*cxz;
    double UW  = +0.25*cxz;
    double UE  = -0.25*cxz;
    // -cyz*u_yz
    double DS  = -0.25*cyz;
    double DN  = +0.25*cyz;
    double US  = +0.25*cyz;
    double UN  = -0.25*cyz;
    
    // Populate stencil
    std::vector<double> stencil = { 
        DS, DW, DO, DE, DN, // z = -1 plane
        MSW, MS, MSE, MW, MO, ME, MNW, MN, MNE, // z = 0 plane
        US, UW, UO, UE, UN // z = +1 plane
    };

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

    // Iterate over local box
    for ( auto i = localBox.first[0]; i <= localBox.last[0]; i++ ) {
        
        // The current row
        size_t dof = gridIndsToDOF( i );

        // Set identity in boundary rows
        if (i == globalBox.first[0] || i == globalBox.last[0]) {
            //localCSRData[dof] = localCSR_identity( dof );
            localCSRData[dof] = localCSR_scaledIdentity( dof, stencil[1] );
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
    for ( auto j = localBox.first[1]; j <= localBox.last[1]; j++ ) {
        for ( auto i = localBox.first[0]; i <= localBox.last[0]; i++ ) {
        
            // The current row
            size_t dof = gridIndsToDOF( i, j );

            // Set identity in boundary rows
            if (j == globalBox.first[1] || j == globalBox.last[1] || i == globalBox.first[0] || i == globalBox.last[0]) {
                //localCSRData[dof] = localCSR_identity( dof );
                localCSRData[dof] = localCSR_scaledIdentity( dof, stencil[4] );
                continue;
            }
            
            // Copy of stencil
            std::vector<double> vals = stencil; 
            // Column indices, ordered consistently with the stencil
            std::vector<size_t> cols = { 
                gridIndsToDOF( i-1, j-1 ), // SW
                gridIndsToDOF( i ,  j-1 ), // S
                gridIndsToDOF( i+1, j-1 ), // SE
                gridIndsToDOF( i-1, j   ), // W
                dof,                       // O
                gridIndsToDOF( i+1, j   ), // E
                gridIndsToDOF( i-1, j+1 ), // NW
                gridIndsToDOF( i ,  j+1 ), // N
                gridIndsToDOF( i+1, j+1 )  // NW
            };

            localCSRData[dof] = { cols, vals };
        }
    }   

    return localCSRData;
}

/* Get CSR structure of anisotropic 3D Laplacian */
std::map<size_t, colsDataPair> PoissonOp::getCSRData3D( ) {    

    // Get 7-point stencil
    auto stencil = getStencil3D( );

    // Get local grid index box w/ zero ghosts
    auto localBox  = getLocalNodeBox( d_BoxMesh );
    auto globalBox = getGlobalNodeBox( d_BoxMesh );

    // Create a map from the DOF to a pair a vectors
    // Map from a DOF to vector of col inds and associated data
    std::map<size_t, colsDataPair> localCSRData;

    // Iterate over local box
    for ( auto k = localBox.first[2]; k <= localBox.last[2]; k++ ) {
        for ( auto j = localBox.first[1]; j <= localBox.last[1]; j++ ) {
            for ( auto i = localBox.first[0]; i <= localBox.last[0]; i++ ) {

                // The current row
                size_t dof = gridIndsToDOF( i, j, k );

                // Set identity in boundary rows
                if (k == globalBox.first[2] || k == globalBox.last[2] || j == globalBox.first[1] || j == globalBox.last[1] || i == globalBox.first[0] || i == globalBox.last[0]) {
                    //localCSRData[dof] = localCSR_identity( dof );
                    localCSRData[dof] = localCSR_scaledIdentity( dof, stencil[9] );
                    continue;
                }

                // Copy of stencil
                std::vector<double> vals = stencil; 
                // Column indices, ordered consistently with the stencil
                // DS, DW, DO, DE, DN, // z = -1 plane
                // MSW, MS, MSE, MW, MO, ME, MNW, MN, MNE, // z = 0 plane
                // US, UW, UO, UE, UN, // z = +1 plane
                std::vector<size_t> cols = { 
                    gridIndsToDOF( i,   j-1, k-1 ), // DS
                    gridIndsToDOF( i-1, j,   k-1 ), // DW
                    gridIndsToDOF( i,   j,   k-1 ), // DO
                    gridIndsToDOF( i+1, j,   k-1 ), // DE
                    gridIndsToDOF( i,   j+1, k-1 ), // DN
                    //
                    gridIndsToDOF( i-1, j-1, k   ), // MSW
                    gridIndsToDOF( i ,  j-1, k   ), // MS
                    gridIndsToDOF( i+1, j-1, k   ), // MSE
                    gridIndsToDOF( i-1, j,   k   ), // MW
                    gridIndsToDOF( i  , j,   k   ), // MO
                    gridIndsToDOF( i+1, j,   k   ), // ME
                    gridIndsToDOF( i-1, j+1, k   ), // MNW
                    gridIndsToDOF( i ,  j+1, k   ), // MN
                    gridIndsToDOF( i+1, j+1, k   ), // MNE
                    //
                    gridIndsToDOF( i,   j-1, k+1 ), // US
                    gridIndsToDOF( i-1, j,   k+1 ), // UW
                    gridIndsToDOF( i,   j,   k+1 ), // UO
                    gridIndsToDOF( i+1, j,   k+1 ), // UE
                    gridIndsToDOF( i,   j+1, k+1 )  // UN
                };
                localCSRData[dof] = { cols, vals };
            }
        } 
    }

    return localCSRData;
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
    } else if ( meshDim == 3 ) {
        localCSRData = this->getCSRData3D( );
    }

    // Create Lambda to return col inds from a given row ind
    auto getColumnIDs = [&](int row) { return localCSRData[row].cols; };
    // Create CSR matrix
    std::shared_ptr<AMP::LinearAlgebra::Matrix> A = AMP::LinearAlgebra::createMatrix( inVec, outVec, "CSRMatrix", getColumnIDs );
    fillMatWithLocalCSRData( A, d_DOFMan, localCSRData );
    
    // Finalize A
    A->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
    // Print statistics of A
    #if 0
    size_t nGlobalRows = A->numGlobalRows();
    size_t nLocalRows  = A->numLocalRows();
    std::cout << "Matrix A: #global rows=" << nGlobalRows << ". #local rows on rank " << A->getComm().getRank() << " (of " << A->getComm().getSize() << ") = " << nLocalRows
              << std::endl;
    #endif
    return A;
}


/* --------------------------------------------------------------
    Implementation of a class describing a diffusion equation 
--------------------------------------------------------------- */
/* Abstract base class representing a linear diffusion problem:
        - grad \dot ( D * grad u ) = s.

    This class is built on top of a PoissonOperator which is a discretization of the term 
        - grad \dot ( D * grad u )
*/
class DiffusionModel {

public:
    bool d_exactSolutionAvailable = false;

    // PoissonOperator on which this class is built
    std::shared_ptr<PoissonOp> d_PoissonOp;

    DiffusionModel( std::shared_ptr<PoissonOp> PoissonOp_) : 
        d_PoissonOp( PoissonOp_ ){ 
        AMP_INSIST( d_PoissonOp,  "Non-null PoissonOp required!" );
    }

    /* Pure virtual functions */
    virtual double sourceTerm( AMP::Mesh::MeshElement &node ) = 0;

    /* Virtual functions */
    virtual double exactSolution( AMP::Mesh::MeshElement & ) {
        AMP_ERROR( "Base class cannot provide an implementation of this function" );
    }
};


/* ------------------------------------------------
    Class implementing basic diffusion equation 
------------------------------------------------- */
/* A zero source term, and no exact solution */
class BasicDiffusionModel : public DiffusionModel {

public:
    // Call base class' constructor
    BasicDiffusionModel( std::shared_ptr<PoissonOp> PoissonOp_) : DiffusionModel( PoissonOp_ ) { }
    
    // Implementation of pure virtual function
    double sourceTerm( AMP::Mesh::MeshElement & ) override {
        return 0.0;
    }
};


/* -------------------------------------------------------
    Class implementing manufactured diffusion equation 
------------------------------------------------------- */
/* A source term and corresponding exact solution are provided */
class ManufacturedDiffusionModel : public DiffusionModel {

public:

    // Call base class' constructor
    ManufacturedDiffusionModel( std::shared_ptr<PoissonOp> PoissonOp_) : DiffusionModel( PoissonOp_ ) { 
        // Set flag indicating this class does provide an implementation of exactSolution
        d_exactSolutionAvailable = true;
    }
    
    // Implementation of pure virtual function
    // Dimesion-agnostic wrapper around the exact source term functions
    double sourceTerm( AMP::Mesh::MeshElement &node ) override {
        int meshDim = d_PoissonOp->getMesh()->getDim();
        if ( meshDim == 1 ) {
            double x = ( node.coord() )[0];
            return sourceTerm_( x );
        } else if ( meshDim == 2 ) {
            double x = ( node.coord() )[0];
            double y = ( node.coord() )[1];
            return sourceTerm_( x, y );
        } else if ( meshDim == 3 ) {
            double x = ( node.coord() )[0];
            double y = ( node.coord() )[1];
            double z = ( node.coord() )[2];
            return sourceTerm_( x, y, z );
        } else {
            AMP_ERROR( "Invalid dimension" );
        }
    }

    // Dimesion-agnostic wrapper around the exact solution functions
    double exactSolution( AMP::Mesh::MeshElement &node ) override {
        int meshDim = d_PoissonOp->getMesh()->getDim();
        if ( meshDim == 1 ) {
            double x = ( node.coord() )[0];
            return exactSolution_( x );
        } else if ( meshDim == 2 ) {
            double x = ( node.coord() )[0];
            double y = ( node.coord() )[1];
            return exactSolution_( x, y );
        } else if ( meshDim == 3 ) {
            double x = ( node.coord() )[0];
            double y = ( node.coord() )[1];
            double z = ( node.coord() )[2];
            return exactSolution_( x, y, z );
        } else {
            AMP_ERROR( "Invalid dimension" );
        }
    }

    // Note that the nodes are atually on the boundary, so just retrun the exact solution evaluate there. 
    double getDirichletValue( int, AMP::Mesh::MeshElement & node ) {
        int meshDim = d_PoissonOp->getMesh()->getDim();
        if ( meshDim == 1 ) {
            double x = ( node.coord() )[0];
            return exactSolution_( x );
        } else if ( meshDim == 2 ) {
            double x = ( node.coord() )[0];
            double y = ( node.coord() )[1];
            return exactSolution_( x, y );
        } else if ( meshDim == 3 ) {
            double x = ( node.coord() )[0];
            double y = ( node.coord() )[1];
            double z = ( node.coord() )[2];
            return exactSolution_( x, y, z );
        } else {
            AMP_ERROR( "Invalid dimension" );
        }
    }

    
private:
    // Exact solution, and corresponding source term
    // 1D
    double exactSolution_(double x) {
        return std::sin(2*M_PI*x - 0.65400000000000003);
    }
    double sourceTerm_(double x) {
        auto d_c   = d_PoissonOp->d_c;
        double cxx = d_c.xx;
        return 4*std::pow(M_PI, 2)*cxx*std::sin(2*M_PI*x - 0.65400000000000003);
    }
    // 2D
    double exactSolution_(double x, double y) {
        return std::sin(2*M_PI*x - 0.32500000000000001)*std::sin(4*M_PI*y + 0.98699999999999999);
    }
    double sourceTerm_(double x, double y) {
        auto d_c   = d_PoissonOp->d_c;
        double cxx = d_c.xx, cyy = d_c.yy, cxy = d_c.xy; 
        return 4*std::pow(M_PI, 2)*(cxx*std::sin(2*M_PI*x - 0.32500000000000001)*std::sin(4*M_PI*y + 0.98699999999999999) - 2*cxy*std::cos(2*M_PI*x - 0.32500000000000001)*std::cos(4*M_PI*y + 0.98699999999999999) + 4*cyy*std::sin(2*M_PI*x - 0.32500000000000001)*std::sin(4*M_PI*y + 0.98699999999999999));
    }
    // 3D
    double exactSolution_(double x, double y, double z) {
        return std::sin(2*M_PI*x - 0.98699999999999999)*std::sin(4*M_PI*y - 0.22500000000000001)*std::sin(6*M_PI*z - 0.47799999999999998);
    }
    double sourceTerm_(double x, double y, double z) {
        auto d_c   = d_PoissonOp->d_c;
        double cxx = d_c.xx, cyy = d_c.yy, czz = d_c.zz, cxy = d_c.xy, cxz = d_c.xz, cyz = d_c.yz; 
        return 4*std::pow(M_PI, 2)*(cxx*std::sin(2*M_PI*x - 0.98699999999999999)*std::sin(4*M_PI*y - 0.22500000000000001)*std::sin(6*M_PI*z - 0.47799999999999998) - 2*cxy*std::sin(6*M_PI*z - 0.47799999999999998)*std::cos(2*M_PI*x - 0.98699999999999999)*std::cos(4*M_PI*y - 0.22500000000000001) - 3*cxz*std::sin(4*M_PI*y - 0.22500000000000001)*std::cos(2*M_PI*x - 0.98699999999999999)*std::cos(6*M_PI*z - 0.47799999999999998) + 4*cyy*std::sin(2*M_PI*x - 0.98699999999999999)*std::sin(4*M_PI*y - 0.22500000000000001)*std::sin(6*M_PI*z - 0.47799999999999998) - 6*cyz*std::sin(2*M_PI*x - 0.98699999999999999)*std::cos(4*M_PI*y - 0.22500000000000001)*std::cos(6*M_PI*z - 0.47799999999999998) + 9*czz*std::sin(2*M_PI*x - 0.98699999999999999)*std::sin(4*M_PI*y - 0.22500000000000001)*std::sin(6*M_PI*z - 0.47799999999999998));
    }
};



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
    auto solver_db = input_db->getDatabase( "LinearSolver" );

    AMP_INSIST( PDE_db,    "A PDE database must be provided" );
    AMP_INSIST( solver_db, "A LinearSolver database must be provided" );

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

    auto myPoissonOp = std::make_shared<PoissonOp>( OpParameters );    
    
    // auto A = myPoissonOp->getMatrix();
    // AMP::IO::AsciiWriter matWriter;
    // matWriter.registerMatrix( A );
    // matWriter.writeFile( "Aout", 0  );



    /****************************************************************
    * Create diffusion equation model                               *
    ****************************************************************/
    std::shared_ptr<DiffusionModel> myDiffusionModel;
    auto model = PDE_db->getWithDefault<std::string>( "model", "" );

    if ( model == "basic" ) {
        AMP::pout << "Using BASIC diffusion model: No exact solution available" << std::endl;
        auto myDiffusionModel_ = std::make_shared<BasicDiffusionModel>( myPoissonOp );
        myDiffusionModel        = myDiffusionModel_;

    } else if ( model == "manufactured" ) {
        AMP::pout << "Using MANUFACTURED diffusion model: Exact solution is available" << std::endl;
        auto myDiffusionModel_ = std::make_shared<ManufacturedDiffusionModel>( myPoissonOp );
        
        // Point the Dirichlet BC values in the PoissonOp to those given by the manufactured problem
        myPoissonOp->setDirichletFunction( std::bind( &ManufacturedDiffusionModel::getDirichletValue, &( *myDiffusionModel_ ), std::placeholders::_1, std::placeholders::_2 ) );
        myDiffusionModel = myDiffusionModel_;        
    } else  {
        AMP_ERROR( "Model not recognised" );
    }

    // Create hassle-free wrappers around source term and exact solution
    auto PDESourceFun = std::bind( &DiffusionModel::sourceTerm, &( *myDiffusionModel ), std::placeholders::_1 );
    auto uexactFun    = std::bind( &DiffusionModel::exactSolution, &( *myDiffusionModel ), std::placeholders::_1 );


    /****************************************************************
    * Set up relevant vectors over the mesh                         *
    ****************************************************************/
    // Create required vectors over the mesh
    auto unumVec    = myPoissonOp->getRightVector();
    auto uexactVec  = myPoissonOp->getRightVector();
    auto rexactVec  = myPoissonOp->getRightVector();
    auto fsourceVec = myPoissonOp->getRightVector();

    // Set RHS vector
    myPoissonOp->fillWithFunction( fsourceVec, PDESourceFun );
    myPoissonOp->applyDirichletCorrectionToVector( fsourceVec );
    //fsourceVec->setToScalar( 1.0 );
    fsourceVec->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
    
    // Initialize unum to random values
    unumVec->setRandomValues();
    unumVec->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
    

    /****************************************************************
    * Compute discrete residual norm on continuous solution (truncation error) *
    ****************************************************************/
    /* Compare numerical solution with manufactured solution */
    if ( myDiffusionModel->d_exactSolutionAvailable ) {
        // Set exact solution 
        myPoissonOp->fillWithFunction( uexactVec, uexactFun );
        uexactVec->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );

        AMP::pout << "\nDiscrete residual of continuous manufactured solution: ";
        myPoissonOp->residual( fsourceVec, uexactVec, rexactVec );
        auto rnorms = getDiscreteNorms( PDE_db->getScalar<double>( "h" ), rexactVec );
        // Print residual norms
        AMP::pout << "||r|| = (" << rnorms[0] << ", " << rnorms[1] << ", " << rnorms[2] << ")" << std::endl << std::endl;

        #if 0
        for ( size_t row = 0; row < myPoissonOp->getMatrix()->numLocalRows(); row++ ) {
            auto f = fsourceVec->getValueByLocalID( row );
            auto u = uexactVec->getValueByLocalID( row );
            auto r = rexactVec->getValueByLocalID( row );
            std::cout << "i=" << row << ": f,u,r=" << f << "," << u << "," << r << std::endl;
        }
        #endif
    }

    
    /****************************************************************
    * Construct linear solver of the LinearOperator and apply it    *
    ****************************************************************/
    // Get the linear solver for operator myPoissonOp
    auto solverStratParams = std::make_shared<AMP::Solver::SolverStrategyParameters>( solver_db );
    solverStratParams->d_comm      = comm;
    solverStratParams->d_pOperator = myPoissonOp;
    auto linearSolver = AMP::Solver::SolverFactory::create( solverStratParams );
    
    // Use zero initial iterate and apply solver
    //linearSolver->setZeroInitialGuess( true );
    linearSolver->apply(fsourceVec, unumVec);
    
    // Compute disretization error
    if ( myDiffusionModel->d_exactSolutionAvailable ) {
        AMP::pout << "\nDiscretization error post linear solve: "; 
        auto e = uexactVec->clone();
        e->axpy(-1.0, *unumVec, *uexactVec); 
        auto enorms = getDiscreteNorms( PDE_db->getScalar<double>( "h" ), e );
        AMP::pout << "||e|| = (" << enorms[0] << ", " << enorms[1] << ", " << enorms[2] << ")" << std::endl;
    }

    // No specific solution is implemented for this problem, so this will just check that the solver converged. 
    checkConvergence( linearSolver.get(), input_db, input_file, *ut );

}
// end of driver()



/*  The input file must contain a "PDE" database, which has 

    dim : the dimension of the problem (1, 2, or 3)
    n   : the number of mesh points (plus 1) in each grid dimension
    
    uj  : j=1,2 for 1D; j=1,2,3,4 for 2D; j=1,2,3,4,5,6 for 3D are Dirichlet boundary conditions on boundary j. Note that the Dirichlet boundaries are not eliminated from the system, with the matrix having scaled identity entries in the corresponding rows. Specifically, each boundary row is scaled by the diagonal connection from the stencil; this is to help with ill-conditioning of the resulting matrix.

    model : either 1. "basic" or 2. "manufactured"
        1. A zero source term is used in the PDE
        2. A manufactured solution are corresponinding source term are available; note that in this case the Dirichlet constants uj in the input file (while still required) are ignored, with the boundary conditions set to the manufactured solution on the boundary.

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
    exeNames.emplace_back( "testLinearSolvers-PoissonFD-BoomerAMG-CG" );

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
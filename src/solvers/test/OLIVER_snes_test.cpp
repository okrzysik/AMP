#include "AMP/IO/PIO.h"
#include "AMP/utils/AMPManager.h"

#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/Vector.h"

#include "AMP/discretization/boxMeshDOFManager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshID.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/mesh/MeshElement.h"
#include "AMP/mesh/structured/BoxMesh.h"

#include "AMP/operators/Operator.h"
#include "AMP/operators/OperatorParameters.h"

#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/solvers/SolverFactory.h"
#include "AMP/solvers/petsc/PetscSNESSolver.h"

#include <iostream>
#include <iomanip>


// Implement the operator A(x) = x^2 - 1 
class QuadraticOper : public AMP::Operator::Operator {

public: 

    // Call base class constructor
    QuadraticOper(std::shared_ptr<const AMP::Operator::OperatorParameters> params) : AMP::Operator::Operator( params ) { };

    // Implementation of pure virtual function
    std::string type() const { return "QuadraticOper"; };
    
    // Implementation of pure virtual function
    // Given input x, evaluate A(x) = x^2 - 1
    void apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> x_vec, std::shared_ptr<AMP::LinearAlgebra::Vector> Ax_vec) 
    {
        std::cout << "\t\tQuadraticOper: entering apply..." << std::endl;
        Ax_vec->multiply(*x_vec, *x_vec);
        Ax_vec->addScalar(*Ax_vec, -1.0);
        std::cout << "\t\tQuadraticOper: exiting apply..." << std::endl;
    }
};


/* Create a basic BoxMesh */
static std::shared_ptr<AMP::Mesh::BoxMesh> createBoxMesh( AMP::AMP_MPI comm, int n )
{   
    auto mesh_db = std::make_shared<AMP::Database>( "Mesh" );
    mesh_db->putScalar<int>( "dim", 1 );
    mesh_db->putScalar<std::string>( "MeshName", "AMP::cube" );
    mesh_db->putScalar<std::string>( "Generator", "cube" );
    mesh_db->putVector<int>( "Size", { n } ); // mesh has n+1 points
    mesh_db->putVector<double>( "Range", { 0.0, 1.0 } );

    // Create MeshParameters
    auto mesh_params = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    mesh_params->setComm( comm );
    
    // Create Mesh
    static std::shared_ptr<AMP::Mesh::BoxMesh> mesh = AMP::Mesh::BoxMesh::generate( mesh_params );
    return mesh;
}


/* Print the vector x that *approximately* solves A(x) = rhs_val */
void print_solution( std::shared_ptr<AMP::LinearAlgebra::Vector> x_vec, double x_exact, double rhs_val, AMP::Mesh::MeshIterator meshIt, std::shared_ptr<AMP::Discretization::DOFManager> dof ) {
    std::cout << std::endl;
    std::cout << std::setprecision(3) << std::showpos;
    std::cout << "x*_i s.t. [A(x*)]_i=" << rhs_val << " (exact sol=" << x_exact << "):" << std::endl;
    for (auto elem = meshIt.begin(); elem != meshIt.end(); elem++ ) {
        auto id = elem->globalID();
        std::vector<size_t> dofs;
        dof->getDOFs(id, dofs);
        std::cout << "i=" << dofs[0] << ": ";
        std::cout << x_vec->getValueByGlobalID(dofs[0]) << ", ";
        std::cout << std::endl;
    }
}


/* Assemble the quadratic operator A(x) = x^2 - 1, and solve A(x) = b for b = 0, b = -1, and b=NULL */
void driver(AMP::AMP_MPI comm, int n ) {

    /****************************************************************
    * Create a mesh, and DOFManager over it                         *
    ****************************************************************/
    std::shared_ptr<AMP::Mesh::BoxMesh> mesh = createBoxMesh( comm, n ); 
    auto meshIt = mesh->getIterator( AMP::Mesh::GeomType::Vertex );
    std::shared_ptr<AMP::Discretization::DOFManager> dof = AMP::Discretization::simpleDOFManager::create(mesh, AMP::Mesh::GeomType::Vertex, 1, 1);

    /****************************************************************
    * Create nonlinear operator over the mesh                       *
    ****************************************************************/
    auto Op_db = std::make_shared<AMP::Database>( "Op_db" );
    auto Op_params = std::make_shared<AMP::Operator::OperatorParameters>( Op_db );
    Op_params->d_Mesh = mesh;
    auto myOper = std::make_shared<QuadraticOper>( Op_params ); 

    /****************************************************************
    * Set up relevant vectors over the mesh                         *
    ****************************************************************/
    auto x_var  = std::make_shared<AMP::LinearAlgebra::Variable>( "x" );
    auto Ax_var = std::make_shared<AMP::LinearAlgebra::Variable>( "Ax" );
    auto x_vec  = AMP::LinearAlgebra::createVector( dof, x_var );
    auto Ax_vec = AMP::LinearAlgebra::createVector( dof, Ax_var );
    x_vec->setRandomValues();

    // Confirm action of operator works and looks as it should
    #if 1
        std::cout << std::endl;
        AMP::pout << "----------------------------------------" << std::endl;
        std::cout << std::setprecision(3) << std::showpos;
        std::cout << "Testing operator apply: i, x_i, [A(x)]_i:" << std::endl;
        myOper->apply( x_vec, Ax_vec );
        for (auto elem = meshIt.begin(); elem != meshIt.end(); elem++ ) {
            auto id = elem->globalID();
            std::vector<size_t> dofs;
            dof->getDOFs(id, dofs);
            std::cout << dofs[0] << ": ";
            std::cout << x_vec->getValueByGlobalID(dofs[0]) << ", ";
            std::cout << Ax_vec->getValueByGlobalID(dofs[0]) << " ";
            std::cout << std::endl;
        }
        AMP::pout << "----------------------------------------\n" << std::endl;
    #endif


    /****************************************************************
    * Create and apply a SNES solver                                *
    ****************************************************************/
    // Some of the code here is based on that in "solvers/test/mechanics/testPetscSNESSolver_NonlinearThermoMechanics_1.cpp"
    // HACK to prevent a double delete on Petsc Vec
    std::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver;

    auto nonlinearSolver_db = std::make_shared<AMP::Database>( "NonlinearSolver" );
    nonlinearSolver_db->putScalar<int>("print_info_level", 1);
    nonlinearSolver_db->putScalar<int>("max_iterations", 8);
    nonlinearSolver_db->putScalar<double>("max_error", 1e-8);
    std::string SNESOptions = "-snes_monitor -snes_type ls -snes_fd -ksp_type gmres -snes_linesearch_monitor";
    nonlinearSolver_db->putScalar<std::string>( "SNESOptions", SNESOptions );
    auto nonlinearSolverParams =
        std::make_shared<AMP::Solver::SolverStrategyParameters>( nonlinearSolver_db );
    
    nonlinearSolverParams->d_comm          = comm;
    nonlinearSolverParams->d_pOperator     = myOper;
    nonlinearSolverParams->d_pNestedSolver = nullptr;
    nonlinearSolverParams->d_pInitialGuess = nullptr;
    nonlinearSolver.reset( new AMP::Solver::PetscSNESSolver( nonlinearSolverParams ) );
    
    //nonlinearSolver->setZeroInitialGuess( false );
    nonlinearSolver->setZeroInitialGuess( true ); // Bug? This is seemingly being ignored...

    // Make rhs vectors of 0s, 1s, and null
    auto rhs_vec0 = x_vec->clone();
    rhs_vec0->zero();
    auto rhs_vec1 = x_vec->clone();
    rhs_vec1->setToScalar(-1.0);
    AMP::LinearAlgebra::Vector::shared_ptr rhs_vecNULL;

    // Solve A(x) = 0 --> x^2 - 1 = 0. For which the solution is x = +1 or x = -1 (we probably converge to the positive one if we initialize from a positive number)
    // Newton converges quadratically here...
    AMP::pout << "----------------------------------------" << std::endl;
    AMP::pout << "Starting Nonlinear Solve with RHS=0.0..." << std::endl;
    nonlinearSolver->apply( rhs_vec0, x_vec );
    print_solution( x_vec, 1.0, 0.0, meshIt, dof );
    AMP::pout << "Finished Nonlinear Solve with RHS=0.0..." << std::endl;
    AMP::pout << "----------------------------------------\n" << std::endl;

    // Solve A(x) = -1 --> x^2 - 1 = -1. For which the solution is x = 0
    // Newton converges linearly here... (consistent with the fact that the jacobian is 2*x, and this vanishes at the root)
    AMP::pout << "-----------------------------------------" << std::endl;
    AMP::pout << "Starting Nonlinear Solve with RHS=-1.0..." << std::endl;
    nonlinearSolver->apply( rhs_vec1, x_vec );
    print_solution( x_vec, 0.0, -1.0, meshIt, dof );
    AMP::pout << "Finished Nonlinear Solve with RHS=-1.0..." << std::endl;
    AMP::pout << "-----------------------------------------\n" << std::endl;

    // Solve A(x) = NULL
    AMP::pout << "-----------------------------------------" << std::endl;
    AMP::pout << "Starting Nonlinear Solve with RHS=NULL..." << std::endl;
    nonlinearSolver->apply( rhs_vecNULL, x_vec );
    //nonlinearSolver->apply( nullptr, x_vec );
    print_solution( x_vec, 1.0, 0.0, meshIt, dof ); // I assume null rhs is interpreted as a 0
    AMP::pout << "Finished Nonlinear Solve with RHS=NULL..." << std::endl;
    AMP::pout << "-----------------------------------------\n" << std::endl;

}
// end of driver()


/* Usage: >> mpirun -n 1 NST 2 */
int main( int argc, char **argv )
{

    AMP::AMPManager::startup( argc, argv );

    // Create a global communicator
    AMP::AMP_MPI comm( AMP_COMM_WORLD );

    // Unpack input
    //int n = atoi(argv[1]); // Number of DOFs is n+1, for n>1. 
    int n = 2;

    // Driver
    driver( comm, n );

    AMP::AMPManager::shutdown();

    return 0;
}
// end of main()
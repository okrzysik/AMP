#include "OxideTimeIntegrator.h"
#include "OxideModel.h"

#include "utils/AMP_MPI.h"
#include "utils/Utilities.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"
#include "vectors/MultiVariable.h"
#include "vectors/VectorSelector.h"


namespace AMP{
namespace TimeIntegrator{


/************************************************************************
* Constructor and destructor for TimeIntegrator.                        *
************************************************************************/
OxideTimeIntegrator::OxideTimeIntegrator( boost::shared_ptr<TimeIntegratorParameters> parameters )
{
   AMP_INSIST(parameters.get()!=NULL, "Null parameter");

   initialize(parameters);   
}

OxideTimeIntegrator::~OxideTimeIntegrator()
{
}


/************************************************************************
* Initialize the time integrator and problem                            *
************************************************************************/
void OxideTimeIntegrator::initialize( boost::shared_ptr<TimeIntegratorParameters> parameters )
{   
    d_current_time = 0.0;
    d_current_dt = 1.0;

    // Get the parameters
    boost::shared_ptr<OxideTimeIntegratorParameters> oxide_parameters = 
        boost::dynamic_pointer_cast<OxideTimeIntegratorParameters>( parameters );
    d_mesh = oxide_parameters->d_mesh;
    AMP_INSIST(d_mesh.get()!=NULL,"Oxide Time Integrator needs a mesh");
    AMP_INSIST((int)d_mesh->getGeomType()<d_mesh->getDim(),
        "Oxide mesh must be a surface mesh (dimension < physical dimension");
    AMP_INSIST(oxide_parameters->d_temp.get()!=NULL,"Oxide Time Integrator needs a temerature vector");
    AMP::LinearAlgebra::VS_Mesh meshSelector("temperature",d_mesh);
    d_temp = (oxide_parameters->d_temp)->select(meshSelector,"temperature");
    AMP_ASSERT(d_temp.get());
    std::vector<size_t> dofs;
    d_temp->getDOFManager()->getDOFs( d_mesh->getIterator(AMP::Mesh::Vertex,0)->globalID(), dofs );
    AMP_INSIST(dofs.size()==1,"Temperature vector must be a nodal scalar vector");
    double total_depth = oxide_parameters->depth;
    AMP_INSIST(total_depth>0,"Surface depth must be > 0");

    // Create the solution vector
    AMP::Discretization::DOFManager::shared_ptr DOF = AMP::Discretization::simpleDOFManager::create(d_mesh,AMP::Mesh::Vertex,1,1,true);
    AMP::LinearAlgebra::Variable::shared_ptr oxide_var( new AMP::LinearAlgebra::Variable("oxide") );
    AMP::LinearAlgebra::Variable::shared_ptr alpha_var( new AMP::LinearAlgebra::Variable("alpha") );
    boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> multivariable( new AMP::LinearAlgebra::MultiVariable("solution") );
    multivariable->add(oxide_var);
    multivariable->add(alpha_var);
    d_solution = AMP::LinearAlgebra::createVector( DOF, multivariable, true );
    d_oxide = d_solution->subsetVectorForVariable( oxide_var );
    d_alpha = d_solution->subsetVectorForVariable( alpha_var );
    d_solution->setToScalar(0.0);
    
    // Create the internal vectors for storing the data
    N_layer = std::vector<int>(3);
    N_layer[0] = 20;    // Number of zones in the oxide layer
    N_layer[1] = 20;    // Number of zones in the alpha layer
    N_layer[2] = 5;     // Number of zones in the zirconium layer
    int N_total = 0;
    for (size_t i=0; i<N_layer.size(); i++)
        N_total += N_layer[i];
    AMP::Discretization::DOFManager::shared_ptr DOF_d = AMP::Discretization::simpleDOFManager::create(d_mesh,AMP::Mesh::Vertex,0,N_layer.size(),true);
    AMP::Discretization::DOFManager::shared_ptr DOF_C = AMP::Discretization::simpleDOFManager::create(d_mesh,AMP::Mesh::Vertex,0,N_total,true);
    AMP::LinearAlgebra::Variable::shared_ptr d_var( new AMP::LinearAlgebra::Variable("depth") );
    AMP::LinearAlgebra::Variable::shared_ptr C_var( new AMP::LinearAlgebra::Variable("C") );
    depth = AMP::LinearAlgebra::createVector( DOF_d, d_var, true );
    conc  = AMP::LinearAlgebra::createVector( DOF_C, C_var, true );

    // Create the initial conditions (500K for 1 day)
    double *C0[10], *C1[10], D[10], Cb[20], x0[11], x1[11], v1[11], depth2[10];   // Allocate enough space for 10 layers
    C0[0] = new double[N_total];
    C1[0] = new double[N_total];
    for (size_t i=1; i<N_layer.size(); i++) {
        C0[i] = &C0[i-1][N_layer[i-1]];
        C1[i] = &C1[i-1][N_layer[i-1]];
    }
    OxideModel::get_equilibrium_concetration( 500, Cb );
    OxideModel::get_diffusion_coefficients( 500, D );
    for (size_t i=0; i<N_layer.size(); i++) {
        x0[i] = i*10.0e-7;        // Set the initial thickness to 10 nm
        for (int j=0; j<N_layer[i]; j++)
            C0[i][j] = Cb[2*i+0] + (Cb[2*i+1]-Cb[2*i+0])*(j+0.5)/((double)N_layer[i]);
    }
    x0[N_layer.size()] = 1e2*total_depth;   // Convert from m to cm
    OxideModel::integrateOxide( 86400, N_layer.size(), &N_layer[0], x0, Cb, C0, D, C1, x1, v1 );
    for (size_t i=0; i<N_layer.size(); i++)
        depth2[i] = x1[i+1]-x1[i];
    // Copy the initial solution to all points in the mesh
    AMP::Discretization::DOFManager::shared_ptr DOF_oxide = d_oxide->getDOFManager();
    AMP::Discretization::DOFManager::shared_ptr DOF_alpha = d_alpha->getDOFManager();
    AMP::Mesh::MeshIterator iterator = d_mesh->getIterator(AMP::Mesh::Vertex,0);
    for (size_t i=0; i<iterator.size(); i++) {
        AMP::Mesh::MeshElementID id = iterator->globalID();
        DOF_C->getDOFs( id, dofs );
        AMP_ASSERT((int)dofs.size()==N_total);
        conc->setLocalValuesByGlobalID( dofs.size(), &dofs[0], C1[0] );
        DOF_d->getDOFs( id, dofs );
        AMP_ASSERT(dofs.size()==N_layer.size());
        depth->setLocalValuesByGlobalID( dofs.size(), &dofs[0], depth2 );
        DOF_oxide->getDOFs( id, dofs );
        AMP_ASSERT(dofs.size()==1);
        d_oxide->setLocalValueByGlobalID( dofs[0], 1e-2*depth2[0] );   // Convert from cm to m
        DOF_alpha->getDOFs( id, dofs );
        AMP_ASSERT(dofs.size()==1);
        d_alpha->setLocalValueByGlobalID( dofs[0], 1e-2*depth2[1] );   // Convert from cm to m
        ++iterator;
    }
    d_solution->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

    // Free the temporary memory
    delete [] C0[0];
    delete [] C1[0];
}


/************************************************************************
* Reset the time integrator                                             *
************************************************************************/
void OxideTimeIntegrator::reset(boost::shared_ptr<TimeIntegratorParameters> parameters )
{
    AMP_ERROR("reset is not programmed for OxideTimeIntegrator");
}


/************************************************************************
* Reset the time integrator                                             *
************************************************************************/
int OxideTimeIntegrator::advanceSolution( const double dt, const bool first_step )
{
    d_current_time += dt;
    d_current_dt = dt;
    // Get the relavent DOF Managers
    AMP::Discretization::DOFManager::shared_ptr DOF_C = conc->getDOFManager();
    AMP::Discretization::DOFManager::shared_ptr DOF_d = depth->getDOFManager();
    AMP::Discretization::DOFManager::shared_ptr DOF_temp = d_temp->getDOFManager();
    AMP::Discretization::DOFManager::shared_ptr DOF_oxide = d_oxide->getDOFManager();
    AMP::Discretization::DOFManager::shared_ptr DOF_alpha = d_alpha->getDOFManager();
    // Allocate memory for the solve
    int N_total = 0;
    for (size_t i=0; i<N_layer.size(); i++)
        N_total += N_layer[i];
    double *C0[10], *C1[10], D[10], Cb[20], x0[11], x1[11], v1[11], depth2[10];   // Allocate enough space for 10 layers
    C0[0] = new double[N_total];
    C1[0] = new double[N_total];
    for (size_t i=1; i<N_layer.size(); i++) {
        C0[i] = &C0[i-1][N_layer[i-1]];
        C1[i] = &C1[i-1][N_layer[i-1]];
    }
    std::vector<size_t> dofs;
    // Loop through the points
    AMP::Mesh::MeshIterator iterator = d_mesh->getIterator(AMP::Mesh::Vertex,0);
    for (size_t i=0; i<iterator.size(); i++) {
        AMP::Mesh::MeshElementID id = iterator->globalID();
        // Get the current temperature
        DOF_temp->getDOFs( id, dofs );
        AMP_ASSERT(dofs.size()==1);
        double T = d_temp->getValueByGlobalID( dofs[0] );
        // Get the equilibrium concentrations and diffusion coefficients
        OxideModel::get_equilibrium_concetration( T, Cb );
        OxideModel::get_diffusion_coefficients( T, D );
        // Get the previous solution's concentration
        DOF_C->getDOFs( id, dofs );
        AMP_ASSERT((int)dofs.size()==N_total);
        conc->getValuesByGlobalID( dofs.size(), &dofs[0], C0[0] );
        // Get the previous solution's coordinates
        DOF_d->getDOFs( id, dofs );
        AMP_ASSERT(dofs.size()==N_layer.size());
        depth->getValuesByGlobalID( dofs.size(), &dofs[0], depth2 );
        x0[0] = 0.0;
        for (size_t i=0; i<N_layer.size(); i++)
            x0[i+1] = x0[i] + depth2[i];
        // Perform the time integration
        double dt2 = dt*3600*24;    // Convert from days to seconds
        AMP_ASSERT(dt2>=0.0);
        OxideModel::integrateOxide( dt2, N_layer.size(), &N_layer[0], x0, Cb, C0, D, C1, x1, v1 );
        for (size_t i=0; i<N_layer.size(); i++)
            depth2[i] = x1[i+1]-x1[i];
        // Save the results
        DOF_C->getDOFs( id, dofs );
        AMP_ASSERT((int)dofs.size()==N_total);
        conc->setLocalValuesByGlobalID( dofs.size(), &dofs[0], C1[0] );
        DOF_d->getDOFs( id, dofs );
        AMP_ASSERT(dofs.size()==N_layer.size());
        depth->setLocalValuesByGlobalID( dofs.size(), &dofs[0], depth2 );
        DOF_oxide->getDOFs( id, dofs );
        AMP_ASSERT(dofs.size()==1);
        d_oxide->setLocalValueByGlobalID( dofs[0], 1e-2*depth2[0] );    // Convert from cm to m
        DOF_alpha->getDOFs( id, dofs );
        AMP_ASSERT(dofs.size()==1);
        d_alpha->setLocalValueByGlobalID( dofs[0], 1e-2*depth2[1] );    // Convert from cm to m
        ++iterator;
    }
    // Update ghost values for the solution
    d_solution->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    // Free the temporary memory
    delete [] C0[0];
    delete [] C1[0];
    return 0;
}


/************************************************************************
* Check the solution                                                    *
************************************************************************/
bool OxideTimeIntegrator::checkNewSolution(void) const
{
    return true;
}


/************************************************************************
* Update the solution                                                   *
************************************************************************/
void OxideTimeIntegrator::updateSolution( void )
{
}


/************************************************************************
* Return time increment for next solution advance.                      *
************************************************************************/
double OxideTimeIntegrator::getNextDt(const bool good_solution)
{
    return 1e10;
}



}
}



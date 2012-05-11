#include "OxideTimeIntegrator.h"
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
    AMP_INSIST(total_depth>0!=NULL,"Surface depth must be > 0");

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
    AMP::Discretization::DOFManager::shared_ptr DOF_t = AMP::Discretization::simpleDOFManager::create(d_mesh,AMP::Mesh::Vertex,0,N_layer.size(),true);
    AMP::Discretization::DOFManager::shared_ptr DOF_C = AMP::Discretization::simpleDOFManager::create(d_mesh,AMP::Mesh::Vertex,0,N_total,true);
    AMP::LinearAlgebra::Variable::shared_ptr t_var( new AMP::LinearAlgebra::Variable("depth") );
    AMP::LinearAlgebra::Variable::shared_ptr C_var( new AMP::LinearAlgebra::Variable("C") );
    depth = AMP::LinearAlgebra::createVector( DOF_t, t_var, true );
    C     = AMP::LinearAlgebra::createVector( DOF_C, C_var, true );

    // Create the initial conditions
    
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



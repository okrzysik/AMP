#include "DTKMapOperator.h"
#include "DTKAMPField.h"

#include <DTK_FieldMultiVector.hpp>
#include <DTK_ConsistentInterpolationOperator.hpp>

#include <functional>
#include <DTK_BasicEntityPredicates.hpp>

namespace AMP {
namespace Operator {

//---------------------------------------------------------------------------//
// Constructor
DTKMapOperator::DTKMapOperator( const AMP::shared_ptr<OperatorParameters> &params )
: d_dtk_operator(nullptr), d_mapOnThisProc(false) 
{
  // Get the operator parameters.
  AMP::shared_ptr<DTKMapOperatorParameters> dtk_op_params =
    AMP::dynamic_pointer_cast<DTKMapOperatorParameters>( params );
  AMP_ASSERT( dtk_op_params );

  int inComm                             = -1;
  if ( dtk_op_params->d_domain_mesh || dtk_op_params->d_range_mesh ){
    d_mapOnThisProc = true;
    inComm = 1;
  }

  // Create a comm spanning the meshes
  d_comm = dtk_op_params->d_globalComm.split( inComm );
  if(d_mapOnThisProc)
    d_TeuchosComm = Teuchos::rcp( new Teuchos::MpiComm<int>( d_comm.getCommunicator() ) );

  // For now, select all elements in a mesh. This needs to be refactored to
  // select elements based on criteria in the parameters.
  DataTransferKit::PredicateFunction select_all = 
    [=]( DataTransferKit::Entity e ) { return true; };

  // Create DTK domain mesh objects.
  d_domain_mesh = std::make_shared<DTKAMPMeshManager>(
      dtk_op_params->d_domain_mesh, dtk_op_params->d_domain_dofs, select_all );

  // Create DTK range mesh objects.
  d_range_mesh = std::make_shared<DTKAMPMeshManager>(
      dtk_op_params->d_range_mesh, dtk_op_params->d_range_dofs, select_all );
  
}

//---------------------------------------------------------------------------//
//! Apply function.
void DTKMapOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
    AMP::LinearAlgebra::Vector::shared_ptr r)
{
  if(d_mapOnThisProc){
    // Create fields.
    AMP::LinearAlgebra::Vector::shared_ptr u_non_const =
      std::const_pointer_cast<AMP::LinearAlgebra::Vector>(u);
    Teuchos::RCP<DataTransferKit::Field> u_field = 
      Teuchos::rcp( new DTKAMPField(u_non_const) );
    Teuchos::RCP<DataTransferKit::Field> r_field = 
      Teuchos::rcp( new DTKAMPField(r) );

    // Create vectors from the field.
    Teuchos::RCP<DataTransferKit::FieldMultiVector::Base> u_vector
      = Teuchos::rcp( 
          new DataTransferKit::FieldMultiVector(d_TeuchosComm, u_field) );
    Teuchos::RCP<DataTransferKit::FieldMultiVector::Base> r_vector
      = Teuchos::rcp( 
          new DataTransferKit::FieldMultiVector(d_TeuchosComm, r_field) );

    // Lazy evaluate the map operator.
    if ( !d_dtk_operator )
    {
      // Build the operator.
      Teuchos::ParameterList dtk_parameters;
      dtk_parameters.sublist("Search");
      dtk_parameters.sublist("Search").set("Track Missed Range Entities", true);
      dtk_parameters.sublist("Consistent Interpolation");
      d_dtk_operator = AMP::shared_ptr<DataTransferKit::MapOperator>(
          new DataTransferKit::ConsistentInterpolationOperator( 
            u_vector->getMap(), r_vector->getMap(), dtk_parameters ) );

      std::fstream fout;
      if(Teuchos::nonnull(d_range_mesh->functionSpace()->entitySet())){
        std::string filename = "debug_DTKMapOp_" + AMP::Utilities::intToString( d_comm.getRank() );
        fout.open( filename.c_str(), std::fstream::out | std::fstream::app  );
        fout<< " Iterator size "<<d_range_mesh->functionSpace()->entitySet()->entityIterator(0).size()<<std::endl;
        DataTransferKit::LocalEntityPredicate pred(d_range_mesh->functionSpace()->entitySet()->communicator()->getRank());
        fout<< " Iterator size with predicate "<<d_range_mesh->functionSpace()->entitySet()->entityIterator(0,pred).size()<<std::endl;
        fout.close();
      }
 
      // Setup the operator.
      d_dtk_operator->setup(
          d_domain_mesh->functionSpace(), d_range_mesh->functionSpace() );
      // get missed range entities
      std::cout<<" Rank "<<d_comm.getRank()<< " Missed Entities "<<std::dynamic_pointer_cast<DataTransferKit::ConsistentInterpolationOperator>(d_dtk_operator)->getMissedRangeEntityIds()<<std::endl; 
    }

    // Do the apply.
    d_dtk_operator->apply( *u_vector, *r_vector );
  }
}

//---------------------------------------------------------------------------//
}
}

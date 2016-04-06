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
: d_dtk_operator(nullptr)
{
  // Get the operator parameters.
  AMP::shared_ptr<DTKMapOperatorParameters> dtk_op_params =
    AMP::dynamic_pointer_cast<DTKMapOperatorParameters>( params );
  AMP_ASSERT( dtk_op_params );

  int inComm                             = -1;
  if ( dtk_op_params->d_domain_mesh || dtk_op_params->d_range_mesh )
    inComm = 1;

  // Create a comm spanning the meshes
  d_comm = dtk_op_params->d_globalComm.split( inComm );
  if(!d_comm.isNull())
    d_TeuchosComm = Teuchos::rcp( new Teuchos::MpiComm<int>( d_comm.getCommunicator() ) );

  // For now, select all elements in a mesh. This needs to be refactored to
  // select elements based on criteria in the parameters.
  DataTransferKit::PredicateFunction select_all = 
    [=]( DataTransferKit::Entity e ) { return true; };

  // Create DTK domain mesh objects.
  d_domain_mesh = std::make_shared<DTKAMPMeshManager>(
      dtk_op_params->d_domain_mesh, dtk_op_params->d_domain_dofs, select_all );
  if(dtk_op_params->d_domain_mesh){    
    DataTransferKit::LocalEntityPredicate localPre (0);
    std::cout << "domain num local element "  <<dtk_op_params->d_domain_mesh->numLocalElements(AMP::Mesh::Volume) << std::endl;
    std::cout << "domain size "<< d_domain_mesh->entitySet()->entityIterator(3, localPre.getFunction() ).size() << " Rank " << d_comm.getRank()  << std::endl;
  } else{
    std::cout << "no domain mesh " << std::endl;
  }

  // Create DTK range mesh objects.
  d_range_mesh = std::make_shared<DTKAMPMeshManager>(
      dtk_op_params->d_range_mesh, dtk_op_params->d_range_dofs, select_all );
  if(dtk_op_params->d_range_mesh){
    DataTransferKit::LocalEntityPredicate localPre (d_comm.getRank());
    std::cout << "range num local vertex "  <<dtk_op_params->d_range_mesh->numLocalElements(AMP::Mesh::Vertex) << std::endl;
    std::cout << "range size "<< d_range_mesh->entitySet()->entityIterator(0, localPre.getFunction() ).size()  << " Rank " << d_comm.getRank() << std::endl;
  }else{
    std::cout << "no range mesh " << std::endl;
  }  
}

//---------------------------------------------------------------------------//
//! Apply function.
void DTKMapOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
    AMP::LinearAlgebra::Vector::shared_ptr r)
{
  if(!d_comm.isNull()){
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
      dtk_parameters.sublist("Consistent Interpolation");
      d_dtk_operator = AMP::shared_ptr<DataTransferKit::MapOperator>(
          new DataTransferKit::ConsistentInterpolationOperator( 
            u_vector->getMap(), r_vector->getMap(), dtk_parameters ) );

      // Setup the operator.
      d_dtk_operator->setup(
          d_domain_mesh->functionSpace(), d_range_mesh->functionSpace() );
    }

    // Do the apply.
    d_dtk_operator->apply( *u_vector, *r_vector );
  }
}

//---------------------------------------------------------------------------//
}
}

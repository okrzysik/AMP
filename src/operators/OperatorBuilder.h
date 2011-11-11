#ifndef included_AMP_OperatorBuilder
#define included_AMP_OperatorBuilder

#include "ampmesh/Mesh.h"
#include "utils/InputDatabase.h"
#include "operators/Operator.h"
#include "ElementPhysicsModel.h"
#include "ElementPhysicsModelFactory.h"
#include "operators/boundary/BoundaryOperator.h"

namespace AMP {
namespace Operator {

class OperatorBuilder{
 public:
  OperatorBuilder(){}
  ~OperatorBuilder(){}

  static boost::shared_ptr<Operator> createOperator(boost::shared_ptr<OperatorParameters>  in_params);

  static boost::shared_ptr<Operator>  createOperator( AMP::Mesh::Mesh::shared_ptr  mesh,
						      std::string operatorName,
						      boost::shared_ptr<AMP::Database>  input_db,
						      boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel,
						      boost::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> localModelFactory = 	boost::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> ());
  
  static boost::shared_ptr<Operator>  createOperator( AMP::Mesh::Mesh::shared_ptr  mesh1,
						      AMP::Mesh::Mesh::shared_ptr  mesh2,
						      boost::shared_ptr<AMP::Database>  input_db);
  
  static boost::shared_ptr<BoundaryOperator>
  createColumnBoundaryOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
				std::string boundaryOperatorName,
				boost::shared_ptr<AMP::InputDatabase> input_db,
				AMP::Operator::Operator::shared_ptr volumeOperator,
				boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel,
				boost::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> localModelFactory );
  
  static boost::shared_ptr<BoundaryOperator>
  createBoundaryOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
			  std::string boundaryOperatorName,
			  boost::shared_ptr<AMP::InputDatabase> input_db,
			  AMP::Operator::Operator::shared_ptr volumeOperator,
			  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel,
			  boost::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> localModelFactory = boost::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> () );
  
 protected:
  
  static boost::shared_ptr<Operator> createFlowFrapconOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                boost::shared_ptr<AMP::InputDatabase> input_db);
  
  static boost::shared_ptr<Operator> createFlowFrapconJacobian( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                boost::shared_ptr<AMP::InputDatabase> input_db);
  
  static boost::shared_ptr<Operator> createNeutronicsRhsOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                      boost::shared_ptr<AMP::InputDatabase> input_db);
  
  static boost::shared_ptr<Operator> createVolumeIntegralOperator(AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                  boost::shared_ptr<AMP::InputDatabase> input_db,
                                  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel);

  static boost::shared_ptr<Operator> createLinearConsMomentumGalWFOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                    boost::shared_ptr<AMP::InputDatabase> input_db,
                                    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel);
  
  static boost::shared_ptr<Operator> createLinearConsMassGalWFOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                    boost::shared_ptr<AMP::InputDatabase> input_db,
                                    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel);
  
  static boost::shared_ptr<Operator> createMassLinearFEOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                 boost::shared_ptr<AMP::InputDatabase> input_db,
                                 boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel);
  
  static boost::shared_ptr<Operator> createLinearDiffusionOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                    boost::shared_ptr<AMP::InputDatabase> input_db,
                                    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel);
  
  static boost::shared_ptr<Operator> createNonlinearDiffusionOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                       boost::shared_ptr<AMP::InputDatabase> input_db,
                                       boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel);
  
  static boost::shared_ptr<Operator> createNonlinearFickSoretOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
								       std::string operatorName,
								       boost::shared_ptr<AMP::InputDatabase> input_db,
								       boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel,
								       boost::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> localModelFactory );

  static boost::shared_ptr<Operator> createGapConductanceOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                       boost::shared_ptr<AMP::InputDatabase> input_db,
                                       boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel);

  static boost::shared_ptr<Operator> createLinearMechanicsOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                    boost::shared_ptr<AMP::InputDatabase> input_db,
                                    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel);

  static boost::shared_ptr<Operator> createNonlinearMechanicsOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                       boost::shared_ptr<AMP::InputDatabase> input_db,
                                       boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel);

  static boost::shared_ptr<Operator>
  createLinearBVPOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
			   std::string operatorName,
			   boost::shared_ptr<AMP::InputDatabase> input_db,
			   boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel,
			   boost::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> localModelFactory );
  
  static boost::shared_ptr<Operator>
  createNonlinearBVPOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
			      std::string operatorName,
			      boost::shared_ptr<AMP::InputDatabase> input_db,
			      boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel,
			      boost::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> localModelFactory);
  
  
  static boost::shared_ptr<BoundaryOperator>
  createDirichletMatrixCorrection( AMP::Mesh::Mesh::shared_ptr meshAdapter,
				   boost::shared_ptr<AMP::InputDatabase> input_db,
				   AMP::Operator::Operator::shared_ptr volumeOperator,
				   boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );

  static boost::shared_ptr<BoundaryOperator>
  createMassMatrixCorrection( AMP::Mesh::Mesh::shared_ptr meshAdapter,
			      boost::shared_ptr<AMP::InputDatabase> input_db,
			      AMP::Operator::Operator::shared_ptr volumeOperator,
			      boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel);

  static boost::shared_ptr<BoundaryOperator>
  createRobinMatrixCorrection( AMP::Mesh::Mesh::shared_ptr meshAdapter,
			       boost::shared_ptr<AMP::InputDatabase> input_db,
			       AMP::Operator::Operator::shared_ptr volumeOperator,
			       boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );

  static boost::shared_ptr<BoundaryOperator>
  createRobinVectorCorrection( AMP::Mesh::Mesh::shared_ptr meshAdapter,
			       boost::shared_ptr<AMP::InputDatabase> input_db,
			       AMP::Operator::Operator::shared_ptr volumeOperator,
			       boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );
  
  static boost::shared_ptr<BoundaryOperator>
  createNeumannVectorCorrection( AMP::Mesh::Mesh::shared_ptr meshAdapter,
				 boost::shared_ptr<AMP::InputDatabase> input_db,
				 AMP::Operator::Operator::shared_ptr volumeOperator,
				 boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );

  static boost::shared_ptr<BoundaryOperator>
  createPressureBoundaryVectorCorrection(AMP::Mesh::Mesh::shared_ptr meshAdapter,
					 boost::shared_ptr<AMP::InputDatabase> input_db,
					 AMP::Operator::Operator::shared_ptr volumeOperator,
					 boost::shared_ptr<AMP::Operator::ElementPhysicsModel> & elementPhysicsModel);

  static boost::shared_ptr<BoundaryOperator>
  createPressureBoundaryVectorCorrection(AMP::Mesh::Mesh::shared_ptr meshAdapter,
					 boost::shared_ptr<AMP::InputDatabase> input_db,
					 boost::shared_ptr<AMP::Operator::ElementPhysicsModel> & elementPhysicsModel);

  static boost::shared_ptr<BoundaryOperator>
  createDirichletVectorCorrection( AMP::Mesh::Mesh::shared_ptr meshAdapter,
				   boost::shared_ptr<AMP::InputDatabase> input_db,
				   AMP::Operator::Operator::shared_ptr volumeOperator,
				   boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );
  
  static boost::shared_ptr<BoundaryOperator>
  createDirichletVectorCorrection(AMP::Mesh::Mesh::shared_ptr meshAdapter,
				  boost::shared_ptr<AMP::InputDatabase> input_db,
				  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> &);


};
 
}
}

#endif

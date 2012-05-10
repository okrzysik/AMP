
#ifndef included_AMP_LinearBVPOperator
#define included_AMP_LinearBVPOperator

#include "LinearOperator.h"
#include "BVPOperatorParameters.h"

namespace AMP {
namespace Operator {

  /**
   * Class LinearBVPOperator is meant to wrap a pointer to a linear volume or interior spatial operator
   * with a pointer to a BoundaryOperator that handles spatial surface boundary conditions. The constructor
   * takes a pointer to a BVPOperatorParameters object.
   */
  class LinearBVPOperator : public LinearOperator 
  {
    public:

      /**
       * Main constructor
         @param [in] parameters The parameters object contains a database object which must contain the
         following fields in addition to the fields expected by the base Operator class:

	 1. name: VolumeOperator, type: string, (required), name of the database associated with the volume operator

	 2. name: BoundaryOperator, type: string, (required), name of the database associated with the boundary operator

	 3. name: name, type: string, (required), must be set to LinearBVPOperator
	 
	 4. name: useSameLocalModelForVolumeAndBoundaryOperators, type: bool, (optional), default value: FALSE, when set to
	    to TRUE the same local model is used for both the volume and boundary operators
        */
      LinearBVPOperator(const boost::shared_ptr<BVPOperatorParameters>& parameters);

      /**
       * virtual destructor which does nothing
       */
      virtual ~LinearBVPOperator() { }

      /**
       * This function is useful for re-initializing/updating an operator
       */
      void reset(const boost::shared_ptr<OperatorParameters>& );

      AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() 
      {
        return d_volumeOperator->getInputVariable();
      }

      AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() 
      {
        return d_volumeOperator->getOutputVariable();
      }

      boost::shared_ptr< LinearOperator > getVolumeOperator(){ return d_volumeOperator; }

      boost::shared_ptr< BoundaryOperator > getBoundaryOperator(){ return d_boundaryOperator; }

      void modifyRHSvector(AMP::LinearAlgebra::Vector::shared_ptr rhs);

    protected :

      /**
       * shared pointer to a nonlinear volume or interior spatial operator
       */
      boost::shared_ptr< LinearOperator > d_volumeOperator;

      /**
       *  shared pointer to a boundary or surface operator that is responsible
       for apply operations on the boundary of the domain
       */
      boost::shared_ptr< BoundaryOperator > d_boundaryOperator;

    private :

  };

}
}

#endif



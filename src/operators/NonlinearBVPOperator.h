
#ifndef included_AMP_NonlinearBVPOperator
#define included_AMP_NonlinearBVPOperator

#include "operators/Operator.h"
#include "operators/boundary/BoundaryOperator.h"
#include "BVPOperatorParameters.h"

namespace AMP {
  namespace Operator {

    /**
     * Class NonlinearBVPOperator is meant to wrap a pointer to a nonlinear volume or interior spatial operator
     * with a pointer to a BoundaryOperator that handles spatial surface boundary conditions. The constructor
     * takes a pointer to a BVPOperatorParameters object
     */
    class NonlinearBVPOperator : public Operator 
    {
      public:

        /**
         * Main constructor
         @param [in] parameters The parameters object contains a database object which must contain the
         following fields in addition to the fields expected by the base Operator class:

         1. name: VolumeOperator, type: string, (required), name of the database associated with the volume operator

         2. name: BoundaryOperator, type: string, (required), name of the database associated with the boundary operator

         3. name: name, type: string, (required), must be either NonlinearBVPOperator

         4. name: useSameLocalModelForVolumeAndBoundaryOperators, type: bool, (optional), default value: FALSE, when set to
         to TRUE the same local model is used for both the volume and boundary operators
         */
        NonlinearBVPOperator(const boost::shared_ptr<BVPOperatorParameters>& parameters);

        /**
         * virtual destructor which does nothing
         */
        virtual ~NonlinearBVPOperator() { }

        /**
          The apply function for this operator performs the following operation:
          r = b*f+a*A(u), if f is not NULL and
          r = a*A(u), if f is NULL
          Here, A represents the action of the composite volume and boundary operator 
          */
        void apply(AMP::LinearAlgebra::Vector::const_shared_ptr, AMP::LinearAlgebra::Vector::const_shared_ptr, 
            AMP::LinearAlgebra::Vector::shared_ptr, double a = -1.0, double b=1.0);

        /**
         * This function is useful for re-initializing/updating an operator
         * \param params
         *        parameter object containing parameters to change
         */
        void reset(const boost::shared_ptr<OperatorParameters>& params);

        /**
         * Returns a shared pointer to a OperatorParameters object that can be used to construct
         * or update an Operator which corresponds to some approximation to the Jacobian or tangent matrix
         * of the NonlinearBVPOperator
         * 
         * \param x
         *        This function takes as input a current vector x at which the Jacobian or tangent matrix needs to be constructed
         */
        boost::shared_ptr<OperatorParameters> getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& x);

        boost::shared_ptr< Operator > getVolumeOperator(){ return d_volumeOperator; }

        boost::shared_ptr< BoundaryOperator > getBoundaryOperator(){ return d_boundaryOperator; }

        void modifyRHSvector(AMP::LinearAlgebra::Vector::shared_ptr rhs);

        void modifyInitialSolutionVector(AMP::LinearAlgebra::Vector::shared_ptr sol);

        AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() 
        {
          return d_volumeOperator->getInputVariable();
        }

        AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() 
        {
          return d_volumeOperator->getOutputVariable();
        }

        bool isValidInput(AMP::LinearAlgebra::Vector::shared_ptr &sol){ return d_volumeOperator->isValidInput(sol); }

      protected :

        /**
         * shared pointer to a nonlinear volume or interior spatial operator
         */
        boost::shared_ptr< Operator > d_volumeOperator;

        /**
         *  shared pointer to a boundary or surface operator that is responsible for apply operations on the boundary of the domain
         */
        boost::shared_ptr< BoundaryOperator > d_boundaryOperator;

      private :

    };

  }
}

#endif



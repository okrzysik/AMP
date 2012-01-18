
#ifndef included_AMP_Operator
#define included_AMP_Operator

#include "boost/shared_ptr.hpp"

#include "operators/OperatorParameters.h"

#include "vectors/Vector.h"

#include "vectors/Variable.h"

#include <string>

#ifdef DEBUG_CHECK_ASSERTIONS
#include <cassert>
#endif


namespace AMP {
  namespace Operator {

    /**
     * Class Operator is an abstract base class for representing
     * a discrete operator which may be linear or nonlinear.  Concrete
     * implementations must include an implementation of the apply() function. 
     * The constructor for the class takes a pointer to
     * a OperatorParameters object. 
     */

    class Operator {

      public :

        typedef boost::shared_ptr <Operator>  shared_ptr;

        /**
          Default constructor
          */
        Operator(void);

        /**
          Constructor
          */
        Operator(const boost::shared_ptr<OperatorParameters> & params);

        /**
          Destructor
          */
        virtual ~Operator() { }

        /**
         * This function is useful for re-initializing/updating an operator
         * \param params
         *        parameter object containing parameters to change
         */
        virtual void reset(const boost::shared_ptr<OperatorParameters>& params);

        /**
          This base class can not give a meaningful definition of apply. See the derived classes for
          how they define apply. Each operator is free to define apply in a way that is appropriate
          for that operator.
          */
        virtual void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, 
            const  AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr  &r,
            const double a = -1.0, const double b = 1.0) = 0;

        /**
         * This function returns a OperatorParameters object
         * constructed by the operator which contains parameters from
         * which the Jacobian or portions of the Jacobian required by
         * solvers and preconditioners can be constructed. Returning
         * a parameter object instead of the Jacobian itself is meant 
         * to give users more flexibility.
         */
        virtual boost::shared_ptr<OperatorParameters> 
          getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& ) {
            //Implemented in derived class.
            boost::shared_ptr<OperatorParameters> emptyPointer;
            return emptyPointer;
          }

        /**
         * Specify level of diagnostic information printed during iterations.
         * \param print_level
         *        zero prints none or minimal information, higher numbers provide increasingly
         *        verbose debugging information.
         */
        virtual void setDebugPrintInfoLevel(int print_level) {
          d_iDebugPrintInfoLevel = print_level;
        }

        virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
          //Implemented in derived classes
          AMP::LinearAlgebra::Variable::shared_ptr emptyPointer;
          return emptyPointer;
        }

        virtual AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() {
          //Implemented in derived classes
          AMP::LinearAlgebra::Variable::shared_ptr emptyPointer;
          return emptyPointer;
        }

        virtual bool isValidInput(boost::shared_ptr<AMP::LinearAlgebra::Vector>&){return true; }

        AMP::Mesh::Mesh::shared_ptr getMesh() {
          return d_Mesh;
        }

      protected :

        void getFromInput(const boost::shared_ptr<AMP::Database>& db);

        int d_iDebugPrintInfoLevel;

        int d_iObject_id;

        static int d_iInstance_id;

        AMP::Mesh::Mesh::shared_ptr d_Mesh;

      private :

    };

  }
}

#endif



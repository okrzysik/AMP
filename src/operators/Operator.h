#ifndef included_AMP_Operator
#define included_AMP_Operator

#include "utils/shared_ptr.h"

#include "operators/OperatorParameters.h"

#include "vectors/Vector.h"
#include "vectors/Variable.h"
#include "vectors/VectorSelector.h"

#include <string>



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

    typedef AMP::shared_ptr<AMP::Operator::Operator>  shared_ptr;

    //! Default constructor
    Operator(void);

    //! Constructor
    Operator(const AMP::shared_ptr<OperatorParameters> & params);

    //! Destructor
    virtual ~Operator() { }

    /**
      * This function is useful for re-initializing/updating an operator
      * \param params
      *    parameter object containing parameters to change
      */
     virtual void reset(const AMP::shared_ptr<OperatorParameters>& params);

    /**
      This base class can not give a meaningful definition of apply. See the derived classes for
      how they define apply. Each operator is free to define apply in a way that is appropriate
      for that operator.
      \param u: shared pointer to const input vector u
      \param f: shared pointer to output vector storing result of applying this operator
      */
    virtual void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u, 
			AMP::LinearAlgebra::Vector::shared_ptr f) = 0;

    /**
     * Default base class implementation of the residual: f-L(u)
     * \param f: shared pointer to const vector rhs
     * \param u: shared pointer to const vector u
     * \param r: shared pointer to vector residual
     */
    virtual void residual(AMP::LinearAlgebra::Vector::const_shared_ptr f, 
			  AMP::LinearAlgebra::Vector::const_shared_ptr u, 
			  AMP::LinearAlgebra::Vector::shared_ptr r);
    
    /**
     * This function returns a OperatorParameters object
     * constructed by the operator which contains parameters from
     * which new operators can be created. Returning
     * a parameter object instead of an Operator itself is meant 
     * to give users more flexibility. Examples of how this functionality
     * might be used would be the construction of Jacobian, frozen Jacobian,
     * preconditioner approximations to the Jacobian, adjoint operators etc
     * \param type: std:string specifying type of return operator parameters
     * being requested. Currently the valid option is Jacobian
     * \param u: const pointer to current solution vector 
     * \param params: pointer to additional parameters that might be required
     * to construct the return parameters
     */
    virtual AMP::shared_ptr<OperatorParameters> 
       getParameters(const std::string &type,
                     AMP::LinearAlgebra::Vector::const_shared_ptr u,
                     AMP::shared_ptr<OperatorParameters> params = NULL ) {

       NULL_USE(params);

       AMP::shared_ptr<OperatorParameters> rPointer;

       if(type=="Jacobian") {
          rPointer = getJacobianParameters(u);
       } else {
          AMP_ERROR("Unknown OperatorParameters type specified");
          // should be implemented in derived class
       }
       return rPointer;
      }

    /**
     * Specify level of diagnostic information printed during iterations.
     * \param print_level
     *    zero prints none or minimal information, higher numbers provide increasingly
     *    verbose debugging information.
     */
    virtual void setDebugPrintInfoLevel(int print_level) {
        d_iDebugPrintInfoLevel = print_level;
    }

    //! Return the output variable
    virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
        //Implemented in derived classes
        AMP::LinearAlgebra::Variable::shared_ptr emptyPointer;
        return emptyPointer;
    }

    //! Return the input variable
    virtual AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() {
        //Implemented in derived classes
        AMP::LinearAlgebra::Variable::shared_ptr emptyPointer;
        return emptyPointer;
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr subsetOutputVector(AMP::LinearAlgebra::Vector::shared_ptr vec);
    virtual AMP::LinearAlgebra::Vector::const_shared_ptr subsetOutputVector(AMP::LinearAlgebra::Vector::const_shared_ptr vec);

    virtual AMP::LinearAlgebra::Vector::shared_ptr subsetInputVector(AMP::LinearAlgebra::Vector::shared_ptr vec);
    virtual AMP::LinearAlgebra::Vector::const_shared_ptr subsetInputVector(AMP::LinearAlgebra::Vector::const_shared_ptr vec);

    virtual bool isValidInput(AMP::shared_ptr<AMP::LinearAlgebra::Vector>&){return true; }

    AMP::Mesh::Mesh::shared_ptr getMesh() {
        return d_Mesh;
    }


protected :

    void getFromInput(const AMP::shared_ptr<AMP::Database>& db);

    /**
     * This function returns a OperatorParameters object
     * constructed by the operator which contains parameters from
     * which new Jacobian operators can be created. Returning
     * a parameter object instead of an Operator itself is meant 
     * to give users more flexibility.
     */
    virtual AMP::shared_ptr<OperatorParameters> 
       getJacobianParameters( AMP::LinearAlgebra::Vector::const_shared_ptr u ) {
       NULL_USE(u);
         //Implemented in derived class.
         AMP::shared_ptr<OperatorParameters> emptyPointer;
         return emptyPointer;
      }

    int d_iDebugPrintInfoLevel;

    int d_iObject_id;

    static int d_iInstance_id;

    AMP::Mesh::Mesh::shared_ptr d_Mesh;

private :

};


}
}

#endif



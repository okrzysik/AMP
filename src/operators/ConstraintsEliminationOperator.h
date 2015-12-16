
#ifndef included_AMP_ConstraintsEliminationOperator
#define included_AMP_ConstraintsEliminationOperator

#include <operators/Operator.h>
#include <vector>
#include <vectors/Variable.h>
#include <vectors/Vector.h>

namespace AMP {
namespace Operator {

/**
  u^s = C u^m + d
  */
class ConstraintsEliminationOperator : public Operator
{

public:
    /**
      Constructor.
      @param [in] params
      */
    ConstraintsEliminationOperator( const AMP::shared_ptr<OperatorParameters> &params )
        : Operator( params )
    {
        AMP_INSIST( params->d_db->keyExists( "InputVariable" ), "key not found" );
        std::string inpVarName = params->d_db->getString( "InputVariable" );
        d_InputVariable.reset( new AMP::LinearAlgebra::Variable( inpVarName ) );

        AMP_INSIST( params->d_db->keyExists( "OutputVariable" ), "key not found" );
        std::string outVarName = params->d_db->getString( "OutputVariable" );
        d_OutputVariable.reset( new AMP::LinearAlgebra::Variable( outVarName ) );
    }

    /**
      Destructor
      */
    virtual ~ConstraintsEliminationOperator() {}

    /**
     * This function is useful for re-initializing/updating an operator
     * \param params
     *        parameter object containing parameters to change
     */
    virtual void reset( const AMP::shared_ptr<OperatorParameters> &params );

    /**
      Calls first addSlaveToMaster(...) and second setSlaveToZero(...) on the residual vector:
      r^m = r^m + C^T r^s
      r^s = 0

      The apply function for this operator, A, performs the following operation:
      r = A(u)
      Here, A(u) is simply a Matrix-Vector multiplication.
      @param [in] u input vector.
      @param [out] f residual/output vector.
      */
    virtual void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                        AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    /**
      @return The variable for the input vector.
      */
    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable();

    /**
      @return The variable for the output vector.
      */
    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable();

    /**
      u^m = u^m + C^T u^s
      */
    virtual void addSlaveToMaster( AMP::LinearAlgebra::Vector::shared_ptr u ) = 0;
    /**
      u^s = C u^m
      */
    virtual void copyMasterToSlave( AMP::LinearAlgebra::Vector::shared_ptr u ) = 0;
    /**
      u^s = 0
      */
    virtual void setSlaveToZero( AMP::LinearAlgebra::Vector::shared_ptr u );
    /**
      u^s = u^s + d
      */
    virtual void addShiftToSlave( AMP::LinearAlgebra::Vector::shared_ptr u );

protected:
    std::vector<size_t> d_SlaveIndices;
    std::vector<double> d_SlaveShift;

private:
    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_InputVariable;  /**< Input variable */
    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_OutputVariable; /**< Output variable */
};
}
}

#endif

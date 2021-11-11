
#ifndef included_AMP_FlowFrapconOperator
#define included_AMP_FlowFrapconOperator

#include "AMP/operators/Operator.h"
#include "AMP/operators/subchannel/FlowFrapconOperatorParameters.h"


namespace AMP {
namespace Operator {

/**
  A class to represent Frapcon Flow operator. This operator acts as a heat sink and
  should be used to compute the 1D flow temperature in the axial direction of the pellet/clad.
  */
class FlowFrapconOperator : public Operator
{
public:
    /**
      Constructor creates a simpleVariables for Input and Output. The reset is called to
      read the flow parameters.
      */
    explicit FlowFrapconOperator( std::shared_ptr<const FlowFrapconOperatorParameters> params );

    /**
      Destructor
      */
    virtual ~FlowFrapconOperator() {}

    std::string type() const override { return "FlowFrapconOperator"; }

    /**
      For this operator we have an in-place apply.
      @param [in]  u input vector.
      @param [out] f output vector.
      */
    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    /**
      This function reads the entries of the database for the flow operator
      and can also be used to change the parameters if required.
      */
    void reset( std::shared_ptr<const OperatorParameters> params ) override;

    /*
          static bool sort_nodes_in_z ( const libMesh::Node *one , const libMesh::Node *two ) {
            return (*one)(2) < (*two)(2);
          }
    */
    std::shared_ptr<AMP::LinearAlgebra::Variable> createInputVariable( const std::string &name,
                                                                       int varId = -1 )
    {
        (void) varId;
        return d_inpVariable->cloneVariable( name );
    }

    std::shared_ptr<AMP::LinearAlgebra::Variable> createOutputVariable( const std::string &name,
                                                                        int varId = -1 )
    {
        (void) varId;
        return d_outVariable->cloneVariable( name );
    }

    std::shared_ptr<AMP::LinearAlgebra::Variable> getInputVariable() override
    {
        return d_inpVariable;
    }

    std::shared_ptr<AMP::LinearAlgebra::Variable> getOutputVariable() override
    {
        return d_outVariable;
    }

    virtual AMP::LinearAlgebra::Vector::shared_ptr
    subsetOutputVector( AMP::LinearAlgebra::Vector::shared_ptr vec ) override;

    virtual AMP::LinearAlgebra::Vector::shared_ptr
    subsetInputVector( AMP::LinearAlgebra::Vector::shared_ptr vec ) override;

    virtual AMP::LinearAlgebra::Vector::const_shared_ptr
    subsetOutputVector( AMP::LinearAlgebra::Vector::const_shared_ptr vec ) override;

    virtual AMP::LinearAlgebra::Vector::const_shared_ptr
    subsetInputVector( AMP::LinearAlgebra::Vector::const_shared_ptr vec ) override;


    /**
      @param [in] zloc is the location vector in z direction.
      */
    void setZLocations( const std::vector<double> &zloc ) { zPoints = zloc; }

    void setVector( AMP::LinearAlgebra::Vector::shared_ptr frozenVec ) { d_cladVec = frozenVec; }

    AMP::LinearAlgebra::Vector::shared_ptr getVector() { return d_cladVec; }

    /**
      This member function returns the 1D locations stl vector.
      */
    std::vector<double> getZLocations() { return zPoints; }

    double getHeatCapacity( double T_b )
    {
        double cp = 0;

        if ( T_b < 544 ) {

            cp = 2.4e5;
        } else if ( 544 <= T_b && T_b < 588 ) {

            cp = 2.4e5 * ( 1.0 + 2.9e-3 * ( 1.8 * T_b - 1031.0 ) );
        } else {

            cp = 2.4e5 * ( 1.0 + 2.9e-3 * ( 1.8 * T_b - 979.4 ) );
        }

        return cp;
    }


    short int d_boundaryId;

    int d_numpoints; /**< Number of points in z direction */

    std::vector<unsigned int> d_dofIds;

    std::vector<double> zPoints; /**< vector to hold z locations */

    AMP::LinearAlgebra::Vector::shared_ptr d_cladVec;

    double d_De; /**< Channel Diameter */

    double Cp; /**< Heat Capacity of Coolant */

    double d_G; /**< Coolant Mass Flux */

    double d_Tin; /**< Coolant Temp Tin */

    double d_K; /**< Coolant conductivity */

    double d_Re; /**< Reynolds Number */

    double d_Pr; /**< Prandtl Number */

    /* Since the map has been taken out the Flow operator
       now expects a SimpleVariable for input & output */
    std::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable; /**< Simple Input Variable */

    std::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable; /**< Simple Output Variable */

protected:
    /**
      This function returns a parameter object that can be used to reset the corresponding
      FlowFrapconOperator operator.
      */
    std::shared_ptr<OperatorParameters>
    getJacobianParameters( AMP::LinearAlgebra::Vector::const_shared_ptr u ) override;

private:
};
} // namespace Operator
} // namespace AMP

#endif

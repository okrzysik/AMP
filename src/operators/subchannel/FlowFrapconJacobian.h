
#ifndef included_AMP_FlowFrapconJacobian
#define included_AMP_FlowFrapconJacobian

#include "operators/Operator.h"
#include "operators/subchannel/FlowFrapconJacobianParameters.h"
#include "vectors/SimpleVector.h"
//#include "operators/map/Map3Dto1D.h"
//#include "operators/map/Map1Dto3D.h"

/* Libmesh files */
#include "libmesh/fe_base.h"
#include "libmesh/fe_type.h"
#include "libmesh/quadrature.h"

namespace AMP {
namespace Operator {

/**
  A class to represent Frapcon Flow operator. This operator acts as a heat sink and
  should be used to compute the 1D flow temperature in the axial direction of the pellet/clad.
  */
class FlowFrapconJacobian : public Operator
{
public:
    /**
      Constructor creates a simpleVariables for Input and Output. The reset is called to
      read the flow parameters.
      */
    explicit FlowFrapconJacobian( const AMP::shared_ptr<FlowFrapconJacobianParameters> &params );

    /**
      Destructor
      */
    virtual ~FlowFrapconJacobian() {}

    /**
      For this operator we have an in-place apply.
      @param [in]  u input vector.
      @param [out] r residual/output vector.
      */
    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr r ) override;

    /**
      This function reads the entries of the database for the flow operator
      and can also be used to change the parameters if required.
      */
    void reset( const AMP::shared_ptr<OperatorParameters> &params ) override;

    /*
          static bool sort_nodes_in_z ( const ::Node *one , const ::Node *two ) {
            return (*one)(2) < (*two)(2);
          }
    */
    AMP::LinearAlgebra::Variable::shared_ptr createInputVariable( const std::string &name,
                                                                  int varId = -1 );

    AMP::LinearAlgebra::Variable::shared_ptr createOutputVariable( const std::string &name,
                                                                   int varId = -1 );

    virtual AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override;

    virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override;

    void setInputVariableName( const std::string &name, int varId = -1 );

    void setOutputVariableName( const std::string &name, int varId = -1 );

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
    void setZLocations( std::vector<double> zloc ) { zPoints = zloc; }

    std::vector<double> getZLocations() { return zPoints; }

    double getCp() { return Cp; }

    double getdCp() { return dCp; }

    void setVector( AMP::LinearAlgebra::Vector::shared_ptr frozenVec ) { d_cladVec = frozenVec; }

    void setFrozenVector( AMP::LinearAlgebra::Vector::shared_ptr frozenVec )
    {
        d_frozenVec = frozenVec;
    }

    double getHeatCapacity( double T_b )
    {
        double cp = 0.0;
        if ( T_b < 544 ) {
            cp = 2.4e5;
        } else if ( 544 <= T_b && T_b < 588 ) {
            cp = 2.4e5 * ( 1.0 + 2.9e-3 * ( 1.8 * T_b - 1031.0 ) );
        } else {
            cp = 2.4e5 * ( 1.0 + 2.9e-3 * ( 1.8 * T_b - 979.4 ) );
        }
        return cp;
    }

    double getHeatCapacityGradient( double T_b )
    {
        double dcp = 0.0;
        if ( T_b >= 544 ) {
            dcp = 2.4e5 * ( 2.9e-3 * 1.8 * T_b );
        }
        return dcp;
    }

    /*
       AMP::shared_ptr< std::vector<double>  > getFlowSolution() {
       return flowSolutionVector;
       }
       */

    AMP::LinearAlgebra::Vector::shared_ptr d_cladVec;
    AMP::LinearAlgebra::Vector::shared_ptr d_frozenVec;

    AMP::LinearAlgebra::Vector::shared_ptr d_localCladVec;
    AMP::LinearAlgebra::Vector::shared_ptr d_localFrozenVec;

    int d_numpoints; /**< Number of points in z direction */

    std::vector<unsigned int> d_dofIds;

    std::vector<double> zPoints; /**< vector to hold z locations */

    double d_De; /**< Channel Diameter */

    double Cp; /**< Heat Capacity of Coolant */

    double dCp; /**< Gradient of Heat Capacity */

    double d_G; /**< Coolant Mass Flux */

    double d_Tin; /**< Coolant Temp Tin */

    double d_K; /**< Coolant conductivity */

    double d_Re; /**< Reynolds Number */

    double d_Pr; /**< Prandtl Number */

    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable; /**< Input Variable */

    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable; /**< Output Variable */

    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_SimpleVariable; /**< Simple Input Variable */

    //      AMP::shared_ptr<AMP::Operator::Map3Dto1D> d_Map3to1;
    AMP::LinearAlgebra::Vector::shared_ptr flowInput;

    //      AMP::shared_ptr<AMP::Operator::Map1Dto3D> d_Map1to3;
    AMP::LinearAlgebra::Vector::shared_ptr flowOutput;

protected:
private:
};
} // namespace Operator
} // namespace AMP

#endif

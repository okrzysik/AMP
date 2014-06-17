
#ifndef included_AMP_FlowFrapconOperator
#define included_AMP_FlowFrapconOperator

#include "operators/Operator.h"
#include "vectors/SimpleVector.h"
#include "operators/subchannel/FlowFrapconOperatorParameters.h"

/* Libmesh files */
#include "libmesh/fe_type.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature.h"

namespace AMP {
namespace Operator {

  /**
    A class to represent Frapcon Flow operator. This operator acts as a heat sink and 
    should be used to compute the 1D flow temperature in the axial direction of the pellet/clad. 
    */
  class FlowFrapconOperator : public Operator
  {
    public :

      /**
        Constructor creates a simpleVariables for Input and Output. The reset is called to 
        read the flow parameters.
        */
      FlowFrapconOperator(const boost::shared_ptr<FlowFrapconOperatorParameters> & params)
        : Operator (params)
      {
        std::string inpVar = params->d_db->getString("InputVariable");
        d_inpVariable.reset(new AMP::LinearAlgebra::Variable(inpVar));

        std::string outVar = params->d_db->getString("OutputVariable");
        d_outVariable.reset(new AMP::LinearAlgebra::Variable(outVar));

        reset(params);
      }

      /**
        Destructor
        */
      virtual ~FlowFrapconOperator() { }

      /**
        For this operator we have an in-place apply.
        @param [in]  f auxillary/rhs vector. 
        @param [in]  u input vector. 
        @param [out] r residual/output vector. 
        @param [in]  a first constant used in the expression: r = a*A(u) + b*f. The default value is -1.
        @param [in]  b second constant used in the expression: r = a*A(u) + b*f. The default value is 1.
        */
      void apply(AMP::LinearAlgebra::Vector::const_shared_ptr f, AMP::LinearAlgebra::Vector::const_shared_ptr u,
          AMP::LinearAlgebra::Vector::shared_ptr r, const double a = -1.0, const double b = 1.0);

      /**
        This function reads the entries of the database for the flow operator
        and can also be used to change the parameters if required.
        */
      void reset(const boost::shared_ptr<OperatorParameters>& params);

/*
      static bool sort_nodes_in_z ( const ::Node *one , const ::Node *two ) {
        return (*one)(2) < (*two)(2);
      }
*/
      AMP::LinearAlgebra::Variable::shared_ptr createInputVariable (const std::string & name, int varId = -1)
      {
        (void) varId;      
        return d_inpVariable->cloneVariable(name);
      }

      AMP::LinearAlgebra::Variable::shared_ptr createOutputVariable (const std::string & name, int varId = -1) 
      {
        (void) varId;      
        return d_outVariable->cloneVariable(name);
      }

      AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() {
        return d_inpVariable;
      }

      AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
        return d_outVariable;
      }

      virtual AMP::LinearAlgebra::Vector::shared_ptr subsetOutputVector(AMP::LinearAlgebra::Vector::shared_ptr vec);

      virtual AMP::LinearAlgebra::Vector::shared_ptr subsetInputVector(AMP::LinearAlgebra::Vector::shared_ptr vec);

      virtual AMP::LinearAlgebra::Vector::const_shared_ptr subsetOutputVector(AMP::LinearAlgebra::Vector::const_shared_ptr vec);

      virtual AMP::LinearAlgebra::Vector::const_shared_ptr subsetInputVector(AMP::LinearAlgebra::Vector::const_shared_ptr vec);


      /**
        @param [in] zloc is the location vector in z direction.
        */
      void setZLocations(std::vector<double> zloc)
      {
        zPoints = zloc;
      }

      void setVector(AMP::LinearAlgebra::Vector::shared_ptr &frozenVec) {
        d_cladVec = frozenVec;
      }
      
      AMP::LinearAlgebra::Vector::shared_ptr getVector() {
        return d_cladVec ;
      }
      
      /**
        This member function returns the 1D locations stl vector.
        */
      std::vector<double> getZLocations()
      {
        return zPoints ;
      }

      double getHeatCapacity(double T_b)
      {
        double cp;

        if(T_b < 544){

          cp = 2.4e5;

        }else if (544 <= T_b && T_b < 588){

          cp = 2.4e5 * (1.0  + 2.9e-3 * (1.8*T_b - 1031.0));

        }else if (T_b >= 588){

          cp = 2.4e5 * (1.0  + 2.9e-3 * (1.8*T_b - 979.4));

        }

        return cp;
      }

      /**
        This function returns a parameter object that can be used to reset the corresponding
        FlowFrapconOperator operator.
        */
      boost::shared_ptr<OperatorParameters> getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& );


      short int d_boundaryId;

      int d_numpoints; /**< Number of points in z direction */

      std::vector<unsigned int> d_dofIds;

      std::vector<double> zPoints; /**< vector to hold z locations */

      AMP::LinearAlgebra::Vector::shared_ptr d_cladVec;

      double d_De; /**< Channel Diameter */

      double Cp; /**< Heat Capacity of Coolant */

      double d_G;  /**< Coolant Mass Flux */

      double d_Tin;/**< Coolant Temp Tin */

      double d_K;  /**< Coolant conductivity */

      double d_Re; /**< Reynolds Number */

      double d_Pr; /**< Prandtl Number */

      /* Since the map has been taken out the Flow operator 
         now expects a SimpleVariable for input & output */
      boost::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable; /**< Simple Input Variable */

      boost::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable; /**< Simple Output Variable */

    private :

  };

}
}

#endif


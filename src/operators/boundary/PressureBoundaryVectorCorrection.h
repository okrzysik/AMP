
#ifndef included_AMP_PressureBoundaryVectorCorrection
#define included_AMP_PressureBoundaryVectorCorrection

#include "BoundaryOperator.h"
#include "PressureBoundaryVectorCorrectionParameters.h"

/* Libmesh files */
#include "fe_type.h"
#include "fe_base.h"
#include "elem.h"
#include "quadrature.h"

#include "enum_order.h"
#include "enum_fe_family.h"
#include "enum_quadrature_type.h"
#include "auto_ptr.h"
#include "string_to_enum.h"

#include <string>
#include <vector>

namespace AMP {
namespace Operator {

  /**
    A class to impose Pressure Boundary Conditions for both Linear and Nonlinear mechanics operator.
    For both the Linear/Nonlinear operator to impose these conditions involves adding the corrections
    to the RHS vector at the appropriate locations. When you do not impose these Pressure condition for 
    the weak formulation, a natural condition is assumed.
    */
  class PressureBoundaryVectorCorrection : public BoundaryOperator
  {
    public :

      /**
        Constructor. This function reads all the parameters required for surface elements.
        */
      PressureBoundaryVectorCorrection(const boost::shared_ptr<PressureBoundaryVectorCorrectionParameters> & params)
        : BoundaryOperator (params)
      {
              d_params = params;

          std::string feTypeOrderName = (params->d_db)->getStringWithDefault("FE_ORDER", "FIRST");

          libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>(feTypeOrderName);

          std::string feFamilyName = (params->d_db)->getStringWithDefault("FE_FAMILY", "LAGRANGE");

          libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>(feFamilyName);

          std::string qruleTypeName = (params->d_db)->getStringWithDefault("QRULE_TYPE", "QGAUSS");

          libMeshEnums::QuadratureType qruleType = Utility::string_to_enum<libMeshEnums::QuadratureType>(qruleTypeName);

              const unsigned int dimension3 = 3;

          d_feType.reset( new ::FEType(feTypeOrder, feFamily) );

          d_fe.reset( (::FEBase::build((dimension3 - 1), (*d_feType))).release() );

          d_fe_3d.reset( (::FEBase::build(dimension3, (*d_feType))).release() );

          std::string qruleOrderName = (params->d_db)->getStringWithDefault("QRULE_ORDER", "DEFAULT");

          libMeshEnums::Order qruleOrder;

          if(qruleOrderName == "DEFAULT") {
              qruleOrder = d_feType->default_quadrature_order();
          } else {
              qruleOrder = Utility::string_to_enum<libMeshEnums::Order>(qruleOrderName);
          }

          d_qrule.reset( (::QBase::build(qruleType, (dimension3 - 1), qruleOrder)).release() );

          d_fe->attach_quadrature_rule( d_qrule.get() );

          d_fe_3d->attach_quadrature_rule( d_qrule.get() );

          d_JxW = &(d_fe->get_JxW());

          d_variable = params->d_variable;

          reset(params);
      }

      /**
        Set the variable for the vector that will be used with this operator.
        */
      void setVariable(const AMP::LinearAlgebra::Variable::shared_ptr & var) {
        d_variable = var;
      }

      /**
        Destructor
        */
      virtual ~PressureBoundaryVectorCorrection() { }

      /**
        Sets pressure values into the appropriate locations of the output vector (r). 
        It also scales the total pressure by an amount "a". It is used for multiple 
        loading steps.
        */
      virtual void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
              AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a = -1.0, const double b = 1.0);

      /**
        This function computes the surface integral for either constant pressure values 
        across the boundary.
        */
      void computeRHScorrection();

      /**
        This function reads parameters related to boundary Ids
        */
      virtual void reset(const boost::shared_ptr<OperatorParameters>& params);

      /**
        Adds a vector to the RHS vector. 
        */
      void addRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr rhs) {
          AMP::LinearAlgebra::Vector::shared_ptr myRhs = rhs->subsetVectorForVariable(d_variable);
          myRhs->add(myRhs, d_rhsCorrectionAdd);
      }

      /**
        This function returns a parameter object that can be used to reset the corresponding
        PressureBoundaryVectorCorrection operator.
        */
      boost::shared_ptr<OperatorParameters> getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& );

      boost::shared_ptr<OperatorParameters> getParameters() {
        return d_params;
      }

      void setVariablePressure(const AMP::LinearAlgebra::Vector::shared_ptr &pressure){
        d_variablePressure = pressure->subsetVectorForVariable(d_variable);
      }

      std::vector<std::vector<std::vector<Point> > > getNormals();

    protected :

      std::vector<short int> d_boundaryIds;

      std::vector<double> d_pressureValues;

      AMP::LinearAlgebra::Vector::shared_ptr d_rhsCorrectionAdd;

      //This must be a simple variable not a dual or multivariable
      AMP::LinearAlgebra::Variable::shared_ptr d_variable;

      bool d_isConstantPressure;

      AMP::LinearAlgebra::Vector::shared_ptr d_Frozen;

      AMP::LinearAlgebra::Vector::shared_ptr d_variablePressure;

      int d_numIds;

      const std::vector<Real> *d_JxW;

      const std::vector<std::vector<Real> > *d_phi;

      const std::vector<std::vector<RealGradient> > *d_dphi;

      const std::vector<Point> *d_normal;

      boost::shared_ptr < ::FEType > d_feType;

      boost::shared_ptr < ::FEBase > d_fe;

      boost::shared_ptr < ::FEBase > d_fe_3d;

      boost::shared_ptr < ::QBase > d_qrule;

      const ::Elem *d_elem;

      const ::Elem *e_elem;

      boost::shared_ptr<PressureBoundaryVectorCorrectionParameters> d_params;

    private :

  };

}
}

#endif


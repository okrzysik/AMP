
#ifndef included_AMP_RobinMatrixCorrection
#define included_AMP_RobinMatrixCorrection

#include "operators/boundary/BoundaryOperator.h"
#include "operators/boundary/libmesh/NeumannVectorCorrection.h"
#include "operators/boundary/libmesh/NeumannVectorCorrectionParameters.h"
#include "operators/boundary/libmesh/RobinMatrixCorrectionParameters.h"

#include "discretization/createLibmeshElements.h"

/* Libmesh files */
#include "libmesh/fe_type.h"
#include "libmesh/fe_base.h"
#include "libmesh/elem.h"
#include "libmesh/quadrature.h"

#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/string_to_enum.h"

#include <string>

namespace AMP {
namespace Operator {

  /**
    A class to impose Robin Boundary conditions for a linear operator. Robin Condition
    is also known as mixed condition of Dirichlet and Neumann Flux conditions. This can 
    be written as \f$\alpha k*\frac{\partial u}{\partial n} + \beta h*u = \gamma*c \f$.
    Imposing this condition would involve:
    1) Imposing a Neumann Flux condition on the RHS Vector
    2) Make appropriate matrix corrections on the boundary nodes.
    */
  class RobinMatrixCorrection : public BoundaryOperator
  {
    public :

    /**
       Constructor. This function reads all the parameters required for surface elements.
       This also constructs new NeumannVectorCorrection parameters and calls it reset.
    */
    RobinMatrixCorrection(const boost::shared_ptr<RobinMatrixCorrectionParameters> & params);
    
    /**
       Set the variable for the vector that will used with this operator.
    */
    void setVariable(const AMP::LinearAlgebra::Variable::shared_ptr & var) {
      d_variable = var;
    }
    
    /**
       Destructor
    */
    virtual ~RobinMatrixCorrection() { }
    
    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr, AMP::LinearAlgebra::Vector::const_shared_ptr,
           AMP::LinearAlgebra::Vector::shared_ptr, const double = -1.0, const double = 1.0)
    {
      //Do Nothing
    }
    
    /**
       This function reads parameters related to boundary Ids. Since this class allow
       for variable flux values, the parameters stores a vector of values. This vector
       is passed to the NeumannVectorCorrection. This function also does a matrix 
       correction on the boundary.
    */
    void reset(const boost::shared_ptr<OperatorParameters>& params);
    
    /**
       Adds a Neumann Correction Vector to the RHS vector. 
    */
    void addRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr rhs)
    {
      d_NeumannCorrection->addRHScorrection(rhs);
    }
    
    
  protected :
    
    Discretization::createLibmeshElements libmeshElements;

    std::vector<short int> d_boundaryIds;
    
    std::vector<std::vector<unsigned int> >d_dofIds;
    
    std::vector<std::vector<double> > d_robinValues;
    
    AMP::LinearAlgebra::Variable::shared_ptr d_variable;
    
    double d_hef;  //Convective Coefficient
    
    double d_alpha;  // pre-factor solid flux
    
    double d_beta;
    
    double d_gamma;
    
    const std::vector<Real> *d_JxW;
    
    const std::vector<std::vector<Real> > *d_phi;
    
    boost::shared_ptr < ::FEType > d_feType;
    
    boost::shared_ptr < ::FEBase > d_fe;
    
    boost::shared_ptr < ::QBase > d_qrule;
    
    std::string         d_qruleOrderName; 

    libMeshEnums::Order d_feTypeOrder;
    libMeshEnums::FEFamily d_feFamily;
    libMeshEnums::QuadratureType d_qruleType;
    libMeshEnums::Order d_qruleOrder;

    AMP::LinearAlgebra::Vector::shared_ptr d_Frozen;
    
    boost::shared_ptr<RobinPhysicsModel> d_robinPhysicsModel;

    AMP::Discretization::DOFManager::shared_ptr d_dofManager; 

  private :

    boost::shared_ptr<NeumannVectorCorrection> d_NeumannCorrection; 
    boost::shared_ptr<NeumannVectorCorrectionParameters> d_NeumannParams; 

  };

}
}

#endif


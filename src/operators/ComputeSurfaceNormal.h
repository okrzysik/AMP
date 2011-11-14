
#ifndef included_AMP_ComputeSurfaceNormal
#define included_AMP_ComputeSurfaceNormal

#include "boost/shared_ptr.hpp"

#include "operators/OperatorParameters.h"
#include "operators/Operator.h"

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

  class ComputeSurfaceNormal 
  {
    public :

      ComputeSurfaceNormal(const boost::shared_ptr<OperatorParameters> & params);

      virtual ~ComputeSurfaceNormal() { }

      AMP::LinearAlgebra::Vector::shared_ptr getNormals(const AMP::LinearAlgebra::Vector::shared_ptr &u);

      void setVariable(const AMP::LinearAlgebra::Variable::shared_ptr & u);

    protected :

      std::vector<short int> d_boundaryIds;

      AMP::LinearAlgebra::Vector::shared_ptr d_inputVec;

      AMP::LinearAlgebra::Vector::shared_ptr d_outputVec;

      int d_numIds;
       
      int d_dimension;

      const std::vector<Real> *d_JxW;

      const std::vector<std::vector<Real> > *d_phi;

      const std::vector<std::vector<RealGradient> > *d_dphi;

      const std::vector<Point> *d_normal;

      boost::shared_ptr < ::FEType > d_feType;

      boost::shared_ptr < ::FEBase > d_fe_face;

      boost::shared_ptr < ::FEBase > d_fe;

      boost::shared_ptr < ::QBase > d_qrule;

      const ::Elem *d_face;

      const ::Elem *d_elem;

      AMP::Mesh::MeshManager::Adapter::shared_ptr d_MeshAdapter;

      AMP::LinearAlgebra::Variable::shared_ptr d_variable;

    private :

  };

}
}

#endif


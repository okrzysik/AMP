#ifndef included_AMP_ContactResidualCorrection
#define included_AMP_ContactResidualCorrection

#include "operators/Operator.h"
#include "vectors/MultiVariable.h"

namespace AMP {
namespace Operator {

  typedef OperatorParameters ContactResidualCorrectionParameters;

  class ContactResidualCorrection : public Operator {
    public:
      ContactResidualCorrection(const boost::shared_ptr<ContactResidualCorrectionParameters> & params)
        : Operator (params) {  }

      ~ContactResidualCorrection() { }

      void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
          AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a = -1.0, const double b = 1.0);

      void reset(const boost::shared_ptr<OperatorParameters>& params) {  }

      void setMasterVariable(const AMP::LinearAlgebra::Variable::shared_ptr & var) {
        d_masterVariable = var;
      }

      void setSlaveVariable(const AMP::LinearAlgebra::Variable::shared_ptr & var) {
        d_slaveVariable = var;
      }

      void setMasterMesh(const AMP::Mesh::Mesh::shared_ptr & mesh) {
        d_Mesh = mesh;
      }

      void setSlaveMesh(const AMP::Mesh::Mesh::shared_ptr & mesh) {
        d_slaveMesh = mesh;
      }

      void setMasterNodes(const std::vector<unsigned int> & vec) {
        d_masterNodes = vec;
      }

      void setSlaveNodes(const std::vector<unsigned int> & vec) {
        d_slaveNodes = vec;
      }

      void setDofs(const std::vector<std::vector<unsigned int> > & vec) {
        d_dofs = vec;
      }

      AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
        boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> retVariable( new AMP::LinearAlgebra::MultiVariable("ContactVariable"));
        retVariable->add(d_masterVariable);
        retVariable->add(d_slaveVariable);
        return retVariable;
      }

      AMP::LinearAlgebra::Variable::shared_ptr getInputVariable(int varId = -1) {
        boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> retVariable( new AMP::LinearAlgebra::MultiVariable("ContactVariable"));
        retVariable->add(d_masterVariable);
        retVariable->add(d_slaveVariable);
        return retVariable;
      }

      AMP::LinearAlgebra::Variable::shared_ptr getMasterVariable() {
        return d_masterVariable;
      }

      AMP::LinearAlgebra::Variable::shared_ptr getSlaveVariable() {
        return d_slaveVariable;
      }

      AMP::Mesh::Mesh::shared_ptr getMasterMesh() {
        return d_Mesh;
      }

      AMP::Mesh::Mesh::shared_ptr getSlaveMesh() {
        return d_slaveMesh;
      }

      std::vector<unsigned int> getMasterNodes() {
        return d_masterNodes;
      }

      std::vector<unsigned int> getSlaveNodes() {
        return d_slaveNodes;
      }

      std::vector<std::vector<unsigned int> > getDofs( ) {
        return d_dofs;
      }

    private:

      AMP::LinearAlgebra::Variable::shared_ptr d_masterVariable;
      AMP::LinearAlgebra::Variable::shared_ptr d_slaveVariable;

      std::vector<unsigned int> d_masterNodes;
      std::vector<unsigned int> d_slaveNodes;

      std::vector<std::vector<unsigned int> > d_dofs;

      AMP::Mesh::Mesh::shared_ptr d_slaveMesh;
  };

}
}

#endif


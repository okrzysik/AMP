
#ifndef included_AMP_OperatorParameters
#define included_AMP_OperatorParameters



#include "boost/shared_ptr.hpp"

#include "utils/Database.h"

#include "utils/ParameterBase.h"

#include "ampmesh/MeshManager.h"

namespace AMP {
namespace Operator {

  /**\class OperatorParameters
   * 
   * OperatorParameters encapsulates parameters used to initialize or reset
   * operators. It is an abstract base class.
   */

  class OperatorParameters: public ParameterBase
  {
    public :
      /**
       * Construct and initialize a parameter list according to input
       * data.  Guess what the required and optional keywords are.
       */
      OperatorParameters(const boost::shared_ptr<AMP::Database> & db)
        : d_db(db) {  }

      /**
       * Destructor.
       */
      virtual ~OperatorParameters() { }

      /**
       *  Database object which needs to be initialized specific to the solver.
       *  Documentation for parameters required by each solver can be found in the
       *  documentation for the solver.
       */
      boost::shared_ptr<AMP::Database> d_db;

      AMP::Mesh::MeshManager::Adapter::shared_ptr d_MeshAdapter;

    protected :

    private :

  };

}
}

#endif


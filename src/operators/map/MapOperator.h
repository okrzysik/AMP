#ifndef included_AMP_MapOperator
#define included_AMP_MapOperator


#include "utils/InputDatabase.h"

#include "boost/shared_ptr.hpp"

#include "operators/Operator.h"
#include "operators/OperatorParameters.h"
#include "operators/map/MapOperatorParameters.h"
#include "vectors/Vector.h"
#include "vectors/Variable.h"

#include <string>

#ifdef DEBUG_CHECK_ASSERTIONS
#include <cassert>
#endif


namespace AMP {
namespace Operator {

  /**
   * Class MapOperator is the base class for various mapping alogorithms. This
   * class stores a pointer to the mapAdapter to which the solution has to be 
   * mapped from Operator's meshAdapter. 
   */

  class MapOperator : public Operator
  {

    public :

      /**
        Constructor calls the reset member which reads the information about the boundary.
        */
      MapOperator (const boost::shared_ptr<OperatorParameters> & params) : Operator (params)
    { 
      boost::shared_ptr<MapOperatorParameters> myparams = 
        boost::dynamic_pointer_cast<MapOperatorParameters>(params);

      reset(myparams);   
    }

      /**
        Destructor
        */
      virtual ~MapOperator() { }

      virtual void reset(const boost::shared_ptr<OperatorParameters>& params);

      virtual AMP::LinearAlgebra::Variable::shared_ptr createInputVariable (const std::string & , int )
      {
        //Implemented in derived classes
        AMP::LinearAlgebra::Variable::shared_ptr emptyPointer;
        return emptyPointer;
      }

      virtual AMP::LinearAlgebra::Variable::shared_ptr createOutputVariable (const std::string & , int ) 
      {
        //Implemented in derived classes
        AMP::LinearAlgebra::Variable::shared_ptr emptyPointer;
        return emptyPointer;
      }

      virtual void setInputVariableName(const std::string & , int ) 
      {
        //Implemented in derived classes. 
      }

      virtual void setOutputVariableName(const std::string & , int )
      {
        //Implemented in derived classes. 
      }

      /*
         boost::shared_ptr<OperatorParameters> getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& )
         {
         boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
         tmp_db->putInteger("BoundaryId", d_boundaryId );
         boost::shared_ptr<MapOperatorParameters> outParams(new MapOperatorParameters(tmp_db));

         return outParams;
         }
         */

    protected :

      unsigned int d_boundaryId;

      AMP::Mesh::Mesh::shared_ptr d_MapMesh;

    private :

  };

}
}

#endif

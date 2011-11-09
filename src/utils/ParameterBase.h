
#ifndef included_AMP_ParameterBase
#define included_AMP_ParameterBase

#include <string>

#include "boost/shared_ptr.hpp"

namespace AMP{

  /**\class ParameterBase
   * ParameterBase is a base class for classes that pass parameters to other classes
   */
  class ParameterBase
  {
    public:

      typedef boost::shared_ptr<ParameterBase> shared_ptr;

      ParameterBase();
      virtual ~ParameterBase();

      std::string d_name;
    protected:

    private:
  };
}

#endif 

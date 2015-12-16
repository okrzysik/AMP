
#ifndef included_AMP_ParameterBase
#define included_AMP_ParameterBase

#include <string>

#include "utils/shared_ptr.h"

namespace AMP {

/**\class ParameterBase
 * ParameterBase is a base class for classes that pass parameters to other classes
 */
class ParameterBase {
public:
    typedef AMP::shared_ptr<ParameterBase> shared_ptr;

    ParameterBase();
    virtual ~ParameterBase();

    std::string d_name;

protected:
private:
};
}

#endif

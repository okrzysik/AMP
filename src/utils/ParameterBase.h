
#ifndef included_AMP_ParameterBase
#define included_AMP_ParameterBase

#include <memory>
#include <string>

#include "AMP/utils/Database.h"

namespace AMP {

/**\class ParameterBase
 * ParameterBase is a base class for classes that pass parameters to other classes
 */
class ParameterBase
{
public:
    typedef std::shared_ptr<ParameterBase> shared_ptr;

    ParameterBase() : d_name( "ParameterBase" ) {}
    virtual ~ParameterBase() = default;

    explicit ParameterBase( std::shared_ptr<AMP::Database> db ) : d_db( db ) {}

    std::string d_name;

    std::shared_ptr<AMP::Database> d_db;

protected:
private:
};
} // namespace AMP

#endif

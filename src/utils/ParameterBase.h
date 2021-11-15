#ifndef included_AMP_ParameterBase
#define included_AMP_ParameterBase

#include <string>


namespace AMP {

/**\class ParameterBase
 * ParameterBase is a base class for classes that pass parameters to other classes
 */
class ParameterBase
{
public:
    ParameterBase();
    virtual ~ParameterBase();

    std::string d_name;

protected:
private:
};

} // namespace AMP

#endif

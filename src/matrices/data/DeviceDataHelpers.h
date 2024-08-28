#ifndef included_AMP_DeviceDataHelpers
#define included_AMP_DeviceDataHelpers

namespace AMP {
namespace LinearAlgebra {

template<typename TYPE>
class DeviceDataHelpers
{
public:
    static void fill_n( TYPE *x, const size_t N, const TYPE alpha );
    static void copy_n( TYPE *x, const size_t N, TYPE *y );
};

} // namespace LinearAlgebra
} // namespace AMP

#endif

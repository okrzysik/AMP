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
    static void exclusive_scan( TYPE *x, const size_t N, TYPE *y, TYPE init );
    static TYPE max_element( TYPE *x, const size_t N );
    static TYPE accumulate( TYPE *x, const size_t N, TYPE init );
};

} // namespace LinearAlgebra
} // namespace AMP

#endif

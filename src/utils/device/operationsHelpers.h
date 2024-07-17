#ifndef included_AMP_OperationsHelpers_h
#define included_AMP_OperationsHelpers_h

namespace AMP {
  namespace LinearAlgebra {
    template <typename TYPE>
    class OperationsHelpers {
    public:
      static void copy_n(const TYPE *x, const size_t N, TYPE *y);
      static void fill_n(TYPE *x, const size_t N, const TYPE data);
    };
  }
}

#endif

#include "utils/Utilities.h"

namespace AMP {
namespace LinearAlgebra {

  inline
  void  Variable::setUnits ( const std::string &t )
  {
    d_Units = t;
  }

  inline
  const std::string &Variable::getUnits () const
  {
    return d_Units;
  }

  inline
  Variable::shared_ptr  Variable::getVariable ( size_t which )
  {
    AMP_ASSERT ( which == 0 );
    return shared_from_this ();
  }

}
}

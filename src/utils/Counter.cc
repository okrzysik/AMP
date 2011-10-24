
#include "Counter.h"

namespace AMP
{
  std::map<std::string , size_t> Counter::d_Data;

  void Counter::initCounter ()
  {
    d_Data["Virtual"] = 0;
    d_Data["FLOPS"] = 0;
  }

}

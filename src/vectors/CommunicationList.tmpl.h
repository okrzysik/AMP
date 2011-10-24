

namespace AMP {
namespace LinearAlgebra {

  template <typename T>
  T * CommunicationList::getBufferToAvoidDebugVectorCrashing ( std::vector<T> &in )
  {
    T *retVal = 0;
    if ( in.size() > 0 ) retVal = &(in[0]);
    return retVal;
  }

  template <typename T>
  const T * CommunicationList::getBufferToAvoidDebugVectorCrashing ( const std::vector<T> &in )
  {
    const T *retVal = 0;
    if ( in.size() > 0 ) retVal = &(in[0]);
    return retVal;
  }

}
}



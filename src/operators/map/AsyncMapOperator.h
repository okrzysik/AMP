#ifndef  included_AMP_AsyncMapOperator
#define  included_AMP_AsyncMapOperator

#include "operators/AsynchronousOperator.h"

namespace AMP {
namespace Operator {

  /** \brief  A base class for asynchronous map operations between meshes.  Maps
   * may impose a serial thread or even deadlock in parallel if implemented synchronously 
   * without great care. 
   */

  class AsyncMapOperator : public AsynchronousOperator
  {
    public:
      /** \brief  Constructor
      */
      AsyncMapOperator ( const boost::shared_ptr <OperatorParameters> & );

      virtual ~AsyncMapOperator ();

      /** \brief  Set a frozen vector for results of the apply operation.
       * \param[in]  p  The vector to set
       */
      virtual void setVector ( AMP::LinearAlgebra::Vector::shared_ptr &p ) = 0;

      bool isMaster() {
        return d_IsMaster;
      }

    protected:
      bool d_IsMaster;
  };

}
}


#endif

#ifndef included_AMP_AsynchronousOperator
#define included_AMP_AsynchronousOperator

#include "utils/AMP_MPI.h"
#include "operators/Operator.h"
#include "AsynchronousOperatorParameters.h"

#include <vector>

namespace AMP {
namespace Operator {


  /** \brief  A class that allows for overlapped communication
    */

  class AsynchronousOperator : public Operator
  {
    protected:

      /** \brief  A list of MPI_Requests for use in derived classes
        */
      std::vector<MPI_Request>         d_RequestList;

      /** \brief  Reserve a number of MPI_Requests for use
        * \param[in] i The number of MPI_Requests that will be used
        */
      void reserveRequests ( size_t i );

      /** \brief  Clear requests
        */
      void clearRequests ();

      /** \brief  Return an iterator to the first MPI_Request
        * \return The iterator
        */
      std::vector<MPI_Request>::iterator    beginRequests ();

      /** \brief  Return an iterator to one past the end of the list of requests
        * \return The iterator
        */
      std::vector<MPI_Request>::iterator    endRequests ();

      /** \brief  Return a specific MPI_Request
        * \return  The request
        */
      MPI_Request           &getRequest ( size_t i );

      /** \brief  Wait for all requests to complete
        */
      void                   waitForAllRequests ();


      /** \brief  Functions to allow std::vector to be used in MPI communications when empty
        * \tparam  T  The type of buffer to translate
        * \param[in] in Buffer in a std::vector
        * \return  Buffer of type T
        */
      template <typename T>
      T * getBufferToAvoidDebugVectorCrashing ( std::vector<T> &in )
      {
        T *retVal = 0;
        if ( in.size() > 0 ) retVal = &(in[0]);
        return retVal;
      }

      /** \brief  Functions to allow std::vector to be used in MPI communications when empty
        * \tparam  T  The type of buffer to translate
        * \param[in] in Buffer in a std::vector
        * \return  Buffer of type T
        */
      template <typename T>
      const T * getBufferToAvoidDebugVectorCrashing ( const std::vector<T> &in )
      {
        const T *retVal = 0;
        if ( in.size() > 0 ) retVal = &(in[0]);
        return retVal;
      }


    public:
      AsynchronousOperator ( const boost::shared_ptr < OperatorParameters > &params );
      virtual ~AsynchronousOperator ();

      /** \brief  Continue the construction of the object requiring asynchronous calls
        * \param[in] params  The parameters to pass to the continuation
        * \return True if construction is complete
        */
      virtual bool continueAsynchronousConstruction ( const boost::shared_ptr < OperatorParameters > &params );

      /** \brief  Start a communicative apply operation
        * \param[in]  f  An input vector
        * \param[in]  u  An input vector
        * \param[out] r  An output vector
        * \param[in]  a  A weight
        * \param[in]  b  A weight
        */
      virtual void applyStart(const AMP::LinearAlgebra::Vector::shared_ptr &f,
             const  AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr  &r,
             const double a = -1.0, const double b = 1.0) = 0;

      /** \brief  Finish a communicative apply operation
        * \param[in]  f  An input vector
        * \param[in]  u  An input vector
        * \param[out] r  An output vector
        * \param[in]  a  A weight
        * \param[in]  b  A weight
        */
      virtual void applyFinish(const AMP::LinearAlgebra::Vector::shared_ptr &f,
             const  AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr  &r,
             const double a = -1.0, const double b = 1.0) = 0;

      virtual void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f,
             const  AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr  &r,
             const double a = -1.0, const double b = 1.0);

  };

}
}


#endif

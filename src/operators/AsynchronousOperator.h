#ifndef included_AMP_AsynchronousOperator
#define included_AMP_AsynchronousOperator

#include "AMP/operators/Operator.h"
#include "AMP/utils/AMP_MPI.h"
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
    std::vector<MPI_Request> d_RequestList;

    /** \brief  Reserve a number of MPI_Requests for use
     * \param[in] i The number of MPI_Requests that will be used
     */
    void reserveRequests( size_t i );

    /** \brief  Clear requests
     */
    void clearRequests();

    /** \brief  Return an iterator to the first MPI_Request
     * \return The iterator
     */
    std::vector<MPI_Request>::iterator beginRequests();

    /** \brief  Return an iterator to one past the end of the list of requests
     * \return The iterator
     */
    std::vector<MPI_Request>::iterator endRequests();

    /** \brief  Return a specific MPI_Request
     * \return  The request
     */
    MPI_Request &getRequest( size_t i );

    /** \brief  Wait for all requests to complete
     */
    void waitForAllRequests();


    /** \brief  Functions to allow std::vector to be used in MPI communications when empty
     * \tparam  T  The type of buffer to translate
     * \param[in] in Buffer in a std::vector
     * \return  Buffer of type T
     */
    template<typename T>
    T *getBufferToAvoidDebugVectorCrashing( std::vector<T> &in )
    {
        T *retVal = 0;
        if ( in.size() > 0 )
            retVal = &( in[0] );
        return retVal;
    }

    /** \brief  Functions to allow std::vector to be used in MPI communications when empty
     * \tparam  T  The type of buffer to translate
     * \param[in] in Buffer in a std::vector
     * \return  Buffer of type T
     */
    template<typename T>
    const T *getBufferToAvoidDebugVectorCrashing( const std::vector<T> &in )
    {
        const T *retVal = 0;
        if ( in.size() > 0 )
            retVal = &( in[0] );
        return retVal;
    }


public:
    explicit AsynchronousOperator( std::shared_ptr<const OperatorParameters> params );

    virtual ~AsynchronousOperator();

    //! Return the name of the operator
    std::string type() const override { return "AsynchronousOperator"; }

    /** \brief  Start a communicative apply operation
     * \details  The specific meaning of applyStart will vary depending on the operator,
     *   but the intended purpose is to start non-blocking communication.
     * \param[in]  u  An input vector
     * \param[out] f  An output vector
     */
    virtual void applyStart( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                             AMP::LinearAlgebra::Vector::shared_ptr f ) = 0;

    /** \brief  Finish a communicative apply operation
     * \details  The specific meaning of applyFinish will vary depending on the operator,
     *   but the intended purpose is to finish non-blocking communication and any remaining
     *   operations for the vector
     * \param[in]  u  An input vector
     * \param[out] f  An output vector
     */
    virtual void applyFinish( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                              AMP::LinearAlgebra::Vector::shared_ptr f ) = 0;

    /** \brief  Apply operation
     * \details  The apply opertion for an asyncronous operator simply calls applyStart and
     * applyFinish.
     * \param[in]  u  An input vector
     * \param[out] f  An output vector
     */
    virtual void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                        AMP::LinearAlgebra::Vector::shared_ptr f ) override;
};
} // namespace Operator
} // namespace AMP


#endif

#ifndef  included_AMP_AsyncMapOperator
#define  included_AMP_AsyncMapOperator

#include "operators/AsynchronousOperator.h"
#include "ampmesh/Mesh.h"

namespace AMP {
namespace Operator {

/** \brief  A base class for asynchronous map operations between meshes.  
 * A map operation involves two meshes and a communicator spanning those meshes.
 * For some processors one of the meshes may be NULL.
 * The constructor may require syncronous communication, but the apply calls
 * should be implimented asynchronously.
 * Note: Maps may impose a serial thread or even deadlock in parallel if 
 * implemented synchronously without great care. 
 */
class AsyncMapOperator : public AsynchronousOperator
{
public:
    //! Constructor
    AsyncMapOperator ( const boost::shared_ptr <OperatorParameters> & );

    virtual ~AsyncMapOperator ();

    /** \brief  Set a frozen vector for results of the apply operation.
     * \param[in]  p  The vector to set
     */
    virtual void setVector ( AMP::LinearAlgebra::Vector::shared_ptr &p ) = 0;

protected:
    AMP_MPI d_comm;
    AMP::Mesh::Mesh::shared_ptr d_mesh1;
    AMP::Mesh::Mesh::shared_ptr d_mesh2;
};

}
}


#endif

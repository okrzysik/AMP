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
    AsyncMapOperator ( const AMP::shared_ptr <OperatorParameters> & );

    virtual ~AsyncMapOperator ();

    /** \brief  Set a frozen vector for results of the apply operation.
     * \param[in]  p  The vector to set
     */
    virtual void setVector ( AMP::LinearAlgebra::Vector::shared_ptr p ) = 0;

    // Overload the apply operator to include makeConsistent
    virtual void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u, 
			AMP::LinearAlgebra::Vector::shared_ptr f );    
    
    // Function to determine if a makeConsistentSet is required
    virtual bool requiresMakeConsistentSet();

    /**
     * Get the meshes this map uses.
     * \param[in] which 1 for d_mesh1 and 2 for d_mesh2
     */
    AMP::Mesh::Mesh::shared_ptr getMesh(int which); 

    virtual AMP::LinearAlgebra::Variable::shared_ptr  getInputVariable () { return d_inpVariable; }
    virtual AMP::LinearAlgebra::Variable::shared_ptr  getOutputVariable () { return d_outVariable; }

protected:

    // Communicator for the Map
    AMP_MPI d_MapComm;

    // Variables to store the individual meshes and the DOFManager
    AMP::Mesh::Mesh::shared_ptr  d_mesh1;
    AMP::Mesh::Mesh::shared_ptr  d_mesh2;
    AMP::Discretization::DOFManager::shared_ptr  d_DOFManager;
    AMP::LinearAlgebra::Variable::shared_ptr      d_inpVariable;
    AMP::LinearAlgebra::Variable::shared_ptr      d_outVariable;

    // Frozen vector for the output results
    AMP::LinearAlgebra::Vector::shared_ptr  d_OutputVector;
};


}
}


#endif

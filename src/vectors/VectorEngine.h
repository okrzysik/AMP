#ifndef included_AMP_VectorEngine
#define included_AMP_VectorEngine

#include "AMP/utils/AMP_MPI.h"
#include "AMP/vectors/data/VectorData.h"
#include "AMP/vectors/operations/VectorOperationsDefault.h"
#include <memory>
#include <vector>


namespace AMP {
namespace LinearAlgebra {

/**
 * \brief The parameters necessary to build a VectorEngine.
 * \details  This is currently empty since there is no default
 * engine...yet
 */

class VectorEngineParameters
{
public:
    typedef std::shared_ptr<VectorEngineParameters> shared_ptr;

    /** \brief Constructor
        \param[in] local_size     The number of elements on this core
        \param[in] global_size    The number of elements in total
        \param[in] comm           Communicator to create the vector on
        \details  This assumes a contiguous allocation of data.
                  Core 0 has global ids \f$(0,1,\ldots,n-1)\f$,
                  core 1 has global ids \f$(n,n+1,n+2,\ldots,m)\f$, etc.
      */
    VectorEngineParameters( size_t local_size, size_t global_size, const AMP_MPI &comm );

    //! destructor
    virtual ~VectorEngineParameters();

    //! Return the local size
    inline size_t getLocalSize() const { return d_end - d_begin; }

    //! Return the local size
    inline size_t getGlobalSize() const { return d_global; }

    //! Return the first DOF on this core
    inline size_t beginDOF() const { return d_begin; }

    /** \brief  Return the Epetra_MpiComm for this engine
     * \return The Epetra_MpiComm
     */
    inline AMP_MPI getComm() const { return d_comm; }

protected:
    size_t d_begin;  // Starting DOF
    size_t d_end;    // Ending DOF
    size_t d_global; // Number of global DOFs
    AMP_MPI d_comm;  // Comm
};


/** \class VectorEngine
 * \brief A class that can perform mathematics on vectors.
 * \see Vector
 * \details  This class will eventually house the mechanics of performing
 * math on vectors.  There is currently a large overlap between this class
 * and Vector.  Eventually, Vector will completely encapsulate data storage
 * and access while this class will completely encapsulate dense kernels on
 * the data.
 */
class VectorEngine : virtual public VectorData, virtual public VectorOperations
{
protected:
    VectorEngineParameters::shared_ptr d_Params;

public:
    /** \brief  Destructor
     */
    virtual ~VectorEngine();

    /** \brief Allocate a new buffer
     * \return A shared pointer to a new buffer
     */
    virtual std::shared_ptr<VectorData> getNewBuffer() = 0;

    /** \brief  Return a copy of this engine
     * \param[in]  p  The buffer to use for the copy.
     * \return  The new engine
     */
    virtual std::shared_ptr<VectorEngine> cloneEngine( std::shared_ptr<VectorData> p ) const = 0;

    /** \brief Swap engines
     * \param[in,out] p  The engine to exchange with
     */
    virtual void swapEngines( std::shared_ptr<VectorEngine> p ) = 0;

    /** \brief  Return the communicator associated with this engine
     * \return  The communicator associated with this engine
     */
    virtual AMP_MPI getComm() const = 0;
};


} // namespace LinearAlgebra
} // namespace AMP

#include "VectorEngine.inline.h"
#endif

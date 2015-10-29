#ifndef included_AMP_VectorEngine
#define included_AMP_VectorEngine


#include "utils/shared_ptr.h"
#include <vector>
#include "VectorOperations.h"
#include "utils/AMP_MPI.h"


namespace AMP {
namespace LinearAlgebra {

/**
 * \brief The parameters necessary to build a VectorEngine.
 * \details  This is currently empty since there is no default
 * engine...yet
*/

class VectorEngineParameters : public Castable
{
public:
    typedef  AMP::shared_ptr<VectorEngineParameters>    shared_ptr;

    /** \brief Constructor
        \param[in] local_size     The number of elements on this core
        \param[in] global_size    The number of elements in total
        \param[in] comm           Communicator to create the vector on
        \details  This assumes a contiguous allocation of data.
                  Core 0 has global ids \f$(0,1,\ldots,n-1)\f$,
                  core 1 has global ids \f$(n,n+1,n+2,\ldots,m)\f$, etc.
      */
     VectorEngineParameters( size_t local_size, size_t global_size, AMP_MPI comm );

    //! destructor
    virtual ~VectorEngineParameters();

    //! Return the local size
    inline size_t getLocalSize() const { return d_end-d_begin; }

    //! Return the local size
    inline size_t getGlobalSize() const { return d_global; }

    //! Return the first DOF on this core
    inline size_t beginDOF() const { return d_begin; }

    //! Return 1 past the last DOF on this core
    inline size_t endDOF() const { return d_end; }

    /** \brief  Return the Epetra_MpiComm for this engine
      * \return The Epetra_MpiComm
      */
    inline AMP_MPI getComm() const { return d_comm; }

protected:
    size_t                     d_begin;     // Starting DOF
    size_t                     d_end;       // Ending DOF
    size_t                     d_global;    // Number of global DOFs
    AMP_MPI                    d_comm;      // Comm
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
class VectorEngine : virtual public VectorOperations
{
protected:
    VectorEngineParameters::shared_ptr        d_Params;

public:
    //! The basic buffer type of vectors
    typedef std::vector<double>             Buffer;

    //! A shared pointer to the buffer
    typedef AMP::shared_ptr<Buffer>         BufferPtr;

    //! A shared pointer to an engine
    typedef AMP::shared_ptr<VectorEngine>     shared_ptr;

    //! A const shared pointer to an engine
    typedef AMP::shared_ptr<const VectorEngine> const_shared_ptr;

    /** \brief  Destructor
      */
    virtual  ~VectorEngine ();

    /** \brief Allocate a new buffer
      * \return A shared pointer to a new buffer
      */
    virtual   BufferPtr   getNewBuffer () = 0;

    /** \brief  True if engines are the same
      * \param[in]  e  Engine to compare against
      * \return True if the engine is the same type as this
      */
    virtual   bool      sameEngine ( VectorEngine &e ) const = 0;

    /** \brief  Return a copy of this engine
      * \param[in]  p  The buffer to use for the copy.
      * \return  The new engine
      */
    virtual   shared_ptr  cloneEngine ( BufferPtr p ) const = 0;

    /** \brief Swap engines
      * \param[in,out] p  The engine to exchange with
      */
    virtual   void      swapEngines ( shared_ptr p ) = 0;

    /** \brief Return the number of contiguous blocks associated with this engine
      * \return  The number of contiguous blocks associated with this engine
      */
    virtual   size_t    numberOfDataBlocks () const = 0;

    /** \brief Get the size of a data block
      * \param[in]  i  Which block to measure size of
      * \return  The size of the ith block
      */
    virtual   size_t    sizeOfDataBlock ( size_t i ) const = 0;

    /** \brief  The number of locally owned elements
      * \return  The number of elements owned by this core
      */
    virtual   size_t    getLocalSize() const = 0;

    /** \brief  The number of elements this engine will use
      * \return The number of elements over all cores
      */
    virtual   size_t    getGlobalSize() const = 0;

    /** \brief  Copy data into the engine's buffer
      * \param[in] in  The data to copy in
      */
    virtual   void      putRawData ( const double *in ) = 0;

    /** \brief Return a contiguous block of data
      */
    virtual   void     *getDataBlock ( size_t i ) = 0;

    /** \brief Return a contiguous block of data
      */
    virtual   const void *getDataBlock ( size_t i ) const = 0;

    /** \brief  Get the parameters used to create this engine
      * \return The parameters
      */
    virtual   VectorEngineParameters::shared_ptr    getEngineParameters() const;

    /** \brief Set values in the engine's buffer
      * \param[in]  i  The number of values to set
      * \param[in]  ndx  The entries in the buffer to set given by local id
      * \param[in]  val  The values to copy in
      * \details  This will not set ghost values
      */
    virtual void setValuesByLocalID(int i, size_t *ndx , const double *val) = 0;

    /** \brief Set values in the engine's buffer
      * \param[in]  i  The number of values to set
      * \param[in]  ndx  The entries in the buffer to set given by global id
      * \param[in]  val  The values to copy in
      * \details  This will not set ghost values
      */
    virtual void setLocalValuesByGlobalID(int i, size_t *ndx , const double *val) = 0;

    /** \brief Add values in the engine's buffer
      * \param[in]  i  The number of values to set
      * \param[in]  ndx  The entries in the buffer to set given by local id
      * \param[in]  val  The values to copy in
      * \details  This will not set ghost values
      */
    virtual void addValuesByLocalID(int i, size_t *ndx , const double *val) = 0;

    /** \brief Add values in the engine's buffer
      * \param[in]  i  The number of values to set
      * \param[in]  ndx  The entries in the buffer to set given by global id
      * \param[in]  val  The values to copy in
      * \details  This will not set ghost values
      */
    virtual void addLocalValuesByGlobalID(int i, size_t *ndx , const double *val) = 0;

    /** \brief Get values in the engine's buffer
      * \param[in]  i  The number of values to set
      * \param[in]  ndx  The entries in the buffer to set given by local id
      * \param[out]  val  The values requested
      * \details  This will not set ghost values
      */
    virtual void getValuesByLocalID(int i, size_t *ndx ,double *val) const = 0;

    /** \brief Get values in the engine's buffer
      * \param[in]  i  The number of values to set
      * \param[in]  ndx  The entries in the buffer to set given by global id
      * \param[out]  val  The values requested
      * \details  This will not set ghost values
      */
    virtual void getLocalValuesByGlobalID(int i, size_t *ndx ,double *val) const = 0;

    /** \brief  Return the communicator associated with this engine
      * \return  The communicator associated with this engine
      */
    virtual AMP_MPI   getComm () const = 0;

    /** \brief  Allocate a buffer of the right size and copy out the data in the engine's buffer
      * \param[out] out The buffer to create and copy to
      * \details  The caller of this method is responsible for deleting the array *out
      */
    virtual void copyOutRawData ( double *out ) const = 0;

};


}
}

#include "VectorEngine.inline.h"
#endif


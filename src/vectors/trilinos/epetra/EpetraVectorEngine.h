#ifndef included_AMP_EpetraVectorEngine
#define included_AMP_EpetraVectorEngine

#include <Epetra_Map.h>
#include <Epetra_Vector.h>

#include "EpetraVector.h"

#include "vectors/VectorEngine.h"
#include "vectors/trilinos/epetra/EpetraVectorData.h"
#include "vectors/trilinos/epetra/EpetraVectorOperations.h"


namespace AMP {
namespace LinearAlgebra {

/** \class EpetraVectorEngineParameters
 * \brief Class that details how to construct an EpetraVectorEngine
 */
class EpetraVectorEngineParameters : public VectorEngineParameters
{
public:
    /** \brief Constructor
        \param[in] local_size     The number of elements on this core
        \param[in] global_size    The number of elements in total
        \param[in] comm         Communicator to create the vector on
        \details  This assumes a contiguous allocation of data.  Core 0 has global ids
       \f$(0,1,\ldots,n-1)\f$, core 1
       has global ids \f$(n,n+1,n+2,\ldots,m)\f$, etc.
        */
    EpetraVectorEngineParameters( size_t local_size, size_t global_size, AMP_MPI comm );

    /** \brief Constructor
     * \param[in]  local_size    The number of elements on this core
     * \param[in]  global_size   The number of elements in total
     * \param[in]  emap        An Epetra_Map for the data
     * \param[in]  ecomm       An Epetra_MpiComm for constructing the vector on
     * \details  This allows construction of an EpetraVectorEngine from handy Epetra objects
     */
    EpetraVectorEngineParameters( size_t local_size,
                                  size_t global_size,
                                  AMP::shared_ptr<Epetra_Map> emap,
                                  AMP_MPI ecomm );

    //! Destructor
    virtual ~EpetraVectorEngineParameters();

    /** \brief  Return the Epetra_Map for this engine
     * \return  The Epetra_Map
     */
    Epetra_Map &getEpetraMap();

private:
    AMP::shared_ptr<Epetra_Map> d_emap; // Epetra map
};


/** \class EpetraVectorEngine
 * \brief A linear algebra engine that uses Epetra
 * \details  Use the Epetra implementation of the L1 BLAS routines.  Unlike other
 * libraries, it is very difficult to separate the data from the engine.  For this
 * reason, the EpetraVectorEngine contains the Epetra_Vector to operate on.
 */
class EpetraVectorEngine : public VectorEngine,
                           public EpetraVectorData,
                           public EpetraVectorOperations
{


public:
    /** \brief Constructor
     * \param[in]  alias  The parameters to construct this engine
     * \param[in]  p  The buffer to use to construct the engine
     */
    EpetraVectorEngine( AMP::shared_ptr<VectorEngineParameters> alias,
        AMP::shared_ptr<VectorData> p = nullptr );

    /** \brief Destructor
     */
    virtual ~EpetraVectorEngine();

    /** \brief  Get the raw Epetra_Vector
     * \return  The Epetra_Vector currently used by this engine
     */
    Epetra_Vector &getEpetra_Vector();

    /** \brief  Get the raw Epetra_Vector
     * \return  The Epetra_Vector currently used by this engine
     */
    const Epetra_Vector &getEpetra_Vector() const;


public: // Functions derived from VectorEngine
    AMP_MPI getComm() const override;
    virtual bool sameEngine( VectorEngine &e ) const override;
    virtual AMP::shared_ptr<VectorEngine> cloneEngine( AMP::shared_ptr<VectorData> p ) const override;
    virtual void swapEngines( AMP::shared_ptr<VectorEngine> ) override;
    virtual AMP::shared_ptr<VectorData> getNewBuffer() override;
};


} // namespace LinearAlgebra
} // namespace AMP


#include "EpetraVectorEngine.inline.h"

#endif

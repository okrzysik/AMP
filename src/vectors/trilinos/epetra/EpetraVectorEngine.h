#ifndef included_AMP_EpetraVectorEngine
#define included_AMP_EpetraVectorEngine

#include <Epetra_Map.h>
#include <Epetra_Vector.h>

#include "AMP/vectors/trilinos/epetra/EpetraVectorData.h"
#include "AMP/vectors/trilinos/epetra/EpetraVectorOperations.h"
#include "EpetraVector.h"


namespace AMP {
namespace LinearAlgebra {


/** \class EpetraVectorEngine
 * \brief A linear algebra engine that uses Epetra
 * \details  Use the Epetra implementation of the L1 BLAS routines.  Unlike other
 * libraries, it is very difficult to separate the data from the engine.  For this
 * reason, the EpetraVectorEngine contains the Epetra_Vector to operate on.
 */
class EpetraVectorEngine : public Vector
{


public:
    /** \brief Constructor
     * \param[in]  alias  The parameters to construct this engine
     * \param[in]  p  The buffer to use to construct the engine
     */
    explicit EpetraVectorEngine( std::shared_ptr<EpetraVectorEngineParameters> alias,
                                 std::shared_ptr<VectorData> p = nullptr );

    /** \brief Destructor
     */
    virtual ~EpetraVectorEngine() {}

    /** \brief  Get the raw Epetra_Vector
     * \return  The Epetra_Vector currently used by this engine
     */
    inline Epetra_Vector &getEpetra_Vector() { return std::dynamic_pointer_cast<EpetraVectorData >(d_VectorData)->getEpetra_Vector(); }

    /** \brief  Get the raw Epetra_Vector
     * \return  The Epetra_Vector currently used by this engine
     */
    inline const Epetra_Vector &getEpetra_Vector() const { return std::dynamic_pointer_cast<EpetraVectorData>(d_VectorData)->getEpetra_Vector(); }

public: // Functions derived from VectorData
    AMP_MPI getComm() const override { return d_Params->getComm(); }

public: // Functions derived from Vector
    using Vector::cloneVector;

    std::string type() const override { return "EpetraVectorEngine"; }
    Vector::shared_ptr cloneVector( const Variable::shared_ptr name ) const override;
    void swapVectors( Vector &other ) override;
    void aliasVector( Vector &other ) override;
    void assemble() override;

protected:
    std::shared_ptr<EpetraVectorEngineParameters> d_Params;
};


} // namespace LinearAlgebra
} // namespace AMP


#endif

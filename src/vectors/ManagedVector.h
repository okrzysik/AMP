#ifndef included_AMP_ManagedVector
#define included_AMP_ManagedVector

#include "AMP/vectors/ManagedVectorData.h"
#include "AMP/vectors/Vector.h"

#include <stdexcept>
#include <vector>


namespace AMP {
namespace LinearAlgebra {

/**
   \brief Class used to control data and kernels of various vector libraries
   \details  A ManagedVector will take an engine and create a buffer, if
   necessary.

   A ManagedVector has two pointers: data and engine.  If the data pointer
   is null, then the engine is assumed to have the data.
*/
class ManagedVector : public Vector
{

public:
    /** \brief Construct a ManagedVector from a set of parameters
     * \param[in] params  The description of the ManagedVector
     */
    explicit ManagedVector( std::shared_ptr<ManagedVectorParameters> params );

    /** \brief Construct a view of an AMP vector
     * \param[in] alias  Vector to view
     */
    explicit ManagedVector( const Vector::shared_ptr alias );

    //! Destructor
    virtual ~ManagedVector();

    /** \brief  If a vector has multiple views to multiple external packages
     * associated with it, this will return the barest version of the vector
     * \return A vector with the fewest views associated with it.
     * \details  A ManagedVector must have an engine and it may have data.
     * If it has an engine with no data, then the engine has must have data.
     * If the engine can be cast to a ManagedVector, it is and getRootVector
     * is called recursively.
     */
    Vector::shared_ptr getRootVector();

    /** \brief  Return the engine associated with this ManagedVector
     * \return The engine
     */
    std::shared_ptr<Vector> getVectorEngine();
    std::shared_ptr<const Vector> getVectorEngine() const;

    virtual bool isAnAliasOf( Vector &rhs );
    virtual bool isAnAliasOf( Vector::shared_ptr rhs );

protected:
    //! The parameters used to create this vector
    std::shared_ptr<ManagedVectorParameters> d_pParameters;

    //! Function that returns a pointer to a managed vector
    virtual ManagedVector *getNewRawPtr() const = 0;

public: // Derived from Vector
    std::string type() const override;
    std::shared_ptr<Vector> cloneVector( const Variable::shared_ptr name ) const override;
    virtual std::shared_ptr<ManagedVectorParameters> getParameters() { return d_pParameters; }
    Vector::shared_ptr subsetVectorForVariable( Variable::const_shared_ptr name ) override;
    Vector::const_shared_ptr
    constSubsetVectorForVariable( Variable::const_shared_ptr name ) const override;
    void swapVectors( Vector &other ) override;

protected: // Derived from Vector
    Vector::shared_ptr selectInto( const VectorSelector & ) override;
    Vector::const_shared_ptr selectInto( const VectorSelector & ) const override;

private:
    ManagedVector();
};


} // namespace LinearAlgebra
} // namespace AMP


#endif

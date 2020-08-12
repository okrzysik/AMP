#ifndef included_AMP_NullVector
#define included_AMP_NullVector

#include "AMP/vectors/Vector.h"
#include "AMP/vectors/data/VectorDataNull.h"
#include <string>


namespace AMP {
namespace LinearAlgebra {


/** \brief An empty vector
 * \details Some operators do not require vectors for application.  In these
 * circumstances, a NullVector is used.  This stores no data and performs no
 * work.
 */
template<class TYPE = double>
class NullVector : public Vector, public VectorDataNull<TYPE>
{
public: // Public constructors
    /**
     *  \brief Create a NullVector
     *  \param[in]  name  Name of variable to associate with this NullVector
     *  \return Vector shared pointer to a NullVector
     */
    static inline Vector::shared_ptr create( const std::string &name )
    {
        return create( std::make_shared<Variable>( name ) );
    }

    /**
     *  \brief Create a NullVector
     *  \param[in]  name  Variable to associate with this NullVector
     *  \return Vector shared pointer to a NullVector
     */
    static inline Vector::shared_ptr create( const Variable::shared_ptr name )
    {
        return std::shared_ptr<NullVector>( new NullVector<TYPE>( name ) );
    }

    virtual ~NullVector() = default;


public: // Functions inherited from Vector
    inline std::string type() const override { return "Null Vector"; }
    inline std::shared_ptr<ParameterBase> getParameters() override
    {
        return std::shared_ptr<ParameterBase>();
    }
    inline shared_ptr cloneVector( const Variable::shared_ptr name ) const override
    {
        return create( name );
    }
    inline void swapVectors( Vector & ) override {}
    inline void aliasVector( Vector & ) override {}
    inline void makeConsistent( ScatterType ) override {}
    inline void assemble() override {}
    using Vector::cloneVector;
    using Vector::dot;


protected:
    virtual Vector::shared_ptr selectInto( const VectorSelector & ) override
    {
        return Vector::shared_ptr();
    }
    virtual Vector::const_shared_ptr selectInto( const VectorSelector & ) const override
    {
        return Vector::const_shared_ptr();
    }


private:
    explicit inline NullVector( Variable::shared_ptr var )
    {
        setVariable( var );
        d_CommList = CommunicationList::createEmpty( 0, AMP_MPI( AMP_COMM_SELF ) );
        d_DOFManager.reset( new AMP::Discretization::DOFManager( 0, AMP_MPI( AMP_COMM_SELF ) ) );
    }
};


} // namespace LinearAlgebra
} // namespace AMP

#endif

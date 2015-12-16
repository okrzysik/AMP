#ifndef included_AMP_NullVector
#define included_AMP_NullVector

#include "Vector.h"
#include <string>

namespace AMP {
namespace LinearAlgebra {


/** \brief An empty vector
  * \details Some operators do not require vectors for application.  In these
  * circumstances, a NullVector is used.  This stores no data and performs no
  * work.
  */
class NullVector : public Vector
{
private:
    explicit NullVector( Variable::shared_ptr );

public:
    /**
      *  \brief Create a NullVector
      *  \param[in]  name  Name of variable to associate with this NullVector
      *  \return Vector shared pointer to a NullVector
      */
    static Vector::shared_ptr create( const std::string &name );

    /**
      *  \brief Create a NullVector
      *  \param[in]  name  Variable to associate with this NullVector
      *  \return Vector shared pointer to a NullVector
      */
    static Vector::shared_ptr create( const Variable::shared_ptr name );

    virtual ~NullVector();

    virtual std::string type() const override { return "Null Vector"; }
    virtual AMP::shared_ptr<ParameterBase> getParameters() override;

    virtual shared_ptr cloneVector( const Variable::shared_ptr name ) const override;

    virtual void copyVector( Vector::const_shared_ptr rhs ) override;
    virtual void swapVectors( Vector & ) override;
    virtual void aliasVector( Vector & ) override;

    virtual void setToScalar( double ) override;
    virtual void scale( double, const VectorOperations & ) override;
    virtual void scale( double ) override;
    virtual void add( const VectorOperations &, const VectorOperations & ) override;
    virtual void subtract( const VectorOperations &, const VectorOperations & ) override;
    virtual void multiply( const VectorOperations &, const VectorOperations & ) override;
    virtual void divide( const VectorOperations &, const VectorOperations & ) override;
    virtual void reciprocal( const VectorOperations & ) override;
    virtual void
    linearSum( double, const VectorOperations &, double, const VectorOperations & ) override;
    virtual void axpy( double, const VectorOperations &, const VectorOperations & ) override;
    virtual void axpby( double, double, const VectorOperations & ) override;
    virtual void abs( const VectorOperations & ) override;
    virtual double min( void ) const override;
    virtual double max( void ) const override;
    virtual void setRandomValues( void ) override;
    virtual void setValuesByLocalID( int, size_t *, const double * ) override;
    virtual void setLocalValuesByGlobalID( int, size_t *, const double * ) override;
    virtual void addValuesByLocalID( int, size_t *, const double * ) override;
    virtual void addLocalValuesByGlobalID( int, size_t *, const double * ) override;
    virtual void getLocalValuesByGlobalID( int, size_t *, double * ) const override;

    virtual void makeConsistent( ScatterType ) override;

    virtual void assemble() override;
    virtual double L1Norm( void ) const override;
    virtual double L2Norm( void ) const override;
    virtual double maxNorm( void ) const override;
    virtual double dot( const VectorOperations & ) const override;

    virtual void putRawData( const double *in ) override;
    virtual void copyOutRawData( double *out ) const override;

    virtual size_t getLocalSize() const override;
    virtual size_t getGlobalSize() const override;
    virtual size_t getGhostSize() const override;

    virtual size_t numberOfDataBlocks() const;
    virtual size_t sizeOfDataBlock( size_t ) const;

    virtual uint64_t getDataID() const override { return 0; }

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

    virtual void *getRawDataBlockAsVoid( size_t ) override;
    virtual const void *getRawDataBlockAsVoid( size_t ) const override;
};
}
}

#endif

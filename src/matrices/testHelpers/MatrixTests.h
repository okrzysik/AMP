#ifndef included_AMP_test_MatrixTests
#define included_AMP_test_MatrixTests

#include "AMP/matrices/Matrix.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"

#include <memory>
#include <string>


namespace AMP::LinearAlgebra {


/**
 * \class MatrixFactory
 * \brief A helper class to generate vectors
 */
class MatrixFactory
{
public:
    virtual ~MatrixFactory() {}
    virtual std::shared_ptr<AMP::Mesh::Mesh> getMesh() const                    = 0;
    virtual std::shared_ptr<AMP::LinearAlgebra::Vector> getVector() const       = 0;
    virtual std::shared_ptr<AMP::LinearAlgebra::Matrix> getMatrix() const       = 0;
    virtual std::shared_ptr<AMP::Discretization::DOFManager> getDOFMap() const  = 0;
    virtual std::shared_ptr<AMP::Discretization::DOFManager> getDOFMapL() const = 0;
    virtual std::string name() const                                            = 0;
    virtual std::string type() const                                            = 0;

protected:
    MatrixFactory() {}
    MatrixFactory( const MatrixFactory & );
};


/**
 * \class vectorTests
 * \brief A helper class to store/run tests for a vector
 */
class MatrixTests
{
public:
    explicit MatrixTests( std::shared_ptr<const MatrixFactory> factory ) : d_factory( factory ) {}

public:
    void InstantiateMatrix( AMP::UnitTest *ut );
    void VerifyGetLeftRightVector( AMP::UnitTest *ut );
    void VerifyGetSetValuesMatrix( AMP::UnitTest *ut );
    void VerifyAXPYMatrix( AMP::UnitTest *ut );
    void VerifyScaleMatrix( AMP::UnitTest *ut );
    void VerifyExtractDiagonal( AMP::UnitTest *ut );
    void VerifyMultMatrix( AMP::UnitTest *ut );
    void VerifyMatMultMatrix( AMP::UnitTest *ut );
    void VerifyAddElementNode( AMP::UnitTest *ut );

private:
    std::shared_ptr<const MatrixFactory> d_factory;
};


} // namespace AMP::LinearAlgebra


// Extra includes
#include "AMP/vectors/testHelpers/VectorTests.inline.h"


#endif

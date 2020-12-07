#ifndef included_AMP_test_MatrixTests
#define included_AMP_test_MatrixTests

#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/matrices/Matrix.h"

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
    virtual AMP::Mesh::Mesh::shared_ptr getMesh() const = 0;
    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const = 0;
    virtual AMP::LinearAlgebra::Matrix::shared_ptr getMatrix() const = 0;
    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const = 0;
    virtual AMP::Discretization::DOFManager::shared_ptr getDOFMapL() const = 0;
    virtual std::string name() const                                 = 0;
    virtual std::string type() const { return "auto"; }
    virtual void initMesh()  = 0;
    virtual void endMesh()  = 0;

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


}


// Extra includes
#include "AMP/vectors/testHelpers/VectorTests.inline.h"


#endif

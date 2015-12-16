
#ifndef included_AMP_NonlinearFEOperator
#define included_AMP_NonlinearFEOperator

/* AMP files */
#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"
#include "libmesh/elem.h"
#include "operators/Operator.h"
#include "operators/libmesh/FEOperatorParameters.h"

#include <vector>

namespace AMP {
namespace Operator {

/**
  An abstract base class for representing a nonlinear finite element (FE) operator.
  This class implements a "shell" for computing the FE residual vector.
  This class only deals with the volume integration, the boundary conditions are
  handled separately by the boundary operators.
  */
class NonlinearFEOperator : public Operator
{
public:
    //! Constructor. This copies the share pointer to the element operation from the input parameter
    //! object.
    explicit NonlinearFEOperator( const AMP::shared_ptr<FEOperatorParameters> &params );

    //! Destructor
    virtual ~NonlinearFEOperator();

    /**
      The apply function for this operator, A, performs the following operation:
      f = A(u)
      @param [in] u input vector.
      @param [out] f residual/output vector.
      */
    virtual void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                        AMP::LinearAlgebra::Vector::shared_ptr f ) override;

protected:
    /**
      This function will be called just before looping over the elements to form the residual
      vector, so if the
      derived classes need to perform some initialization operations just before looping over the
      elements they can
      do so by implementing these operations in this function.
      Also, the derived classes can access the input (u) and output (r) vectors
      passed to the apply function by implementing this function.
      @param [in] u Input vector
      @param [out] r Output vector
      */
    virtual void preAssembly( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                              AMP::LinearAlgebra::Vector::shared_ptr r ) = 0;

    /**
      This function will be called just after looping over all the elements to form the residual
      vector, so if the
      derived classes need to perform any operations such as freeing any temporary memory
      allocations they
      can do so by implementing these operations in this function.
      */
    virtual void postAssembly() = 0;

    /**
      This function will be called once for each element, just before performing the element
      operation.
      Ideally, the element operation should not deal with global mesh related information such as
      DOFMap and global vectors and matrices. This function typically extracts the local information
      from
      these global objects and passes them to the element operation.
      */
    virtual void preElementOperation( const AMP::Mesh::MeshElement & ) = 0;

    /**
      This function will be called once for each element, just after performing the element
      operation.
      Typically, the result of the element operation is added to the global output vector in this
      function.
      */
    virtual void postElementOperation() = 0;

    void createLibMeshElementList();

    void destroyLibMeshElementList();

    std::vector<::Elem *> d_currElemPtrs;

    size_t d_currElemIdx;

    AMP::shared_ptr<ElementOperation> d_elemOp; /**< Shared pointer to the element operation */
};
}
}

#endif

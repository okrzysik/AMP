
#ifndef included_AMP_NonlinearFEOperator
#define included_AMP_NonlinearFEOperator

/* AMP files */
#include "operators/Operator.h"
#include "FEOperatorParameters.h"
#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"
#include "elem.h"

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
      public :

        /**
          Constructor. This copies the share pointer to the element operation from the input parameter object.
          */
        NonlinearFEOperator(const boost::shared_ptr<FEOperatorParameters>& params)
          : Operator(params)
        {
          d_elemOp = (params->d_elemOp);
          createLibMeshElementList();
          d_currElemIdx = static_cast<unsigned int>(-1);
        }

        /**
          Destructor
          */
        virtual ~NonlinearFEOperator() { 
          destroyLibMeshElementList();
        }

        /**
          The apply function for this operator, A, performs the following operation:
          r = b*f+a*A(u), if f is not NULL and r = a*A(u), if f is NULL.
          @param [in] f auxillary/rhs vector. 
          @param [in] u input vector. 
          @param [out] r residual/output vector. 
          @param [in] a first constant used in the expression: r = a*A(u) + b*f. The default value is -1.
          @param [in] b second constant used in the expression: r = a*A(u) + b*f. The default value is 1.
          */
        virtual void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
            AMP::LinearAlgebra::Vector::shared_ptr &r, const double a = -1.0, const double b = 1.0);

      protected :

        /**
          This function will be called just before looping over the elements to form the residual vector, so if the
          derived classes need to perform some initialization operations just before looping over the elements they can 
          do so by implementing these operations in this function.
          Also, the derived classes can access the input (u) and output (r) vectors 
          passed to the apply function by implementing this function.
          @param [in] u Input vector
          @param [out] r Output vector
          */
        virtual void preAssembly(const boost::shared_ptr<AMP::LinearAlgebra::Vector> &u, boost::shared_ptr<AMP::LinearAlgebra::Vector>  &r)=0;

        /**
          This function will be called just after looping over all the elements to form the residual vector, so if the
          derived classes need to perform any operations such as freeing any temporary memory allocations they
          can do so by implementing these operations in this function.
          */
        virtual void postAssembly()=0;

        /**
          This function will be called once for each element, just before performing the element operation.
          Ideally, the element operation should not deal with global mesh related information such as 
          DOFMap and global vectors and matrices. This function typically extracts the local information from
          these global objects and passes them to the element operation.
          */
        virtual void preElementOperation(const AMP::Mesh::MeshElement &)=0;

        /**
          This function will be called once for each element, just after performing the element operation.
          Typically, the result of the element operation is added to the global output vector in this function.
          */
        virtual void postElementOperation()=0;

        void createLibMeshElementList();

        void destroyLibMeshElementList();

        std::vector< ::Elem* > d_currElemPtrs;

        size_t d_currElemIdx;

        boost::shared_ptr<ElementOperation> d_elemOp; /**< Shared pointer to the element operation */

    };

  } 
}

#endif


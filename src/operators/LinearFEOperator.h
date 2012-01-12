
#ifndef included_AMP_LinearFEOperator
#define included_AMP_LinearFEOperator

/* AMP files */
#include "LinearOperator.h"
#include "FEOperatorParameters.h"
#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"

namespace AMP {
  namespace Operator {

    /**
      An abstract base class for representing a linear finite element (FE) operator.
      This class implements a "shell" for building the FE stiffness matrix. 
      This class only deals with the volume integration, the boundary conditions are 
      handled separately by the boundary operators. 
      */
    class LinearFEOperator : public LinearOperator 
    {
      public :

        /**
          Constructor. This copies the share pointer to the element operation from the input parameter object.
          */
        LinearFEOperator(const boost::shared_ptr<FEOperatorParameters>& params)
          : LinearOperator (params) 
        { 
          d_elemOp = (params->d_elemOp);
        }

        /**
          Destructor
          */
        virtual ~LinearFEOperator() { }

        /**
          This function will be called just before looping over the elements to 
          build the stiffness matrix, so if the derived classes need to perform 
          some initialization operations just before looping over the elements they can 
          do so by implementing these operations in this function.
          Also, the derived classes can access the parameters passed to the reset
          function by implementing this function.
          */
        virtual void preAssembly(const boost::shared_ptr<OperatorParameters>& )
        {
          //Implemented in derived classes. 
        }

        /**
          This function will be called just after looping over all the elements to
          build the stiffness matrix, so if the derived classes need to perform any
          operations such as freeing any temporary memory allocations they
          can do so by implementing these operations in this function.
          */
        virtual void postAssembly()
        {
          //Implemented in derived classes. 
        }

        /**
          This function will be called once for each element, just before performing
          the element operation. Ideally, the element operation should not deal with 
          global mesh related information such as DOFMap and global vectors and matrices. 
          This function typically extracts the local information from
          these global objects and passes them to the element operation.
          */
        virtual void preElementOperation(const AMP::Mesh::MeshElement &)
        {
          //Implemented in derived classes. 
        }

        /**
          This function will be called once for each element, just after performing the element operation.
          Typically, the element stiffness matrix is added to the global stiffness matrix in this function.
          */
        virtual void postElementOperation()
        {
          //Implemented in derived classes. 
        }

        /**
          This function creates the stiffness matrix and uses virtual
          function calls for setting values into the matrix.
          */
        void reset(const boost::shared_ptr<OperatorParameters>& );

      protected :

        boost::shared_ptr<ElementOperation> d_elemOp; /**< Shared pointer to the element operation */

      private :

    };

  }
}

#endif


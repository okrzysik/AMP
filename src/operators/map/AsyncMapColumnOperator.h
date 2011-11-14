#ifndef  included_AMP_AsyncMapColumnOperator
#define  included_AMP_AsyncMapColumnOperator

#include "operators/AsynchronousColumnOperator.h"
#include "operators/AsynchronousColumnOperatorParameters.h"

#include "ampmesh/Mesh.h"

namespace AMP {
namespace Operator {

  /** \brief   AsyncMapColumnOperator does not require parameters beyond those of AsynchronousColumnOperator
    */
  typedef AsynchronousColumnOperatorParameters  AsyncMapColumnOperatorParameters;


  /** \brief  A column of map operators used for Mesh to Mesh interactions through a mapping.
    * \details  Mappings often require overlapping communication to prevent serial threads
    * from arising during the computation.  Mappings often require results to be stored
    * in auxiliary or "frozen" vectors.  Therefore, the AsyncMapColumnOperator supports a
    * setVector method that will call setVector on all maps in the column.
    *
    * \see AsyncMapOperator
    */
  class  AsyncMapColumnOperator : public AsynchronousColumnOperator
  {
    public:
      /** \brief Constructor
        */
      AsyncMapColumnOperator ( const boost::shared_ptr<OperatorParameters> & params );

      /** \brief  Call setVector on all vectors in the column
        * \param[in] p  The auxiliary or "frozen" vector to store results in
        */
      void  setVector ( AMP::LinearAlgebra::Vector::shared_ptr &p );

      virtual void append ( boost::shared_ptr < Operator > op );


      /** \brief  A factory method.
        * \param[in]  manager  The mesh manager to search for maps to add to the column
        * \param[in]  db  The global database of the simulation
        * \tparam  MAP_TYPE  The type of map to create
        * \return  A column of map operators of type MAP_TYPE
        * \details  In order to use this builder, the MAP_TYPE must have 3 things:  a static member 
        * (or enum) called CommTagBase, a static method:  bool validMapType ( const std::string & ), and
        * a typedef for the parameters of the class called Parameters.
        *  
        */
      template <typename MAP_TYPE>
      static boost::shared_ptr<AsyncMapColumnOperator>  build ( AMP::Mesh::Mesh::shared_ptr mesh, boost::shared_ptr<Database> db );
  };


}
}

#include "AsyncMapColumnOperator.tmpl.h"

#endif

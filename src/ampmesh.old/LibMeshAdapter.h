#ifndef included_AMP_LibMeshAdapter
#define included_AMP_LibMeshAdapter

#include <boost/shared_ptr.hpp>
#include <map>
#include "utils/Database.h"
#include "vectors/Vector.h"
#include "vectors/VectorEntryMap.h"
#include "matrices/Matrix.h"

#include "mesh.h"
#include "mesh_data.h"
#include "equation_systems.h"

#include "MeshIterators.h"
#include "vectors/Variable.h"
#include "LibMeshPoint.h"
#include "LibMeshNode.h"
#include "LibMeshElement.h"


#include "LibMeshBoundarySet.h"
#include "LibMeshNodeElemMap.h"

#include "DOFMap.h"
#include "NodalRowMap.h"

namespace AMP { 
namespace Mesh {

  /** \class LibMeshAdapter
    * \brief  An adapter for libMesh that conforms to the AMP mesh interface.  DO NOT USE DIRECTLY, USE MeshManager::MeshAdapter.
    * \details  This is a lightweight wrapper class around the libMesh library.
    */
  class LibMeshAdapter
  {
    public:
      /** \brief  Convenience typedef
        */
      typedef   boost::shared_ptr<LibMeshAdapter>                                      shared_ptr;

      /** \brief  Underlying mesh data structure
        */
      typedef   ::Mesh                                                                 MeshDB;

      /** \brief  An element
        */
      typedef   LibMeshElement                                                         Element;

      /** \brief  A facet
        */
      typedef   LibMeshSide                                                            Facet;

      /** \brief  A node
        */
      typedef   LibMeshNode                                                            Node;

      /** \brief  A point
        */
      typedef   LibMeshPoint                                                           Point;

      /** \brief  Iterator for elements
        */
      typedef   MeshIteratorWrapper< ::Mesh::element_iterator,LibMeshElement>          ElementIterator;

      /** \brief  Iterator for elements
        */
      typedef   MeshIteratorWrapper< ::Mesh::const_element_iterator,LibMeshElement>    ConstElementIterator;

      /** \brief  Iterator for nodes
        */
      typedef   MeshIteratorWrapper< ::Mesh::node_iterator,LibMeshNode>                NodeIterator;

      /** \brief  Iterator for nodes
        */
      typedef   MeshIteratorWrapper< ::Mesh::const_node_iterator,LibMeshNode>          ConstNodeIterator;

      /** \brief  Iterator for nodes owned by this core
        */
      typedef   OwnedMeshIteratorWrapper< ::Mesh::node_iterator,LibMeshNode>                OwnedNodeIterator;

      /** \brief  Iterator for nodes on a boundary
        */
      typedef   MeshIteratorWrapper< LibMeshBoundarySet::NodeListIterator , LibMeshNode>     BoundaryNodeIterator;

      /** \brief  Iterator for facets on a boundary
        */
      typedef   MeshIteratorWrapper< LibMeshBoundarySet::SideListIterator , LibMeshElement>  BoundarySideIterator;

      /** \brief  Iterator for the nodes in an element in the "correct" order
        */
      typedef   MeshIteratorWrapper< LibMeshNodeElemMap::ElemListIterator , LibMeshElement>  NodeElementIterator;

      /** \brief  Iterator for nodes on a boundary owned by this core
        */
      typedef   OwnedMeshIteratorWrapper< LibMeshBoundarySet::NodeListIterator , LibMeshNode>  OwnedBoundaryNodeIterator;

      /** \brief  Iterator for the data registered with this mesh
        */
      typedef   std::map<std::string , AMP::LinearAlgebra::Vector::shared_ptr>::iterator   DataIterator;

      /** \brief  A mapping for computing sparsity structure in parallel
        */
      typedef   NodalRowMapTmpl<LibMeshAdapter>             NodalRowMap;

    private:
      boost::shared_ptr< ::Mesh>          d_libMesh;
      boost::shared_ptr< ::MeshData>      d_libMeshData;

      LibMeshBoundarySet               d_BoundarySet;
      LibMeshNodeElemMap               d_NodeElemMap;

      std::map < size_t , AMP::LinearAlgebra::Vector::shared_ptr >   d_vVectorCache;
      std::map < size_t , AMP::LinearAlgebra::Matrix::shared_ptr >   d_vMatrixCache;
      std::map < size_t , DOFMap::shared_ptr >   d_vDOFMapCache;

      std::map < std::string , AMP::LinearAlgebra::Vector::shared_ptr >   d_mData;

      void buildDOFMap ( AMP::LinearAlgebra::Variable::shared_ptr );
      void buildVector ( AMP::LinearAlgebra::Variable::shared_ptr );
      void buildMatrix ( AMP::LinearAlgebra::Variable::shared_ptr );
      void buildMatrix ( AMP::LinearAlgebra::Variable::shared_ptr , AMP::LinearAlgebra::Variable::shared_ptr );

      static LibMeshInit  *context;

      AMP::LinearAlgebra::ObjectSorter::shared_ptr  d_DefaultSorter;

      std::string     d_MeshName;
      bool            d_UseMeshNameInVariable;

      boost::shared_ptr<Database>    d_MeshDatabase;

    public:
      /** \brief Given a variable, this will return a variable specific to this mesh
        * \param[in]  var  The variable to associate with this mesh
        * \return  The variable of the same type as var associated with this mesh
        */
      AMP::LinearAlgebra::Variable::shared_ptr  getMeshSpecificVariable ( AMP::LinearAlgebra::Variable::shared_ptr var );

      /** \brief  Add the name of the mesh to the variable
        * \param[in]  var  The variable to add the name of this mesh to
        */
      void  appendMeshNameToVariable ( AMP::LinearAlgebra::Variable::shared_ptr var );

      /** \brief  Return the name of this mesh
        * \return  The name of the mesh
        */
      std::string     getMeshName ();

      /** \brief  Return a database to use with this mesh
        * \return The database
        */
      boost::shared_ptr<Database>  getDatabase () { return d_MeshDatabase; }

      /** \brief  Call makeConsistent on all registered AMP::LinearAlgebra::Vectors
        */
      void   makeDataConsistent ();

      /** \brief  Initialize libMesh
        * \param[in]  argc  The length of argv
        * \param[in]  argv  Parameters to initialize libMesh with
        * \param[in]  comm  The communicator to initialize libMesh on
        */
      static  void  initAdapter ( int argc , char **argv , AMP_MPI comm );

      /** \brief  Finalize libMesh
        */
      static  void  finalizeAdapter ();

      /** \brief  Set the name of the mesh
        * \param[in]  name  The new name of the mesh
        * \param[in]  t  Whether or not the mesh name should be used in variables associated
        * with this mesh
        */
      void setName ( const std::string &name , bool t = true );

      LibMeshAdapter ( const boost::shared_ptr< ::Mesh> &in );
      LibMeshAdapter ( boost::shared_ptr<Database>   d_MeshDatabase = boost::shared_ptr<Database> () );
      virtual ~LibMeshAdapter ();

      virtual void  readExodusIIFile ( std::string );
      virtual void  writeExodusIIFile ( std::string );

      virtual void  readIDEASFile ( std::string );
      virtual void  writeIDEASFile ( std::string );

      virtual void  generateCube ( size_t );

      void    registerVectorAsData ( AMP::LinearAlgebra::Vector::shared_ptr vec , const std::string & s = "");
      void    clearData ();

      DataIterator      beginData ();
      DataIterator      endData ();
      size_t            sizeData ();

      ElementIterator   beginElement();
      ElementIterator   endElement();

      ConstElementIterator  beginElement () const;
      ConstElementIterator  endElement ()   const;

      NodeElementIterator  beginElementForNode ( const Node &n );
      NodeElementIterator  endElementForNode ( const Node &n );

      Element  getElement ( size_t i );
      Node     getNode    ( size_t i );

      NodeIterator   beginNode();
      NodeIterator   endNode();

      NodeIterator   beginUsedNode ();
      NodeIterator   endUsedNode ();

      OwnedNodeIterator   beginOwnedNode();
      OwnedNodeIterator   endOwnedNode();

      ConstNodeIterator  beginNode () const;
      ConstNodeIterator  endNode ()   const;

      BoundarySideIterator  beginSideBoundary ( short int which );
      BoundarySideIterator  endSideBoundary ( short int which );

      BoundaryNodeIterator  beginBoundary ( short int which );
      BoundaryNodeIterator  endBoundary ( short int which );

      OwnedBoundaryNodeIterator  beginOwnedBoundary ( short int which );
      OwnedBoundaryNodeIterator  endOwnedBoundary ( short int which );

      Element  getElementOfFacet ( short int bid , Facet &f )
      {
        return  d_BoundarySet.getElementOfFacet ( bid , f.globalID () );
      }

      size_t    numLocalElements ();
      size_t    numLocalNodes ();

      size_t    numUsedNodes ();
      size_t    numUsedElements ();

      size_t    numTotalElements ();
      size_t    numTotalNodes ();

            ::Mesh  &getMesh();
      const ::Mesh  &getMesh() const;

      AMP_MPI      getComm();

      AMP::LinearAlgebra::Vector::shared_ptr   getPositionVector ( std::string name );

      AMP::LinearAlgebra::Vector::shared_ptr   createVector ( AMP::LinearAlgebra::Variable::shared_ptr variable );
      AMP::LinearAlgebra::Matrix::shared_ptr   createMatrix ( AMP::LinearAlgebra::Variable::shared_ptr operand , AMP::LinearAlgebra::Variable::shared_ptr result = AMP::LinearAlgebra::Variable::shared_ptr() );
      DOFMap::shared_ptr   getDOFMap ( AMP::LinearAlgebra::Variable::shared_ptr variable );


      int      dim ();

      std::vector<short int>     getBoundaryIds ( Node );
      bool                       isOnBoundary ( Node , short int );
      const std::set<short int>  &getBoundaryIds ();

      void     translate ( double x , double y , double z );
      void     scale  ( double );

      template <typename SIDE>
      LibMeshElement   getElementFromSide ( SIDE &in );

      AMP::LinearAlgebra::CommunicationList::shared_ptr  getNodalCommunicationList ();

      size_t   geometricDimension ();

      void   displaceMesh ( const AMP::LinearAlgebra::Vector::shared_ptr p );
  };

  template <typename SIDE>
  LibMeshElement   LibMeshAdapter::getElementFromSide ( SIDE &in )
  {
    LibMeshNode  root_node = getNode ( in.getNodeID ( 0 ) );
    LibMeshAdapter::NodeElementIterator  curElement = beginElementForNode ( root_node );
    while ( curElement != endElementForNode ( root_node ) )
    {
      bool found = true;
      for ( size_t j = 1 ; j != in.numNodes () ; j++ )
      {
        LibMeshAdapter::NodeElementIterator othElement = beginElementForNode ( getNode ( in.getNodeID ( j ) ) );
        bool  thisFound = false;
        while ( othElement != endElementForNode ( getNode ( in.getNodeID ( j ) ) ) )
        {
          if ( othElement->globalID() == curElement->globalID() )
            thisFound = true;
          othElement++;
        }
        found &= thisFound;
      }
      if ( found ) return *curElement;
      curElement++;
    }
    AMP_ERROR( "Can't find element" );
    return LibMeshElement();
  }


}
}

#include "LibMeshAdapter.inline.h"

#endif

#ifndef included_AMP_MeshManager
#define included_AMP_MeshManager

#include <vector>
#include <map>
#include <sstream>

#include <boost/shared_ptr.hpp>

#include "utils/AMP_MPI.h"
#include "utils/Utilities.h"
#include "utils/Database.h"

#include "vectors/MultiVariable.h"
#include "vectors/MultiVector.h"
#include "MeshVariable.h"


namespace AMP { 
namespace Mesh {

  namespace BoundaryIDs
  {
    enum { Inside = 1 , Outside = 2 , Top = 4 , Bottom = 8 };
  }

  class MeshManagerParameters : public ParameterBase
  {
    public:
      typedef  boost::shared_ptr<MeshManagerParameters>    shared_ptr;

      boost::shared_ptr<AMP::Database>  d_db;
      int     d_argc;
      char  **d_argv;
      MeshManagerParameters ( const boost::shared_ptr<AMP::Database>  & db );
     ~MeshManagerParameters ();
  };

  /**
    * \class MeshManagerTmpl
    * \brief A class used to abstract away mesh management from an application
    *
    * \details  This class provides routines for reading, accessing and writing meshes
    * from an input database.  Also, it handles initialization of the chosen
    * mesh database library.  The database format is of the following:
    *\code
    NumberOfMeshes = n
    Mesh_1 {
      Filename = "filename"
      MeshName = "name of mesh used in rest of file"
      MeshGroup = "name of mesh group"
      x_offset = x
      y_offset = y
      z_offset = z
      NumberOfElements = number of elements
    }
      .
      .
      .
    Mesh_n {
      Filename = "filename"
      MeshName = "name of mesh used in rest of file"
      MeshGroup = "name of mesh group"
      x_offset = x
      y_offset = y
      z_offset = z
      NumberOfElements = number of elements
    }
    \endcode
    */

  template <typename ADAPTER>
  class MeshManagerTmpl
  {
    public:
      /**
        *\typedef shared_ptr
        *\brief  Name for the shared pointer.
        *\details  Use this typedef for a reference counted pointer to a mesh manager object.
        */
      typedef          boost::shared_ptr<MeshManagerTmpl<ADAPTER> >    shared_ptr;

      enum DecompositionType { COARSE_DECOMP , DOMAIN_DECOMP };

      /**
        *\typedef Adapter
        *\brief  Name for the template parameter
        *\details  Use this typedef for the adapter to the mesh database.  This adapter
        * abstracts the mesh library away from amp code.
        */
      typedef          ADAPTER                                         Adapter;

      /**
        *\typedef AdapterPtr
        *\brief  Name for a shared pointer to an adapter
        *\details  Use this typedef for a reference counted pointer to the Adapter.  This
        *typedef is used to reduce the use of the <em>typename</em> keyword in the template
        *code
        */
      typedef typename ADAPTER::shared_ptr                             AdapterPtr;

      /**
        *\typedef MeshCollection
        *\brief Vector of mesh adapter shared pointers used to store meshes on this MPI communicator
        *\details  Depending on the type of parallelism used in the simulation, this collection may
        * be limited to one mesh or it may store many meshes.  For applications that require storing
        * several mesh communicators, such as contact, there may be more than one MeshCollection in
        * the manager.
        */
      typedef          std::vector< AdapterPtr >                       MeshCollection;

      /**
        *\typedef MeshIterator
        *\brief Name of the iterator over a mesh collection
        *\details A convenience typedef to reduce the need for <em>typename</em> in the template
        * implemenation
        */
      typedef typename std::vector< AdapterPtr >::iterator             MeshIterator;

      /**
        *\typedef ConstMeshIterator
        *\brief Name of the iterator over a mesh collection
        *\details A convenience typedef to reduce the need for <em>typename</em> in the template
        * implemenation
        */
      typedef typename std::vector< AdapterPtr >::const_iterator       ConstMeshIterator;

      /**
        *\typedef MeshNameIterator
        *\brief A type that maps mesh names to indices into a mesh collection.
        *\details A type that maps mesh names to indices into a mesh collection.
        */
      typedef          std::map<std::string , size_t>::iterator        MeshNameIterator;

      enum  MapDominance { Master , Slave };
      enum  MapConstructionParam { Synchronous , Asynchronous };

    private:
      std::map<std::string , size_t>   d_mMeshNameLookup;
      MeshCollection                   d_vMeshes;
      unsigned int                     d_iIterationCount;
      DecompositionType                d_DecompositionType;

      std::vector<boost::shared_ptr<Database> >  d_MeshDatabases;

      std::vector<AMP_MPI>                       d_MapComms;
      std::vector<std::string>                   d_MapMeshName;
      std::vector<short int>                     d_MapBoundaryId;
      std::vector<boost::shared_ptr<Database> >  d_MapDBs;
      std::vector<MapDominance>                  d_MapDominance;
      std::vector<MapConstructionParam>          d_MapConstructionParam;
      std::vector<std::string>                   d_MapType;
      std::vector<size_t>                        d_MapTagOffset;



      std::string  getMeshNameInDB ( int i , boost::shared_ptr<Database> & );
      int   getNumberOfMeshes ( boost::shared_ptr<Database> &db );
      void  readMesh ( int i , boost::shared_ptr<Database> & );
      std::string      d_MyMeshName;
      size_t           d_MyMeshIndex;
      size_t           d_FirstMeshIndex;
      size_t           d_NumMeshes;
      AdapterPtr       d_MyMesh;
      AMP_MPI          d_DecompComm;

      void  computeRoughDecomposition ( std::vector<size_t> &procs_per_mesh , boost::shared_ptr<Database> &db );
      void  computeMeshIndex ( std::vector<size_t> & );
      void  decompSetup ( const MeshManagerParameters::shared_ptr & params );
      void  nondecompSetup ( const MeshManagerParameters::shared_ptr & params );
      void  makeDataConsistent ();
      AMP::LinearAlgebra::Vector::shared_ptr  backwardsCompatibleCreateVector ( AMP::LinearAlgebra::Variable::shared_ptr  in_var );
      void  setupMeshToMesh ( const MeshManagerParameters::shared_ptr & );
      void  buildMeshToMesh ( boost::shared_ptr<Database> , size_t);

    protected:
      MeshManagerTmpl ( const MeshManagerTmpl & );

      /**
        * \brief Destroy the mesh database environment.
        * \details  Even though this is the 21st century, some people believe that their work is so important
        * the work of others.  This function cleans up these global variables.
        */
      void finalize();

      /**
        * \brief Free the reference counted pointers
        * \details  Sets internal pointers to 0 and resizes any MeshCollection to 0
        */
      void   free ();

      boost::shared_ptr<Database>  d_ThisCommDB;

    public:
      //!  @name Constructors/destructors

      /**
        * \fn MeshManagerTmpl ( const MeshManagerParameters::shared_ptr &params )
        * \param params Parameters for constructing a mesh manager.
        * \brief Read in mesh files, partition domain, and prepare environment for simulation
        * \details  For trivial parallelsim, this method reads in the meshes on each processor.  Each
        * processor contains a piece of each mesh.  For massive parallelism, each mesh is on its own
        * communicator.  As such, some math libraries must be initialized accordingly.
        */
      MeshManagerTmpl ( const MeshManagerParameters::shared_ptr &params );

      /**
        * \fn ~MeshManagerTmpl ()
        * \brief Alias for free ()
        * \details
        */
     ~MeshManagerTmpl ();


      //!   @name Mesh and data access
      /**
       * \return Iterator pointing to first mesh in the main mesh collection.
       * \brief  Returns the first mesh in the main mesh collection
       * \details In trivial parallelism, the mesh collection is a set of meshes.  This
       * will return a pointer to the reference counted pointer of the first mest in the
       * collection.  For massive parallelism, this will return the mesh on this processor.
       */
      MeshIterator       beginMeshes ();

      /**
        * \return Iterator point past the last mesh in the main mesh collection.
        * \brief  Returns a pointer past the end of the main mesh collection
        * \details As per the iterator idiom, this returns an iterator on past
        * the end of the main mesh collection
        */
      MeshIterator       endMeshes();

      /**
       * \return Iterator pointing to first mesh in the main mesh collection.
       * \brief  Returns the first mesh in the main mesh collection
       * \details In trivial parallelism, the mesh collection is a set of meshes.  This
       * will return a pointer to the reference counted pointer of the first mest in the
       * collection.  For massive parallelism, this will return the mesh on this processor.
       */
      ConstMeshIterator  beginMeshes () const;

      /**
        * \return Iterator point past the last mesh in the main mesh collection.
        * \brief  Returns a pointer past the end of the main mesh collection
        * \details As per the iterator idiom, this returns an iterator on past
        * the end of the main mesh collection
        */
      ConstMeshIterator  endMeshes () const;

      /**
        * \return Iterator pointing to the name of the first mesh
        * \brief An iterable list of names of meshes in the same order of the meshes.
        */
      MeshNameIterator   beginMeshNames();

      /**
        * \return Iterator pointing past the name of the last mesh
        * \brief An iterable list of names of meshes in the same order of the meshes.
        */
      MeshNameIterator   endMeshNames();

      /**
        * \return Number of meshes on this communicator
        * \brief  Returns the number of meshes on this communicator
        */
      size_t             numMeshes() { return d_vMeshes.size(); }

      /**
        * \return  A pointer to the mesh
        * \brief   Return the mesh named mesh_name, null if the mesh does not exist
        * \param[in]  mesh_name  The name of the mesh to return
        */
      typename Adapter::shared_ptr   getMesh ( const std::string &mesh_name );

      /**
        * \return  The mesh on this communicator
        * \brief  If the communicator has but one mesh on it, this method will return it
        * \deprecated  Use the beginMeshes() interface
        */
      typename Adapter::shared_ptr   getMesh ();

      /**
        * \return  A Vector
        * \brief  Create a vector based on variable
        * \param[in]  variable  The variable that describes the vector to be created
        */
      AMP::LinearAlgebra::Vector::shared_ptr    createVector ( AMP::LinearAlgebra::Variable::shared_ptr variable );

      /**
        * \return A Nodal3Vector that stores the positions of the nodes
        * \brief  Create a vector that stores the positions of all nodes on this communicator
        * \param[in]  in  The name to give the vector
        */
      AMP::LinearAlgebra::Vector::shared_ptr    createPositionVector ( std::string in );

      /**
        * \brief  When reading or writing the meshes in this manager, also read/write the vector v
        * \param[in] v  The vector to read/write
        * \param[in] s  The name to give the vector if different from the variable name
        */
      void                  registerVectorAsData ( AMP::LinearAlgebra::Vector::shared_ptr v , const std::string &s = "" );

      AMP_MPI         getMeshComm () { return d_DecompComm; }


      size_t          getNumMaps () const { return d_MapComms.size(); }
      AMP_MPI         getMapComm ( size_t i )  const{ return d_MapComms[i]; }
      std::string     getMapMeshName ( size_t i ) const { return d_MapMeshName[i]; }
      short int       getMapBoundaryId ( size_t i ) const { return d_MapBoundaryId[i]; }
      boost::shared_ptr<Database>  getMapDB ( size_t i ) const { return d_MapDBs[i]; }
      MapDominance    getMapDominance ( size_t i ) const { return d_MapDominance[i]; }
      MapConstructionParam  getMapConstructionParam ( size_t i ) const { return d_MapConstructionParam[i]; }
      const std::string    &getMapType ( size_t i ) const { return d_MapType[i]; }
      size_t          getMapTagOffset ( size_t i ) const { return d_MapTagOffset[i]; }

      //!   @name File I/O and mesh creation

      template <template<typename MANAGER> class IO>
      void writeFile ( const std::string &fname , size_t iteration_count );

      template <template<typename MANAGER> class IO>
      void updateFile ( const char *fname , double time );

      template <template<typename MANAGER> class IO>
      void readFile ( const std::string &fname , size_t iteration_count );

      typename Adapter::shared_ptr   readMeshFromExodusFile ( std::string fname , std::string , boost::shared_ptr<Database> );
      void                           addMesh ( AdapterPtr , const std::string & );
      void     clearData ();
      typename Adapter::shared_ptr   generateCube ( size_t numNodesPerSide , std::string );

      boost::shared_ptr<Database>  getMeshDatabase ( typename Adapter::shared_ptr mesh ) { return getMeshDatabase ( mesh->getMeshName() ); }
      boost::shared_ptr<Database>  getMeshDatabase ( size_t i ) { return d_MeshDatabases[i]; }
      boost::shared_ptr<Database>  getMeshDatabase ( const std::string &s ) { return d_MeshDatabases[d_mMeshNameLookup[s]]; }

      boost::shared_ptr<Database>  getThisCommDB () { return d_ThisCommDB; }
      AMP::LinearAlgebra::Variable::shared_ptr  getMeshVariable ( AMP::LinearAlgebra::Variable::shared_ptr  var , const std::string &mesh_name );
  };

}
}

#include "MeshManager.tmpl.h"

#include "MeshType.h"
#include "LibMeshAdapter.h"
#include "MeshVectorSelector.h"
namespace AMP { 
namespace Mesh {

  typedef  MeshManagerTmpl<MESH_TYPE>     MeshManager;
  typedef  VS_ByMeshNameTmpl<MESH_TYPE>   VS_ByMeshName;
  typedef  VS_ByMeshTmpl<MESH_TYPE>       VS_ByMesh;

}
}

#endif

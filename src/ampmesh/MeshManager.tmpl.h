
#include <vector>
#include <map>
#include <sstream>
#include <string.h>

#include <boost/shared_ptr.hpp>

#include "utils/Utilities.h"

#include "utils/AMP_MPI.h"

#include "vectors/MultiVariable.h"
#include "vectors/MultiVector.h"
#include "vectors/NullVector.h"
#include "vectors/CommCollectVector.h"

#include "MeshCollectionVector.h"


namespace AMP { 
namespace Mesh {

  template <typename ADAPTER>
  void MeshManagerTmpl<ADAPTER>::free ()
  {
    d_MyMesh = AdapterPtr();
    d_vMeshes.resize ( 0 );
    finalize ();
  }

  template <typename ADAPTER>
  MeshManagerTmpl<ADAPTER>::~MeshManagerTmpl ()
  {
    free ();
  }

  template <typename ADAPTER>
  std::string  MeshManagerTmpl<ADAPTER>::getMeshNameInDB ( int i , boost::shared_ptr<Database> &db )
  {
    char id[1000];
    sprintf(id, "Mesh_%d", i);

    AMP_ASSERT ( db->keyExists ( id ) );
    AMP_ASSERT ( db->isDatabase ( id ) );
    return std::string (id);
  }

  template <typename ADAPTER>
  int MeshManagerTmpl<ADAPTER>::getNumberOfMeshes ( boost::shared_ptr<Database> &db )
  {
    AMP_ASSERT ( db->keyExists ( "NumberOfMeshes" ) );
    AMP_ASSERT ( db->isInteger ( "NumberOfMeshes" ) );
    return db->getInteger ( "NumberOfMeshes" );
  }

  template <typename ADAPTER>
  void MeshManagerTmpl<ADAPTER>::computeRoughDecomposition ( std::vector<size_t> &first_node , boost::shared_ptr<Database>  &db )
  {
    int numToRead = getNumberOfMeshes ( db );
    std::vector<size_t>  mesh_sizes ( numToRead );
    std::vector<size_t>  num_procs ( numToRead );
    first_node.resize ( numToRead );
    int total_mesh_size = 0;
    AMP_MPI globalComm(AMP_COMM_WORLD);
    int size = globalComm.getSize();

    if ( numToRead == 0 ) return;

    for ( int i = 1 ; i <= numToRead ; i++ )
    {
      std::string  mesh_name = getMeshNameInDB ( i , db );
      boost::shared_ptr<Database>  curMeshDB = db->getDatabase ( mesh_name );
      AMP_ASSERT ( curMeshDB->keyExists ( "NumberOfElements" ) );
      AMP_ASSERT ( curMeshDB->isInteger ( "NumberOfElements" ) );
      total_mesh_size += mesh_sizes[i-1] = curMeshDB->getInteger ( "NumberOfElements" );
      num_procs[i-1] = 1;
    }
    int  part_size = total_mesh_size / size;
    int  left_to_find = size - numToRead;
    while ( left_to_find > 0 )
    {
      left_to_find--;
      std::vector<size_t>::iterator  which_mesh = std::max_element ( mesh_sizes.begin() , mesh_sizes.end() );
      *which_mesh -= part_size;
      int  mesh_id = which_mesh - mesh_sizes.begin();
      num_procs[mesh_id]++;
    }
    size_t cur_left = 0;
    for ( int i = 0 ; i != numToRead ; i++ )
    {
      first_node[i] = cur_left;
      cur_left += num_procs[i];
    }
  }

  template <typename ADAPTER>
  void MeshManagerTmpl<ADAPTER>::computeMeshIndex ( std::vector<size_t> &partition )
  {
    d_MyMeshIndex = partition.size() - ( 1 - d_FirstMeshIndex );
    // This loop has crazy fencepost issues since meshes in files maybe numbered from 1
    AMP_MPI globalComm(AMP_COMM_WORLD);
    int rank = globalComm.getRank();
    for ( size_t i = 1 ; i != partition.size() ; i++ )
    {
      if ( (int)partition[i] > rank )
      {
        d_MyMeshIndex = i - (1 - d_FirstMeshIndex);
        break;
      }
    }
  }


  template <typename ADAPTER>
  void MeshManagerTmpl<ADAPTER>::decompSetup ( const MeshManagerParameters::shared_ptr &params )
  {

    int numToRead = 0;
    if ( params->d_db->keyExists ( "NumberOfMeshes" ) )
    {
      AMP_ASSERT ( params->d_db->isInteger ( "NumberOfMeshes" ) );
      numToRead = params->d_db->getInteger ( "NumberOfMeshes" );
    }

    d_FirstMeshIndex = params->d_db->getIntegerWithDefault ( "FirstMeshIndex" , 1 );
    
    AMP_MPI globalComm(AMP_COMM_WORLD);
    int comm_size = globalComm.getSize();
    bool doDecompSetup = false;
    if ( params->d_db->keyExists ( "DomainDecomposition" ) )
    {
      doDecompSetup = params->d_db->getInteger ( "DomainDecomposition" ) == 1 ? true : false;
    }
    if ( doDecompSetup )
    {
      AMP_ASSERT ( numToRead <= comm_size );
      d_DecompositionType = DOMAIN_DECOMP;
      std::vector<size_t>  decomposition;
      computeRoughDecomposition ( decomposition , params->d_db );
      computeMeshIndex ( decomposition );
      int  rank = globalComm.getRank();
      d_DecompComm = globalComm.split(d_MyMeshIndex,rank);
      ADAPTER::initAdapter ( params->d_argc , params->d_argv , d_DecompComm );
      readMesh ( d_MyMeshIndex , params->d_db );
      std::string MeshNameInDB = getMeshNameInDB ( d_MyMeshIndex , params->d_db );
      if ( params->d_db->keyExists ( MeshNameInDB ) )
        {
	  if(params->d_db->getDatabase ( MeshNameInDB )->keyExists ( "DatabaseName" ))
	    {
	      d_ThisCommDB = params->d_db->getDatabase ( params->d_db->getDatabase ( MeshNameInDB )->getString ( "DatabaseName" ) );
	    }
	  else
	    {
	      d_ThisCommDB = params->d_db;
	    }
        }
      else
        {
        d_ThisCommDB = params->d_db;
        }
    }
    else
    {
      d_DecompositionType = COARSE_DECOMP;
      nondecompSetup ( params );
      d_MyMesh = AdapterPtr ();
      d_MyMeshName = "";
      if ( numToRead == 1 )
      {
        d_MyMesh = *beginMeshes();
       // std::string MeshNameInDB = getMeshNameInDB ( d_MyMeshIndex , params->d_db );
       // d_MyMeshName = params->d_db->getDatabase ( MeshNameInDB )->getString ( "MeshName" );
      }
    }

  }

  template <typename ADAPTER>
  void MeshManagerTmpl<ADAPTER>::readMesh ( int i , boost::shared_ptr<Database> &db )
  {
    std::string mesh_name = getMeshNameInDB ( i , db );
    boost::shared_ptr<Database>  curMeshDB = db->getDatabase ( mesh_name );
    boost::shared_ptr<Database> meshOpDb;
    if ( curMeshDB->keyExists ( "DatabaseName" ) )
    {
      meshOpDb = db->getDatabase ( curMeshDB->getString ( "DatabaseName" ) );
    }

    AMP_ASSERT ( curMeshDB->keyExists ( "Filename" ) );
    AMP_ASSERT ( curMeshDB->isString ( "Filename" ) );
    AMP_ASSERT ( curMeshDB->keyExists ( "MeshName" ) );
    AMP_ASSERT ( curMeshDB->isString ( "MeshName" ) );
    d_MyMeshName = curMeshDB->getString ( "MeshName" );
    d_MyMesh = readMeshFromExodusFile ( curMeshDB->getString ( "Filename" ) ,
                                        d_MyMeshName , meshOpDb );
    d_MyMesh->setName ( d_MyMeshName , true );

    AMP_ASSERT ( curMeshDB->keyExists ( "x_offset" ) );
    AMP_ASSERT ( curMeshDB->keyExists ( "y_offset" ) );
    AMP_ASSERT ( curMeshDB->keyExists ( "z_offset" ) );
    AMP_ASSERT ( curMeshDB->isDouble ( "x_offset" ) );
    AMP_ASSERT ( curMeshDB->isDouble ( "y_offset" ) );
    AMP_ASSERT ( curMeshDB->isDouble ( "z_offset" ) );
    getMesh ( d_MyMeshName )->translate ( curMeshDB->getDouble ( "x_offset" ) ,
                                          curMeshDB->getDouble ( "y_offset" ) ,
                                          curMeshDB->getDouble ( "z_offset" ) );
  }

  template <typename ADAPTER>
  void MeshManagerTmpl<ADAPTER>::nondecompSetup ( const MeshManagerParameters::shared_ptr &params )
  {
    AMP_MPI globalComm(AMP_COMM_WORLD);
    d_DecompComm = globalComm.dup();
    ADAPTER::initAdapter ( params->d_argc , params->d_argv , d_DecompComm );
    int numToRead = getNumberOfMeshes ( params->d_db );

    for ( int i = 1 ; i <= numToRead ; i++ )
    {
      readMesh ( i , params->d_db );
    }
    d_ThisCommDB = params->d_db;
    d_DecompositionType = COARSE_DECOMP;
    if ( numToRead > 1 )
      d_MyMesh.reset();
  }

  template <typename ADAPTER>
  void MeshManagerTmpl<ADAPTER>::finalize ()
  {
    ADAPTER::finalizeAdapter ();
  }

  template <typename ADAPTER>
  MeshManagerTmpl<ADAPTER>::MeshManagerTmpl ( const MeshManagerParameters::shared_ptr & params )
  {
    decompSetup ( params );
    setupMeshToMesh ( params );
  }

  template <typename ADAPTER>
  void MeshManagerTmpl<ADAPTER>::buildMeshToMesh ( boost::shared_ptr<Database> db , size_t mapID )
  {
	  //std::string mname1 = db->getString ( "Mesh1" );
	  //std::string mname2 = db->getString ( "Mesh2" );
	  //short int  surfID1 = db->getInteger ( "Surface1" );
	  //short int  surfID2 = db->getInteger ( "Surface2" );

	  int inComm = -1;
	  if ( getMesh ( db->getString ( "Mesh1")) || getMesh ( db->getString ( "Mesh2" )))
	  {
		  inComm = 1;
	  }
      AMP_MPI globalComm(AMP_COMM_WORLD);
      int rank = globalComm.getRank();
      AMP_MPI newComm = globalComm.split(inComm,rank);
	  if ( getMesh ( db->getString ( "Mesh1" )))
	  {
		  d_MapBoundaryId.push_back ( db->getInteger ( "Surface1" ));
		  d_MapComms.push_back ( newComm );
		  if ( getMesh ( db->getString ( "Mesh2")))
		  {
			  d_MapConstructionParam.push_back ( Asynchronous );
		  }
		  else
		  {
			  d_MapConstructionParam.push_back ( Synchronous );
		  }
		  d_MapDBs.push_back ( db );
		  d_MapDominance.push_back ( Master );
		  d_MapMeshName.push_back ( db->getString ( "Mesh1" ) );
		  d_MapType.push_back ( db->getString ( "MapType" ) );
      d_MapTagOffset.push_back ( mapID );
	  }
	  if ( getMesh ( db->getString ( "Mesh2" )))
	  {
		  d_MapBoundaryId.push_back ( db->getInteger ( "Surface2" ));
		  d_MapComms.push_back ( newComm );
		  if ( getMesh ( db->getString ( "Mesh1")))
		  {
			  d_MapConstructionParam.push_back ( Asynchronous );
		  }
		  else
		  {
			  d_MapConstructionParam.push_back ( Synchronous );
		  }
		  d_MapDBs.push_back ( db );
		  d_MapDominance.push_back ( Slave );
		  d_MapMeshName.push_back ( db->getString ( "Mesh2" ) );
		  d_MapType.push_back ( db->getString ( "MapType" ) );
      d_MapTagOffset.push_back ( mapID );
	  }
  }

  template <typename ADAPTER>
  void MeshManagerTmpl<ADAPTER>::setupMeshToMesh ( const MeshManagerParameters::shared_ptr &params )
  {
    if ( params->d_db->keyExists ( "NumberOfMeshToMeshMaps" ) )
    {
      AMP_ASSERT ( params->d_db->isInteger ( "NumberOfMeshToMeshMaps" ) );
      int numMapsToSearch = params->d_db->getInteger ( "NumberOfMeshToMeshMaps" );
      for ( int i = 1 ; i <= numMapsToSearch ; i++ )
      {
        std::stringstream dbName;
        dbName << "MeshToMeshMap_" << i;
        buildMeshToMesh ( params->d_db->getDatabase ( dbName.str() ) , 2*i );
      }
    }
  }

  template <typename ADAPTER>
  typename MeshManagerTmpl<ADAPTER>::MeshIterator  MeshManagerTmpl<ADAPTER>::beginMeshes ()
  {
    return d_vMeshes.begin();
  }

  template <typename ADAPTER>
  typename MeshManagerTmpl<ADAPTER>::MeshIterator  MeshManagerTmpl<ADAPTER>::endMeshes()
  {
    return d_vMeshes.end();
  }

  template <typename ADAPTER>
  typename MeshManagerTmpl<ADAPTER>::ConstMeshIterator  MeshManagerTmpl<ADAPTER>::beginMeshes () const
  {
    return d_vMeshes.begin();
  }

  template <typename ADAPTER>
  typename MeshManagerTmpl<ADAPTER>::ConstMeshIterator  MeshManagerTmpl<ADAPTER>::endMeshes () const
  {
    return d_vMeshes.end();
  }

  template <typename ADAPTER>
  typename MeshManagerTmpl<ADAPTER>::MeshNameIterator   MeshManagerTmpl<ADAPTER>::beginMeshNames()
  {
    return d_mMeshNameLookup.begin();
  }

  template <typename ADAPTER>
  typename MeshManagerTmpl<ADAPTER>::MeshNameIterator   MeshManagerTmpl<ADAPTER>::endMeshNames()
  {
    return d_mMeshNameLookup.end();
  }

  template <typename ADAPTER>
  AMP::LinearAlgebra::Variable::shared_ptr  MeshManagerTmpl<ADAPTER>::getMeshVariable ( AMP::LinearAlgebra::Variable::shared_ptr  var , const std::string &mesh_name )
  {
    std::stringstream  new_name;
    new_name << var->getName() << "_" << mesh_name;
    return var->cloneVariable ( new_name.str() );
  }

  template <typename ADAPTER>
  AMP::LinearAlgebra::Vector::shared_ptr MeshManagerTmpl<ADAPTER>::backwardsCompatibleCreateVector ( AMP::LinearAlgebra::Variable::shared_ptr in_var )
  {
    AMP_MPI globalComm(AMP_COMM_WORLD);
    AMP::LinearAlgebra::Vector::shared_ptr ret_val;
    if ( d_MyMesh )  // Check for dcomp
    {
//      AMP::LinearAlgebra::Variable::shared_ptr  new_var = getMeshVariable ( in_var , d_MyMeshName );
//      ret_val = MultiMeshVector::view ( d_MyMesh->createVector ( in_var ) , globalComm );
      ret_val = d_MyMesh->createVector ( in_var );
    }
    else
    {
      AMP::LinearAlgebra::Variable::shared_ptr  new_var ( new AMP::LinearAlgebra::MultiVariable ( in_var->getName() ) );
      ret_val = AMP::LinearAlgebra::Vector::shared_ptr ( AMP::LinearAlgebra::MultiVector::create ( new_var , globalComm ) );
      MeshIterator curMesh = beginMeshes ();
      while ( curMesh != endMeshes () )
      {
        AMP::LinearAlgebra::Variable::shared_ptr meshVar = in_var->cloneVariable();
        new_var->castTo<AMP::LinearAlgebra::MultiVariable>().add ( meshVar->cloneVariable() );
        meshVar->castTo<MeshVariable<Adapter> >().setMesh ( *curMesh );
        ret_val->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( (*curMesh)->createVector ( meshVar ) );
        curMesh++;
      }
    }
//    std::cout << "WARNING: using a soon-to-be-deprecated createVector case" << std::endl;
    return ret_val;
  }

  template <typename ADAPTER>
  AMP::LinearAlgebra::Vector::shared_ptr MeshManagerTmpl<ADAPTER>::createVector ( AMP::LinearAlgebra::Variable::shared_ptr in_var )
  {
    AMP_MPI globalComm(AMP_COMM_WORLD);
    AMP::LinearAlgebra::Vector::shared_ptr retVal;
    if ( in_var->isA<AMP::LinearAlgebra::MultiVariable> () )
    {
      AMP::LinearAlgebra::MultiVariable  &mvar = in_var->castTo<AMP::LinearAlgebra::MultiVariable>();
      retVal = AMP::LinearAlgebra::MultiVector::create ( in_var , getMeshComm() );
      AMP::LinearAlgebra::MultiVariable::iterator  curVar = mvar.beginVariable();
      while ( curVar != mvar.endVariable() )
      {
        retVal->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( createVector ( *curVar ) );
        curVar++;
      }
    }
    else if ( in_var->isA<MeshVariable<Adapter> > () )
    {
      MeshVariable<Adapter> &var = in_var->castTo<MeshVariable<Adapter> >();
      if ( var.getMesh() )
      {
        retVal = var.getMesh()->createVector ( in_var );
      }
      else
      {
        retVal = MeshCollectionVector::create<ADAPTER> ( in_var , globalComm );
        MeshIterator curMesh = beginMeshes ();
        while ( curMesh != endMeshes () )
        {
          AMP::LinearAlgebra::Variable::shared_ptr meshVar = in_var->cloneVariable();
          meshVar->castTo<MeshVariable<Adapter> >().setMesh ( *curMesh );
          retVal->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( (*curMesh)->createVector ( meshVar ) );
          curMesh++;
        }
      }
    }
    else
    {
      AMP_ERROR( "Unsupported variable type" );
    }
    return retVal;
  }

  template <typename ADAPTER>
  AMP::LinearAlgebra::Vector::shared_ptr  MeshManagerTmpl<ADAPTER>::createPositionVector ( std::string name )
  {
    AMP_MPI globalComm(AMP_COMM_WORLD);
    AMP::LinearAlgebra::Vector::shared_ptr    ret_val;
    if ( d_MyMesh )
    {
//      std::stringstream new_name;
//      new_name << name << "_" << d_MyMeshName;
//      ret_val = MultiMeshVector::view ( d_MyMesh->getPositionVector ( name ), globalComm );
    ret_val = d_MyMesh->getPositionVector ( name );
    }
    else
    {
      AMP::LinearAlgebra::Variable::shared_ptr  new_var ( new AMP::LinearAlgebra::MultiVariable ( name ) );
      ret_val = AMP::LinearAlgebra::Vector::shared_ptr ( AMP::LinearAlgebra::MultiVector::create ( new_var , globalComm ) );
      MeshIterator curMesh = beginMeshes ();
      while ( curMesh != endMeshes () )
      {
        std::stringstream  new_name;
        new_name << name << "_" << (*curMesh)->getMeshName();
        AMP::LinearAlgebra::Vector::shared_ptr  subvector ( (*curMesh)->getPositionVector ( new_name.str() ) );
        new_var->castTo<AMP::LinearAlgebra::MultiVariable>().add ( subvector->getVariable () );
        ret_val->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( subvector );
        curMesh++;
      }
    }
    return ret_val;
  }

  template <typename ADAPTER>
  void MeshManagerTmpl<ADAPTER>::registerVectorAsData ( AMP::LinearAlgebra::Vector::shared_ptr vec , const std::string &name )
  {
    if ( vec->isA<AMP::LinearAlgebra::CommCollectVector>() )
    {
      vec = vec->castTo<AMP::LinearAlgebra::CommCollectVector>().getSmallCommVector();
    }

    if ( vec->getVariable()->isA<MeshVariable<ADAPTER> >() )
    {
      if ( vec->getVariable()->castTo<MeshVariable<ADAPTER> >().getMesh() )
      {
        vec->getVariable()->castTo<MeshVariable<ADAPTER> >().getMesh()->registerVectorAsData ( vec , name );
        return;
      }
    }

    if ( ( d_DecompositionType == COARSE_DECOMP ) && d_MyMesh )
    {
      d_MyMesh->registerVectorAsData ( vec , name );
    }
    else
    {
      AMP::LinearAlgebra::MultiVector &mv = vec->castTo<AMP::LinearAlgebra::MultiVector> ();
      AMP::LinearAlgebra::MultiVector::vector_iterator  curVector = mv.beginVector();
      while ( curVector != mv.endVector() )
      {
        if ( (*curVector)->getVariable()->isA<MeshVariable<Adapter> > () )
        {
          MeshVariable<Adapter>  &curVar = (*curVector)->getVariable()->castTo<MeshVariable<Adapter> >();
          if ( curVar.getMesh() )
            curVar.getMesh()->registerVectorAsData ( (*curVector) , name );
        }
        curVector++;
      }
    }
  }

  template <typename ADAPTER>
  void MeshManagerTmpl<ADAPTER>::makeDataConsistent ()
  {
    MeshIterator curMesh = d_vMeshes.begin();
    while ( curMesh != d_vMeshes.end() )
    {
      (*curMesh)->makeDataConsistent();
      curMesh++;
    }
  }

  template <typename ADAPTER>
  template <template<typename MANAGER> class IO>
  void MeshManagerTmpl<ADAPTER>::updateFile ( const char *fname , double time )
  {
    DEPRECATED("updateFile","writeFile");

    IO<MeshManagerTmpl<ADAPTER> >  writer;
    writer.updateFile ( *this , fname , time , d_iIterationCount );
    d_iIterationCount++;
  }

  template <typename ADAPTER>
  template <template<typename MANAGER> class IO>
  void MeshManagerTmpl<ADAPTER>::writeFile ( const std::string &fname , size_t iteration_count )
  {
    makeDataConsistent ();
    std::stringstream name, mastername;
    IO<MeshManagerTmpl<ADAPTER> >  writer;

    if ( d_MyMesh )
    {
      name << fname << "_" << d_MyMeshName << "_" << iteration_count << "." << writer.getExtension ();
      mastername << fname << "_" << iteration_count << "." << writer.getExtension ();
      writer.decompWriteFile ( *this , *d_MyMesh , name.str() , mastername.str() );
    }
    else
    {
      name << fname << "_" << iteration_count << "." << writer.getExtension ();
      writer.writeFile ( *this , name.str() );
    }
  }

  template <typename ADAPTER>
  template <template<typename MANAGER> class IO>
  void MeshManagerTmpl<ADAPTER>::readFile ( const std::string &fname , size_t iteration_count )
  {
    std::stringstream name, mastername;
    IO<MeshManagerTmpl<ADAPTER> >  writer;

    if ( d_MyMesh )
    {
      name << fname << "_" << d_MyMeshName << "_" << iteration_count << "." << writer.getExtension ();
      mastername << fname << "_" << iteration_count << "." << writer.getExtension ();
      writer.decompReadFile ( *this , *d_MyMesh , name.str() , mastername.str() );
    }
    else
    {
      name << fname << "_" << iteration_count << "." << writer.getExtension ();
      writer.readFile ( *this , name.str() );
    }
    makeDataConsistent ();
  }

  template <typename ADAPTER>
  typename ADAPTER::shared_ptr   MeshManagerTmpl<ADAPTER>::getMesh ( const std::string &mesh_name )
  {
    std::map<std::string , size_t>::iterator  mesh = d_mMeshNameLookup.find ( mesh_name );
    if ( mesh == d_mMeshNameLookup.end() )
    {
      return typename ADAPTER::shared_ptr ();
    }
    return d_vMeshes[mesh->second];
  }

  template <typename ADAPTER>
  typename ADAPTER::shared_ptr   MeshManagerTmpl<ADAPTER>::getMesh ()
  {
    AMP_ASSERT ( d_MyMesh );
    return d_MyMesh;
  }

  template <typename ADAPTER>
  void  MeshManagerTmpl<ADAPTER>::addMesh ( AdapterPtr ptr , const std::string &s )
  {
    d_vMeshes.push_back ( ptr );
    d_mMeshNameLookup[s] = d_vMeshes.size();
    d_MyMesh = ptr;
  }

  template <typename ADAPTER>
  typename ADAPTER::shared_ptr   MeshManagerTmpl<ADAPTER>::readMeshFromExodusFile ( std::string fname , std::string mesh_name , boost::shared_ptr<Database> db )
  {
    typename ADAPTER::shared_ptr  retVal ( new ADAPTER ( db ) );
    d_vMeshes.push_back ( retVal );

    if ( d_mMeshNameLookup.find ( mesh_name ) != d_mMeshNameLookup.end() )
    {
      AMP_ERROR( "Mesh already exists" );
    }
    d_mMeshNameLookup[mesh_name] = d_vMeshes.size()-1;

    retVal->readExodusIIFile ( fname );
    return retVal;
  }


  template <typename ADAPTER>
  typename ADAPTER::shared_ptr   MeshManagerTmpl<ADAPTER>::generateCube ( size_t numNodesPerSide , std::string mesh_name )
  {
    typename ADAPTER::shared_ptr  retVal ( new ADAPTER );

    d_vMeshes.push_back ( retVal );

    if ( mesh_name.size() > 0 )
    {
      d_mMeshNameLookup[mesh_name] = d_vMeshes.size()-1;
    }

    retVal->generateCube ( numNodesPerSide );
    return retVal;
  }

}
}


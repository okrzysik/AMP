
//Libmesh includes
#include "mesh_generation.h"
#include "petsc_matrix.h"
#include "petsc_vector.h"
#include "linear_implicit_system.h"
#include "exodusII_io.h"
#include "boundary_info.h"

//AMP includes
#include "matrices/ManagedPetscMatrix.h"

#include "utils/AMP_MPI.h"
#include "utils/Utilities.h"
#include "utils/AMPManager.h"

#include "vectors/NativePetscVector.h"
#include "vectors/ManagedPetscVector.h"
#include "vectors/trilinos/EpetraVectorEngine.h"
#include "vectors/DualVector.h"
#include "vectors/DualVariable.h"
#include "vectors/MultiVector.h"
#include "vectors/MultiVariable.h"
#include "vectors/NullVariable.h"
#include "vectors/NullVector.h"

#include "LibMeshAdapter.h"
#include "DefaultMeshObjectSorter.h"


namespace AMP { 
namespace Mesh {

  // I love LibMesh.  Everytime I turn around, I see a wonderful
  // design decision.  This ridiculousness is due to LibMeshInit
  // being a variable that, when destroyed, destroys a LibMesh
  // context.  I can see the need for this, given the ridiculous
  // decisions made elsewhere in the package.
  //
  // I second this opinion

void  LibMeshAdapter::initAdapter ( int argc , char **argv , AMP_MPI comm )
{
    AMPManager::initializeLibmesh(comm);
}


void  LibMeshAdapter::finalizeAdapter ()
{
}


LibMeshAdapter::LibMeshAdapter ( const boost::shared_ptr< ::Mesh > &in )
           : d_libMesh ( in )
           , d_libMeshData ( new ::MeshData ( *d_libMesh ) )
           , d_BoundarySet ()
           , d_MeshDatabase ()
{
    d_BoundarySet.build ( *this );
    d_NodeElemMap.build ( *this );
    d_DefaultSorter = AMP::LinearAlgebra::ObjectSorter::shared_ptr ( new DefaultMeshObjectSorter<NodeIterator> ( beginUsedNode() , endUsedNode() ) );
    buildDOFMap ( AMP::LinearAlgebra::Variable::shared_ptr ( new NodalScalarVariable ( "build_comm_list" ) ) );
    buildDOFMap ( AMP::LinearAlgebra::Variable::shared_ptr ( new SingleGaussPointVariable ( "build_comm_list" ) ) );
}

  LibMeshAdapter::LibMeshAdapter ( boost::shared_ptr<Database>  db )
           : d_libMesh ( new ::Mesh ( 3 ) )
           , d_libMeshData ( new ::MeshData ( *d_libMesh ) )
           , d_BoundarySet ()
           , d_MeshDatabase ( db )
  {
  }

  void LibMeshAdapter::readExodusIIFile ( std::string fname )
  {
    d_libMesh->read ( fname );
    d_BoundarySet.build ( *this );
    d_NodeElemMap.build ( *this );
    d_DefaultSorter = AMP::LinearAlgebra::ObjectSorter::shared_ptr ( new DefaultMeshObjectSorter<NodeIterator> ( beginUsedNode() , endUsedNode() ) );
    buildDOFMap ( AMP::LinearAlgebra::Variable::shared_ptr ( new NodalScalarVariable ( "build_comm_list" ) ) );
    buildDOFMap ( AMP::LinearAlgebra::Variable::shared_ptr ( new SingleGaussPointVariable ( "build_comm_list" ) ) );
  }

  void LibMeshAdapter::writeExodusIIFile ( std::string fname )
  {
    d_libMesh->write ( fname );
  }

  void LibMeshAdapter::readIDEASFile ( std::string fname )
  {
    d_libMeshData->activate ();
    d_libMesh->read( fname , d_libMeshData.get() );
    d_BoundarySet.build ( *this );
    d_NodeElemMap.build ( *this );
    d_DefaultSorter = AMP::LinearAlgebra::ObjectSorter::shared_ptr ( new DefaultMeshObjectSorter<NodeIterator> ( beginUsedNode() , endUsedNode() ) );
    buildDOFMap ( AMP::LinearAlgebra::Variable::shared_ptr ( new NodalScalarVariable ( "build_comm_list" ) ) );
    buildDOFMap ( AMP::LinearAlgebra::Variable::shared_ptr ( new SingleGaussPointVariable ( "build_comm_list" ) ) );
  }

  void LibMeshAdapter::writeIDEASFile ( std::string fname )
  {
    d_libMesh->write ( fname );
  }

  void LibMeshAdapter::generateCube ( size_t size )
  {
    MeshTools::Generation::build_cube ( *d_libMesh , size , size , size , -1. , 1 , -1 , 1 , -1 , 1 , HEX8 );
    d_libMesh->prepare_for_use ( false );
    d_NodeElemMap.build ( *this );
    d_DefaultSorter = AMP::LinearAlgebra::ObjectSorter::shared_ptr ( new DefaultMeshObjectSorter<NodeIterator> ( beginUsedNode() , endUsedNode() ) );
    buildDOFMap ( AMP::LinearAlgebra::Variable::shared_ptr ( new NodalScalarVariable ( "build_comm_list" ) ) );
    buildDOFMap ( AMP::LinearAlgebra::Variable::shared_ptr ( new SingleGaussPointVariable ( "build_comm_list" ) ) );
  }

  void LibMeshAdapter::buildDOFMap ( AMP::LinearAlgebra::Variable::shared_ptr var )
  {
    // Make DOFMap:
    AMP::LinearAlgebra::CommunicationList::shared_ptr  commList;
    if ( var->isA<NodalVariable> () )
    {
      NodalRowMap::Parameters  *newDofMapParams = new NodalRowMap::Parameters ( beginElement() , endElement ()  );
      newDofMapParams->d_mesh = this;
      newDofMapParams->d_comm = getComm();
      newDofMapParams->d_localsize = d_libMesh->n_local_nodes() * var->DOFsPerObject();
      NodalRowMap::Parameters::shared_ptr nrmParams ( newDofMapParams );
      commList = AMP::LinearAlgebra::CommunicationList::shared_ptr ( new NodalRowMap ( nrmParams ) );
    }
    else
    {
      size_t numLocal = d_libMesh->n_local_elem() * var->DOFsPerObject();
      commList = AMP::LinearAlgebra::CommunicationList::createEmpty ( numLocal , getComm() );
    }

    DOFMap::Parameters  *dofMapParams = new DOFMap::Parameters;
    dofMapParams->d_SortedObjects = d_DefaultSorter;
    dofMapParams->d_Variable = var;
    dofMapParams->d_Comm = getComm();
    dofMapParams->d_CommList = commList;

    DOFMap::Parameters::shared_ptr dofParams ( dofMapParams );
    DOFMap::shared_ptr pDofMap ( new DOFMap ( dofParams ) );
    d_vDOFMapCache[var->variableID()] = pDofMap;
  }

  void LibMeshAdapter::buildVector ( AMP::LinearAlgebra::Variable::shared_ptr var )
  {
    if ( d_vDOFMapCache.find ( var->variableID() ) == d_vDOFMapCache.end() )
    {
      buildDOFMap ( var );
    }

    DOFMap::shared_ptr   pDofMap = d_vDOFMapCache[var->variableID()];

    AMP::LinearAlgebra::ManagedPetscVectorParameters  *mvparams = new AMP::LinearAlgebra::ManagedPetscVectorParameters ();
    AMP::LinearAlgebra::EpetraVectorEngineParameters  *eveparams = new AMP::LinearAlgebra::EpetraVectorEngineParameters ( pDofMap->numLocalElements() ,
                                                                                  pDofMap->numGlobalElements () ,
                                                                                  getComm() );

    int i = 0;
    for ( size_t local_start = pDofMap->beginDOF() ; local_start != pDofMap->endElement() ; local_start++ , i++ )
    {
      eveparams->addMapping ( i , local_start );
    }
    AMP::LinearAlgebra::VectorEngine::BufferPtr t_buffer ( new AMP::LinearAlgebra::VectorEngine::Buffer ( pDofMap->numLocalElements() ) );
    mvparams->d_Engine = AMP::LinearAlgebra::VectorEngine::shared_ptr ( new AMP::LinearAlgebra::EpetraVectorEngine ( AMP::LinearAlgebra::VectorEngineParameters::shared_ptr ( eveparams ) , t_buffer ) );
    mvparams->d_CommList = d_vDOFMapCache[var->variableID()]->getCommunicationList();

    d_vVectorCache[var->variableID()] = AMP::LinearAlgebra::Vector::shared_ptr ( new AMP::LinearAlgebra::ManagedPetscVector ( AMP::LinearAlgebra::VectorParameters::shared_ptr ( mvparams ) ) );
  }

  void LibMeshAdapter::buildMatrix ( AMP::LinearAlgebra::Variable::shared_ptr  var_operand , AMP::LinearAlgebra::Variable::shared_ptr var_result )
  {
    if ( d_vVectorCache.find ( var_operand->variableID() ) == d_vVectorCache.end() )
    {
      buildVector ( var_operand );
    }

    if ( d_vVectorCache.find ( var_result->variableID() ) == d_vVectorCache.end() )
    {
      buildVector ( var_result );
    }



    DOFMap::shared_ptr  pOperandDofMap = d_vDOFMapCache[var_operand->variableID()];
    DOFMap::shared_ptr  pResultDofMap = d_vDOFMapCache[var_result->variableID()];

    AMP::LinearAlgebra::Vector::shared_ptr  pOperandVec = d_vVectorCache[var_operand->variableID()];
    AMP::LinearAlgebra::Vector::shared_ptr  pResultVec = d_vVectorCache[var_result->variableID()];

    AMP::LinearAlgebra::ManagedPetscMatrixParameters *params =
                new AMP::LinearAlgebra::ManagedPetscMatrixParameters ( pResultDofMap->numLocalElements() ,
                                                   pResultDofMap->numGlobalElements() ,
                                                   0 ,
                                                   pOperandDofMap->numGlobalElements() ,
                                                   0 ,
                                                   getComm() );

    AMP_INSIST ( var_operand->isA<NodalVariable>() , "Non-square matrices are only implemented for nodal variables" );
    AMP_INSIST ( var_result->isA<NodalVariable>() , "Non-square matrices are only implemented for nodal variables" );

    int multiplier = var_operand->DOFsPerObject();
    int divisor = var_result->DOFsPerObject();
    NodalRowMap rowMap = pResultDofMap->getCommunicationList()->castTo<NodalRowMap>();
    size_t numLocalElements = pResultDofMap->numLocalElements();
    for ( unsigned int i = 0 ; i != numLocalElements; i++ )
    {
      params->addMapping ( i , pResultDofMap->beginDOF() + i );
      size_t nnz = rowMap.getNNZ ( i );
      params->setEntriesInRow ( i , nnz * multiplier/divisor );
      params->addColumns ( nnz * multiplier/divisor , (int *)rowMap.getColumns ( i*multiplier/divisor ) );
    }

    params->d_CommListLeft = d_vDOFMapCache[var_result->variableID()]->getCommunicationList();
    params->d_CommListRight = d_vDOFMapCache[var_operand->variableID()]->getCommunicationList();
    AMP::LinearAlgebra::Matrix::shared_ptr  newMatrix = AMP::LinearAlgebra::Matrix::shared_ptr ( new AMP::LinearAlgebra::ManagedPetscMatrix ( AMP::LinearAlgebra::Matrix::ParametersPtr ( params ) ) );
    size_t  mat_id = (var_result->variableID() << 10) + var_operand->variableID();
    d_vMatrixCache[mat_id] = newMatrix;

    double  values[1000];  // A little bit of a hack...
    for ( size_t i = 0 ; i != 1000 ; i++ )
      values[i] = 0.0;
    NodalRowMap rowMap2 = pOperandDofMap->getCommunicationList()->castTo<NodalRowMap>();
    for ( size_t i=0 ; i!=numLocalElements; i++ )
    {
      int cur_row_id = (int)i + pResultDofMap->beginDOF();
      int related_col = i * multiplier/divisor;
      newMatrix->castTo<AMP::LinearAlgebra::ManagedMatrix>().createValuesByGlobalID ( 1 ,
                                          rowMap2.getNNZ (related_col) ,
                                         &cur_row_id ,
                                   (int *)rowMap2.getColumns(related_col) ,
                                          values );
    }
    newMatrix->castTo<AMP::LinearAlgebra::EpetraMatrix>().setEpetraMaps ( pResultVec , pOperandVec );
    newMatrix->makeConsistent ();
  }

  void LibMeshAdapter::makeDataConsistent ()
  {
    std::map<std::string , AMP::LinearAlgebra::Vector::shared_ptr>::iterator curVec = d_mData.begin();
    while ( curVec != d_mData.end() )
    {
      curVec->second->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
      curVec++;
    }
  }

  void LibMeshAdapter::buildMatrix ( AMP::LinearAlgebra::Variable::shared_ptr  var )
  {
    if ( d_vVectorCache.find ( var->variableID() ) == d_vVectorCache.end() )
    {
      buildVector ( var );
    }


    DOFMap::shared_ptr  pDofMap = d_vDOFMapCache[var->variableID()];

    AMP::LinearAlgebra::ManagedPetscMatrixParameters *params =
                new AMP::LinearAlgebra::ManagedPetscMatrixParameters ( pDofMap->numLocalElements() ,
                                                   pDofMap->numGlobalElements() ,
                                                   0 ,
                                                   getComm() );

    NodalRowMap rowMap =  pDofMap->getCommunicationList()->castTo<NodalRowMap>();
    size_t beginDOF = pDofMap->beginDOF();
    size_t numLocalElements = pDofMap->numLocalElements();
    for (unsigned int i=0; i!=numLocalElements; i++)
    {
      params->addMapping ( i, beginDOF + i );
      params->setEntriesInRow ( i, rowMap.getNNZ( i ) );
    }
    params->d_CommListLeft = pDofMap->getCommunicationList();
    params->d_CommListRight = pDofMap->getCommunicationList();

    AMP::LinearAlgebra::Matrix::shared_ptr  newMatrix = AMP::LinearAlgebra::Matrix::shared_ptr ( new AMP::LinearAlgebra::ManagedPetscMatrix ( AMP::LinearAlgebra::Matrix::ParametersPtr ( params ) ) );
    d_vMatrixCache[var->variableID()] = newMatrix;

    double  values[1000];  // A little bit of a hack...
    for ( size_t i = 0 ; i != 1000 ; i++ )
      values[i] = 0.0;
    for ( size_t i = 0 ; i != pDofMap->numLocalElements() ; i++ )
    {
      int cur_row_id = (int) (i + beginDOF);
      int  nnz = rowMap.getNNZ(i);
      int *cols = reinterpret_cast<int *> (rowMap.getColumns(i));
      newMatrix->castTo<AMP::LinearAlgebra::ManagedMatrix>().createValuesByGlobalID ( 1 ,nnz , &cur_row_id , cols , values );
    }
    newMatrix->makeConsistent ();
  }


  AMP::LinearAlgebra::CommunicationList::shared_ptr  LibMeshAdapter::getNodalCommunicationList ()
  {
    AMP::LinearAlgebra::Variable::shared_ptr t ( new NodalScalarVariable ( "temp" ) );
    AMP_ASSERT ( d_vDOFMapCache.find ( t->variableID() ) != d_vDOFMapCache.end() );
    return d_vDOFMapCache [t->variableID()]->getCommunicationList();
  }

  AMP::LinearAlgebra::Vector::shared_ptr   LibMeshAdapter::createVector ( AMP::LinearAlgebra::Variable::shared_ptr variable )
  {
    AMP_ASSERT ( variable );
    AMP::LinearAlgebra::Vector::shared_ptr  retVal;
    if ( variable->isA<AMP::LinearAlgebra::NullVariable>() )
    {
      retVal = AMP::LinearAlgebra::NullVector::create ( variable );
    }
    else if ( variable->isA<NodalVariable>() || variable->isA<IntegrationPointVariable>() )
    {
      if ( d_vVectorCache.find ( variable->variableID () ) == d_vVectorCache.end() )
      {
        buildVector ( variable );
      }

      retVal = d_vVectorCache[variable->variableID()]->cloneVector ( variable );
    }
    else if ( variable->isA<AMP::LinearAlgebra::DualVariable> () )
    {
      retVal = AMP::LinearAlgebra::DualVector::create ( createVector ( variable->castTo<AMP::LinearAlgebra::DualVariable>().first() ) ,
                                    createVector ( variable->castTo<AMP::LinearAlgebra::DualVariable>().second() ) ,
                                    variable );
    }
    else if ( variable->isA<AMP::LinearAlgebra::MultiVariable> () )
    {
      AMP::LinearAlgebra::MultiVariable &in_var = variable->castTo<AMP::LinearAlgebra::MultiVariable> ();
      retVal = AMP::LinearAlgebra::MultiVector::create ( variable , getComm() );
      for ( size_t i = 0 ; i != in_var.numVariables () ; i++ )
      {
        retVal->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( createVector ( in_var.getVariable ( i ) ) );
      }
    }
    else
    {
      AMP_ERROR( "Oh well, I don't know what to do here..." );
    }

    return retVal;
  }

  AMP::LinearAlgebra::Matrix::shared_ptr   LibMeshAdapter::createMatrix ( AMP::LinearAlgebra::Variable::shared_ptr operand , AMP::LinearAlgebra::Variable::shared_ptr result )
  {
    AMP::LinearAlgebra::Matrix::shared_ptr  retVal;
    bool bSquare = true;
    if ( result )
      if ( result->variableID() != operand->variableID() )
        bSquare = false;
    if ( bSquare )
    {
      if ( operand->isA<NodalVariable>() )
      {
        if ( d_vMatrixCache.find ( operand->variableID () ) == d_vMatrixCache.end() )
        {
          buildMatrix ( operand );
        }

        retVal = d_vMatrixCache[operand->variableID()]->cloneMatrix ();
      }
      else
      {
        AMP_ERROR( "Matrices on integration points?!?!?!?" );
      }
    }
    else
    {
      size_t  mat_id = (result->variableID() << 10) + operand->variableID();
      if ( d_vMatrixCache.find ( mat_id ) == d_vMatrixCache.end() )
      {
        buildMatrix ( operand , result );
      }

      retVal = d_vMatrixCache[mat_id]->cloneMatrix();
    }
    return retVal;
  }

  DOFMap::shared_ptr   LibMeshAdapter::getDOFMap ( AMP::LinearAlgebra::Variable::shared_ptr variable )
  {
    DOFMap::shared_ptr  retVal;
    if ( variable->isA<NodalVariable>() || variable->isA<IntegrationPointVariable>() )
    {
      if ( d_vDOFMapCache.find ( variable->variableID () ) == d_vDOFMapCache.end() )
      {
        buildDOFMap ( variable );
      }

      retVal = d_vDOFMapCache[variable->variableID()];
    }
    else
    {
      AMP_ERROR( "No DOFMap on integration points just yet..." );
    }
    return retVal;
  }


  std::vector<short int>    LibMeshAdapter::getBoundaryIds ( Node n )
  {
    ::Node  &t = n.getNode();
    return  d_libMesh->boundary_info->boundary_ids ( &t );
  }

  bool        LibMeshAdapter::isOnBoundary ( Node n , short int b_id )
  {
    std::vector<short int>  b_ids = getBoundaryIds ( n );
    for ( size_t i = 0 ; i != b_ids.size() ; i++ )
    {
      if ( b_ids[i] == b_id )
      {
        return true;
      }
    }
    return false;
  }

  void    LibMeshAdapter::translate ( double x , double y , double z )
  {
    NodeIterator  curNode = beginUsedNode();
    while ( curNode != endUsedNode () )
    {
      curNode->translate ( x , y , z );
      curNode++;
    }
  }

  void    LibMeshAdapter::scale ( double alpha )
  {
    NodeIterator  curNode = beginUsedNode();
    while ( curNode != endUsedNode () )
    {
      curNode->x() *= alpha;
      curNode->y() *= alpha;
      curNode->z() *= alpha;
      curNode++;
    }
  }

  const std::set<short int> &   LibMeshAdapter::getBoundaryIds ()
  {
    return d_libMesh->boundary_info->get_boundary_ids();
  }

  AMP::LinearAlgebra::Vector::shared_ptr   LibMeshAdapter::getPositionVector ( std::string name )
  {
    AMP::LinearAlgebra::Variable::shared_ptr  var ( new Nodal3VectorVariable ( name ) );
    AMP::LinearAlgebra::Vector::shared_ptr  retVal = createVector ( AMP::LinearAlgebra::Variable::shared_ptr ( new Nodal3VectorVariable ( name ) ) );
    DOFMap::shared_ptr  dofmap = getDOFMap ( var );
    NodeIterator  cur_node = beginNode ();
    NodeIterator  end_node = endNode ();
    std::vector<unsigned int> ndx(100,0);   // Initialize the vector to allocate memory
    std::vector<unsigned int> vn(100,0);    // Initialize the vector to allocate memory
    while ( cur_node != end_node )
    {
      if ( cur_node->isOwned() )
      {
        ndx.clear();  // Clear the vector (does not deallocate the memory)
        vn.clear();   // Clear the vector (does not deallocate the memory)
        dofmap->getDOFs ( *cur_node , ndx , vn );
        if ( ndx.size() > 0 ) // ??Why would getDOFs fail??
        {
          retVal->setValueByGlobalID ( ndx[0] , cur_node->x() );
          retVal->setValueByGlobalID ( ndx[1] , cur_node->y() );
          retVal->setValueByGlobalID ( ndx[2] , cur_node->z() );
        }
      }
      cur_node++;
    }
    retVal->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    return retVal;
  }

  void LibMeshAdapter::displaceMesh ( const AMP::LinearAlgebra::Vector::shared_ptr disp )
  {
    AMP_ASSERT ( disp->getVariable()->isA<Nodal3VectorVariable> () );
    // Make sure the displacement is consistent
    disp->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
    // Update the locally owned nodes
    AMP::LinearAlgebra::Vector::const_iterator  curData = disp->begin();
    unsigned int  curNode = disp->getCommunicationList()->getStartGID() / 3;
    while ( curData != disp->end() )
    {
      getNode ( curNode ).x() += (*curData);
      curData++;
      getNode ( curNode ).y() += (*curData);
      curData++;
      getNode ( curNode ).z() += (*curData);
      curData++;
      curNode++;
    }
    // Update the ghost nodes
    const std::vector<unsigned int> &ghosts = disp->getCommunicationList()->getGhostIDList ();
    std::vector<double> ghost_disp(ghosts.size(),0);
    if ( ghosts.size()>0 )
        disp->getValuesByGlobalID(ghosts.size(),(int*)&ghosts[0],&ghost_disp[0]);
    AMP_ASSERT(ghosts.size()%3==0);
    for (size_t i=0; i<ghosts.size()/3; i++) {
        curNode = ghosts[3*i]/3;
        getNode ( curNode ).x() += ghost_disp[3*i+0];
        getNode ( curNode ).y() += ghost_disp[3*i+1];
        getNode ( curNode ).z() += ghost_disp[3*i+2];
    }
    /*std::vector<unsigned int>::const_iterator curGhost = ghosts.begin();
    while ( curGhost != ghosts.end() )
    {
      curNode = *curGhost / 3;
      getNode ( curNode ).x() += (*curGhost);
      curGhost++;
      getNode ( curNode ).y() += (*curGhost);
      curGhost++;
      getNode ( curNode ).z() += (*curGhost);
      curGhost++;
    }*/
  }

}
}


#include "utils/AMP_MPI.h"
// \cond EXLUDE_CODE_FROM_DOXYGEN

namespace AMP { 
namespace Mesh {


template <typename MANAGER>
void SiloIO<MANAGER>::computeNodeRenumber ( MeshAdapter &m )
{
    size_t i = 0;
    d_NodeRenumber.clear();
    AMP::LinearAlgebra::CommunicationList::shared_ptr commList = m.getNodalCommunicationList ();
    const std::vector<unsigned int> &ghosts = commList->getGhostIDList ();
    std::vector<unsigned int>::const_iterator curGhost = ghosts.begin();
    while ( curGhost != ghosts.end() )
    {
      d_NodeRenumber[*curGhost] = i;
      i++;
      curGhost++;
    }

    for ( size_t j = 0 ; j != commList->numLocalRows() ; j++ )
    {
      d_NodeRenumber [ j + commList->getStartGID() ] = i;
      i++;
    }
}


template <typename MANAGER>
void  SiloIO<MANAGER>::readNodalVector ( MeshAdapter &m , AMP::LinearAlgebra::Vector::shared_ptr  p , const std::string & , const std::string &v_name )
{
    DBucdvar *data = DBGetUcdvar ( d_FileHandle , v_name.c_str() );
    AMP_INSIST ( data , "Cannot find data" );

    AMP_INSIST ( data->nels == (int)d_NodeRenumber.size() , "Data is wrong size" );
    AMP_INSIST ( data->nvals == (int)p->getVariable()->DOFsPerObject() , "Data is wrong dimension" );

    DOFMap::shared_ptr  dof_map = m.getDOFMap ( p->getVariable() );
    NodeIterator  curNode = m.beginUsedNode();
    std::vector<unsigned int> empty_vec;
    double **d_ptr = reinterpret_cast<double **> ( data->vals );
    while ( curNode != m.endUsedNode() )
    {
      if ( d_NodeRenumber.find(curNode->globalID()) != d_NodeRenumber.end() )
      {
        std::vector<unsigned int> dofs;
        dof_map->getDOFs ( *curNode , dofs , empty_vec );
        for ( unsigned int i = 0 ; i != dofs.size() ; i++ )
        {
          p->setValueByGlobalID ( dofs[i] , d_ptr[i][d_NodeRenumber[curNode->globalID()]] );
        }
      }
      curNode++;
    }
    DBFreeUcdvar ( data );
}


template <typename MANAGER>
void  SiloIO<MANAGER>::writeNodalVector ( MeshAdapter &m , AMP::LinearAlgebra::Vector::shared_ptr  p , const std::string &mesh_name , const std::string &v_name )
{
    DBoptlist  *opts = DBMakeOptlist ( 3 );
    DBAddOption ( opts , DBOPT_DTIME , &d_Time );
    DBAddOption ( opts , DBOPT_CYCLE , &d_IterCount );
    DBAddOption ( opts , DBOPT_UNITS , const_cast<char *> (p->getVariable()->getUnits().c_str()) );
    int dimLen = p->getVariable()->DOFsPerObject();
    int vecLen = d_NodeRenumber.size();
    double **values = new double * [p->getVariable()->DOFsPerObject()];
    for ( unsigned int i = 0 ; i != p->getVariable()->DOFsPerObject() ; i++ )
      values[i] = new double [vecLen];

    DOFMap::shared_ptr  dof_map = m.getDOFMap ( p->getVariable() );
    NodeIterator  curNode = m.beginUsedNode();
    std::vector<unsigned int> empty_vec;
    while ( curNode != m.endUsedNode() )
    {
      if ( d_NodeRenumber.find(curNode->globalID()) != d_NodeRenumber.end() )
      {
        std::vector<unsigned int> dofs;
        dof_map->getDOFs ( *curNode , dofs , empty_vec );
        for ( unsigned int i = 0 ; i != dofs.size() ; i++ )
        {
          values[i][d_NodeRenumber[curNode->globalID()]] = p->getValueByGlobalID ( dofs[i] );
        }
      }
      curNode++;
    }

    char **coordnames = new char *[3];
    coordnames[0] = strdup ( "X" );
    coordnames[1] = strdup ( "Y" );
    coordnames[2] = strdup ( "Z" );

    DBPutUcdvar ( d_FileHandle ,
                  v_name.c_str() ,
                  mesh_name.c_str() ,
                  dimLen ,
                  coordnames ,
                  values ,
                  vecLen ,
                  0 ,
                  0 ,
                  DB_DOUBLE ,
                  DB_NODECENT ,
                  opts );
    DBFreeOptlist ( opts );
    for ( unsigned int i = 0 ; i != p->getVariable()->DOFsPerObject() ; i++ )
      delete [] values[i];
    delete [] values;
    free ( coordnames[0] );
    free ( coordnames[1] );
    free ( coordnames[2] );
    delete [] coordnames;
}


template <typename MANAGER>
void  SiloIO<MANAGER>::readCellVector ( MeshAdapter &m , AMP::LinearAlgebra::Vector::shared_ptr  p , const std::string & , const std::string &v_name )
{
    DBucdvar *data = DBGetUcdvar ( d_FileHandle , v_name.c_str() );
    AMP_INSIST ( data , "Cannot find data" );

    AMP_INSIST ( data->nels == (int)p->getLocalSize() , "Data is wrong size" );
    int dimLen = p->getVariable()->DOFsPerObject();
    AMP_INSIST ( data->nvals == dimLen , "Data is wrong dimension" );

    DOFMap::shared_ptr  dof_map;
    if ( p->getVariable()->isA<IntegrationPointVariable>() )
      dof_map = m.getDOFMap ( p->getVariable() );
    else
      dof_map = m.getDOFMap ( AMP::LinearAlgebra::Variable::shared_ptr ( new RunTimeIntegrationPointVariable ( "..." , p->getVariable()->DOFsPerObject() ) ) );
    ElementIterator  curElem = m.beginElement();
    std::vector<unsigned int> empty_vec;
    size_t curElemID = 0;
    double ** d_ptr = reinterpret_cast<double **> ( data->vals );
    while ( curElem != m.endElement() )
    {
      for ( int i = 0 ; i != dimLen ; i++ )
      {
        int ID = dof_map->getGlobalID ( curElem->globalID() , i );
        p->setValueByGlobalID ( ID , d_ptr[i][curElemID] );
      }
      curElemID++;
      curElem++;
    }
    DBFreeUcdvar ( data );
}


template <typename MANAGER>
void  SiloIO<MANAGER>::writeCellVector ( MeshAdapter &m , AMP::LinearAlgebra::Vector::shared_ptr  p , const std::string &mesh_name , const std::string &v_name )
{
    DBoptlist  *opts = DBMakeOptlist ( 3 );
    DBAddOption ( opts , DBOPT_DTIME , &d_Time );
    DBAddOption ( opts , DBOPT_CYCLE , &d_IterCount );
    DBAddOption ( opts , DBOPT_UNITS , const_cast<char *> (p->getVariable()->getUnits().c_str()) );
    int dimLen = p->getVariable()->DOFsPerObject();
    int vecLen = p->getLocalSize();
    double **values = new double * [p->getVariable()->DOFsPerObject()];
    for ( unsigned int i = 0 ; i != p->getVariable()->DOFsPerObject() ; i++ )
      values[i] = new double [vecLen];

    DOFMap::shared_ptr dof_map;
    if ( p->getVariable()->isA<IntegrationPointVariable>() )
      dof_map = m.getDOFMap ( p->getVariable() );
    else
      dof_map = m.getDOFMap ( AMP::LinearAlgebra::Variable::shared_ptr ( new RunTimeIntegrationPointVariable ( "..." , p->getVariable()->DOFsPerObject() ) ) );
    ElementIterator  curElem = m.beginElement();
    std::vector<unsigned int> empty_vec;
    size_t curElemID = 0;
    while ( curElem != m.endElement() )
    {
      for ( int i = 0 ; i != dimLen ; i++ )
      {
        values[i][curElemID] = p->getValueByGlobalID (dof_map->getGlobalID ( curElem->globalID() , i ));
      }
      curElemID++;
      curElem++;
    }

    char **coordnames = new char *[3];
    coordnames[0] = strdup ( "X" );
    coordnames[1] = strdup ( "Y" );
    coordnames[2] = strdup ( "Z" );

    DBPutUcdvar ( d_FileHandle ,
                  v_name.c_str() ,
                  mesh_name.c_str() ,
                  dimLen ,
                  coordnames ,
                  values ,
                  vecLen ,
                  0 ,
                  0 ,
                  DB_DOUBLE ,
                  DB_ZONECENT ,
                  opts );
    DBFreeOptlist ( opts );

    free ( coordnames[0] );
    free ( coordnames[1] );
    free ( coordnames[2] );
    delete [] coordnames;

    for ( unsigned int i = 0 ; i != p->getVariable()->DOFsPerObject() ; i++ )
      delete [] values[i];
    delete [] values;
}


template <typename MANAGER>
void  SiloIO<MANAGER>::readVector ( MeshAdapter &m , AMP::LinearAlgebra::Vector::shared_ptr  p , const std::string &mesh_name , const std::string &v_name )
{
    if ( p->getVariable()->isA<NodalVariable>() )
    {
      readNodalVector ( m , p , mesh_name , v_name );
    }
    else if ( p->getVariable()->isA<IntegrationPointVariable> () )
    {
      readCellVector ( m , p , mesh_name , v_name );
    }
    else if ( p->getVariable()->isA<AMP::LinearAlgebra::StridedVariable> () )
    {
      readCellVector ( m , p , mesh_name , v_name );
    }
    else
    {
      AMP_ERROR( "Cannot write vector of this type" );
    }
}


template <typename MANAGER>
void  SiloIO<MANAGER>::writeVector ( MeshAdapter &m , AMP::LinearAlgebra::Vector::shared_ptr  p , const std::string &mesh_name , const std::string &v_name )
{
    if ( p->getVariable()->isA<NodalVariable>() )
    {
      writeNodalVector ( m , p , mesh_name , v_name );
    }
    else if ( p->getVariable()->isA<IntegrationPointVariable> () )
    {
      writeCellVector ( m , p , mesh_name , v_name );
    }
    else if ( p->getVariable()->isA<AMP::LinearAlgebra::StridedVariable> () )
    {
      writeCellVector ( m , p , mesh_name , v_name );
    }
    else
    {
      AMP_ERROR( "Cannot write vector of this type" );
    }
}


template <typename MANAGER>
int  SiloIO<MANAGER>::writeConnectivityList ( MeshAdapter &m , const std::string &zone_name )
{
    int   nodes_per_elem = m.beginElement()->numNodes();
    int   nodelist_max = m.numLocalElements()*nodes_per_elem;
    int  *nodelist = new int [ nodelist_max ] ;
    int   cur_spot = 0;
    int   num_elems = 0;
    ElementIterator  curElem = m.beginElement();
    while (( curElem != m.endElement()) || (cur_spot != nodelist_max ))
    {
      bool keep_element = true;
      for ( int i = 0 ; i != nodes_per_elem ; i++ )
      {
        if ( d_NodeRenumber.find ( curElem->getNodeID(i) ) == d_NodeRenumber.end() )
        {
          keep_element = false;
        }
        else
        {
          nodelist[cur_spot] = d_NodeRenumber[curElem->getNodeID(i)];
        }
        cur_spot++;
      }
      if ( !keep_element )
      {
        cur_spot -= curElem->numNodes();
      }
      else
      {
        num_elems++;
      }
      curElem++;
    }
    int t0 = -1;
    switch ( nodes_per_elem )
    {
      case 8:  t0 = DB_ZONETYPE_HEX;
               break;
               // Dimensionality of libMesh 2D meshes isn't working properly
      case 4:  t0 = DB_ZONETYPE_QUAD;
               break;
    }
    int t1 = nodes_per_elem;
    int t2 = num_elems;
    DBPutZonelist2 ( d_FileHandle ,
                     zone_name.c_str() ,
                     num_elems ,
                     m.geometricDimension() ,
                     nodelist ,
                     cur_spot ,
                     0 ,
                     0 ,
                     0 ,
                     &t0 ,
                     &t1 ,
                     &t2 ,
                     1 ,
                     0 );
    delete [] nodelist;
    return num_elems;
}


template <typename MANAGER>
void  SiloIO<MANAGER>::readMesh ( MeshAdapter &m , const std::string &mesh_name )
{
    computeNodeRenumber ( m );
    DataIterator curVec = m.beginData ();

    while ( curVec != m.endData() )
    {
      std::string varname;
      varname = curVec->first;
      varname += "_";
      varname += mesh_name;
      readVector ( m , curVec->second , mesh_name , varname );
      curVec++;
    }
}


template <typename MANAGER>
void  SiloIO<MANAGER>::writeMesh ( MeshAdapter &m , const std::string &mesh_name )
{
    AMP_MPI globalComm(AMP_COMM_WORLD);
    std::stringstream  zone_name;
    zone_name << mesh_name << "_zone";
    if ( d_SyncForMasterFile ) {
        if ( globalComm.getRank() == 0 ) {
            std::stringstream zn;
            zn << d_FileName << ":" << d_SubdirName << "/" << mesh_name;
            char *zname = new char [ zn.str().size() + 1];
            strcpy ( zname , zn.str().c_str() );
            d_MeshName.push_back ( zname );
        } else {
            d_OutstandingRequests.push_back( globalComm.Isend(mesh_name.c_str(),mesh_name.size(),0,43682) );            
        }
    }

    computeNodeRenumber ( m );
    int numWrittenElements = writeConnectivityList ( m , zone_name.str() );
    double **nodes = new double * [3];
    nodes[0] = new double [d_NodeRenumber.size()];
    nodes[1] = new double [d_NodeRenumber.size()];
    nodes[2] = new double [d_NodeRenumber.size()];
    NodeIterator curNode = m.beginUsedNode();
    while ( curNode != m.endUsedNode() )
    {
        std::vector<unsigned int>  node_id;
        if ( d_NodeRenumber.find ( curNode->globalID() ) == d_NodeRenumber.end () )
        {
            curNode++;
            continue;
        }
        nodes[0][d_NodeRenumber[curNode->globalID()]] = curNode->x();
        nodes[1][d_NodeRenumber[curNode->globalID()]] = curNode->y();
        nodes[2][d_NodeRenumber[curNode->globalID()]] = curNode->z();
        curNode++;
    }
    char **coordnames = new char *[3];
    coordnames[0] = strdup ( "X" );
    coordnames[1] = strdup ( "Y" );
    coordnames[2] = strdup ( "Z" );
    DBPutUcdmesh ( d_FileHandle ,
                   mesh_name.c_str() ,
                   3 ,
                   coordnames ,
                   nodes ,
                   d_NodeRenumber.size() ,
                   numWrittenElements ,
                   zone_name.str().c_str() ,
                   0 ,
                   DB_DOUBLE ,
                   0 );
    DataIterator curVec = m.beginData ();
    if ( d_SyncForMasterFile ) {
        int t = m.sizeData();
        if ( globalComm.getRank() != 0 ) {
            d_OutstandingRequests.push_back( globalComm.Isend(&t,1,0,43683) );
        }
    }

    while ( curVec != m.endData() )
    {
        std::string varname;
        varname = curVec->first;
        varname += "_";
        varname += mesh_name;
        if ( d_SyncForMasterFile ) {
            if ( globalComm.getRank() == 0 ) {
                std::stringstream zn;
                zn << d_FileName << ":" << d_SubdirName << "/" << varname;
                char *zname = new char [ zn.str().size() + 1 ];
                strcpy ( zname , zn.str().c_str() );
                d_VecNames[varname].push_back ( zname );
            } else {
                d_OutstandingRequests.push_back( globalComm.Isend(varname.c_str(),varname.size(),0,43684) );
            }
        }
        writeVector ( m , curVec->second , mesh_name , varname );
        curVec++;
    }
    delete[] nodes[0];
    delete[] nodes[1];
    delete[] nodes[2];
    delete[] nodes;
    free ( coordnames[0] );
    free ( coordnames[1] );
    free ( coordnames[2] );
    delete[] coordnames;
}


template <typename MANAGER>
void  SiloIO<MANAGER>::updateFile ( MeshManager &m , std::string &fname , double t , unsigned int ic )
{
    d_Time = t;
    d_IterCount = ic;
    d_FileHandle = DBOpen ( fname.c_str() , DB_HDF5 , DB_APPEND );

    MeshNameIterator  curName = m.beginMeshNames();
    while ( curName != m.endMeshNames() )
    {
      DataIterator curVec = m.getMesh ( curName->first )->beginData ();
      while ( curVec != m.getMesh( curName->first )->endData() )
      {
        writeVector ( *m.getMesh ( curName->first ) , *curVec , curName->first.c_str() , curVec->first );
        curVec++;
      }
      curName++;
    }

    DBClose ( d_FileHandle );
}


template <typename MANAGER>
void  SiloIO<MANAGER>::readFile ( MeshManager &m , const std::string &fname )
{
    d_Time = 0.0;
    d_IterCount = 0;
    d_SyncForMasterFile = false;

    AMP_MPI meshComm(m.getMeshComm());
    AMP_MPI globalComm(AMP_COMM_WORLD);

    d_FileHandle = DBOpen ( fname.c_str() , DB_HDF5 , DB_APPEND );

    std::stringstream dir_name;
    dir_name << "rank_" << globalComm.getRank();
    DBSetDir ( d_FileHandle , dir_name.str().c_str() );

    MeshNameIterator  curName = m.beginMeshNames();
    while ( curName != m.endMeshNames() )
    {
      readMesh ( *m.getMesh ( curName->first ) , curName->first );
      curName++;
    }

    DBClose ( d_FileHandle );
}


template <typename MANAGER>
void  SiloIO<MANAGER>::writeFile ( MeshManager &m , const std::string &fname )
{
    d_Time = 0.0;
    d_IterCount = 0;
    d_SyncForMasterFile = false;
    AMP_MPI meshComm(m.getMeshComm());
    AMP_MPI globalComm(AMP_COMM_WORLD);

    for (int i=0; i<meshComm.getSize(); i++) {
        if ( meshComm.getRank() == i ) {
            if ( meshComm.getRank() == 0 ) {
                d_FileHandle = DBCreate ( fname.c_str() , DB_CLOBBER , DB_LOCAL , NULL , DB_HDF5 );
            } else {
                d_FileHandle = DBOpen ( fname.c_str() , DB_HDF5 , DB_APPEND );
            }

            std::stringstream dir_name;
            dir_name << "rank_" << globalComm.getRank();
            DBMkDir ( d_FileHandle , dir_name.str().c_str() );
            DBSetDir ( d_FileHandle , dir_name.str().c_str() );

            MeshNameIterator  curName = m.beginMeshNames();
            while ( curName != m.endMeshNames() ) {
                writeMesh ( *m.getMesh ( curName->first ) , curName->first );
                curName++;
            }

            DBClose ( d_FileHandle );
        }
        meshComm.barrier();
    }
}


template <typename MANAGER>
void SiloIO<MANAGER>::flushMPIRequests ()
{
    // Wait for all MPI_Requests to finish
	if ( d_OutstandingRequests.size() > 0 ) {
        AMP::AMP_MPI::waitAll(d_OutstandingRequests.size(),&(d_OutstandingRequests[0]));
        d_OutstandingRequests.resize ( 0 );
    }
}


template <typename MANAGER>
void  SiloIO<MANAGER>::decompReadFile ( MeshManager &m , MeshAdapter &adapt , const std::string &fname , const std::string & )
{
    d_Time = 0.0;
    d_IterCount = 0;
    d_SyncForMasterFile = true;
    //int gsize;
    AMP_MPI globalComm(AMP_COMM_WORLD);

    d_FileHandle = DBOpen ( fname.c_str() , DB_HDF5 , DB_READ );
    std::stringstream dir_name;
    dir_name << "rank_" << globalComm.getRank();
    DBSetDir ( d_FileHandle , dir_name.str().c_str() );
    readMesh ( adapt , "mesh" );
    DBClose ( d_FileHandle );
}


template <typename MANAGER>
void  SiloIO<MANAGER>::decompWriteFile ( MeshManager &m , MeshAdapter &adapt , const std::string &fname , const std::string &mastername )
{
    d_Time = 0.0;
    d_IterCount = 0;
    d_SyncForMasterFile = true;
    //int gsize;
    AMP_MPI meshComm(m.getMeshComm());
    AMP_MPI globalComm(AMP_COMM_WORLD);

    for ( int i = 0 ; i != meshComm.getSize() ; i++ ) {
        if ( meshComm.getRank() == i ) {
            if ( meshComm.getRank() == 0 ) {
                d_FileHandle = DBCreate ( fname.c_str() , DB_CLOBBER , DB_LOCAL , NULL , DB_HDF5 );
            } else {
                d_FileHandle = DBOpen ( fname.c_str() , DB_HDF5 , DB_APPEND );
            }

            if ( globalComm.getRank() == 0 ) {
                d_FileName = fname;
            } else {
                d_OutstandingRequests.push_back( globalComm.Isend(fname.c_str(),fname.size(),0,43680) );
            }

            std::stringstream dir_name;
            dir_name << "rank_" << globalComm.getRank();
            DBMkDir ( d_FileHandle , dir_name.str().c_str() );
            DBSetDir ( d_FileHandle , dir_name.str().c_str() );
            if ( globalComm.getRank() == 0 ) {
                d_SubdirName = dir_name.str();
            } else {
                d_OutstandingRequests.push_back( globalComm.Isend(dir_name.str().c_str(),dir_name.str().size(),0,43681) );
            }

            writeMesh ( adapt , "mesh" );
            DBClose ( d_FileHandle );

        }
        meshComm.barrier();
    }

    if ( globalComm.getRank() == 0 ) {
        writeMasterFile ( mastername );
    }

    flushMPIRequests ();
}


template <typename MANAGER>
void  SiloIO<MANAGER>::writeMasterFile ( const std::string &master_fname )
{
    AMP_MPI globalComm(AMP_COMM_WORLD);
    int totalSize = globalComm.getSize();
    
    for ( int i = 1 ; i != totalSize ; i++ ) {

        int char_size;
        char_size = (int) globalComm.probe(i,43680)/sizeof(char);
        char *fname = new char[char_size+1];
        globalComm.recv(fname,char_size,i,false,43680);
        fname[char_size] = 0;

        char_size = (int) globalComm.probe(i,43681)/sizeof(char);
        char *dirname = new char[char_size+1];
        globalComm.recv(dirname,char_size,i,false,43681);
        dirname[char_size] = 0;

        char_size = (int) globalComm.probe(i,43682)/sizeof(char);
        char *zonename = new char[char_size+1];
        globalComm.recv(zonename,char_size,i,false,43682);
        zonename[char_size] = 0;

        std::stringstream  mesh_path;
        mesh_path << fname << ":" << dirname << "/" << zonename;
        char  *mpath_str = new char[mesh_path.str().size()+1];
        d_MeshName.push_back ( mpath_str );
        strcpy ( mpath_str , mesh_path.str().c_str() );

        int numVecs;
        int size_tmp = 1;
        globalComm.recv(&numVecs,size_tmp,i,false,43683);
        for ( int j = 0 ; j != numVecs ; j++ ) {
            char_size = (int) globalComm.probe(i,43684)/sizeof(char);
            char *vecname = new char[char_size+1];
            globalComm.recv(vecname,char_size,i,false,43684);
            vecname[char_size] = 0;

            std::stringstream  vec_path;
            vec_path << fname << ":" << dirname << "/" << vecname;
            char *vpath_str = new char[vec_path.str().size()+1];
            strcpy ( vpath_str , vec_path.str().c_str() );
            d_VecNames[vecname].push_back ( vpath_str );

            delete [] vecname;
        }

        delete [] zonename;
        delete [] fname;
        delete [] dirname;
    }

    DBfile  *master_handle = DBCreate ( master_fname.c_str() , DB_CLOBBER , DB_LOCAL , NULL , DB_HDF5 );
    std::vector<int>  mesh_types ( d_MeshName.size() , DB_UCDMESH );
    DBPutMultimesh ( master_handle , "whole_mesh" , d_MeshName.size() , &(d_MeshName[0]) , &(mesh_types[0]) , 0 );
    std::map<std::string,std::vector<char *> >::iterator  cur_vec = d_VecNames.begin();
    while ( cur_vec != d_VecNames.end() )
    {
        std::vector<int> vec_types ( cur_vec->second.size() , DB_UCDVAR );
        DBPutMultivar ( master_handle , cur_vec->first.c_str() , cur_vec->second.size() , &(cur_vec->second[0]) , &(vec_types[0]) , 0 );
        cur_vec++;
    }

    DBClose ( master_handle );
    for ( size_t i = 0 ; i != d_MeshName.size() ; i++ )
        delete [] d_MeshName[i];

    cur_vec = d_VecNames.begin();
    while ( cur_vec != d_VecNames.end() )
    {
        for ( size_t i = 0 ; i != cur_vec->second.size() ; i++ )
            delete [] cur_vec->second[i];
        cur_vec++;
    }
}


} // Mesh namespace
} // AMP namespace

// \endcond EXLUDE_CODE_FROM_DOXYGEN


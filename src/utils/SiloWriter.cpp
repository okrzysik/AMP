#include "AMP/utils/SiloWriter.h"
#include "AMP/ampmesh/MultiMesh.h"
#include "AMP/utils/Utilities.h"

#include "ProfilerApp.h"

#include <chrono>


namespace AMP::Utilities {


static inline size_t find_slash( const std::string &filename )
{
    size_t i1 = filename.find_last_of( 47 );
    size_t i2 = filename.find_last_of( 92 );
    size_t i  = std::string::npos;
    if ( i1 == std::string::npos )
        i = i2;
    else if ( i2 == std::string::npos )
        i = i1;
    else if ( i1 != std::string::npos && i2 != std::string::npos )
        i = std::max( i1, i2 );
    return i;
}


// Function to replace all instances of a string with another
static inline void strrep( std::string &str, const std::string &s, const std::string &r )
{
    size_t i = 0;
    while ( i < str.length() ) {
        i = str.find( s, i );
        if ( i == std::string::npos ) {
            break;
        }
        str.replace( i, s.length(), r );
        i += r.length();
    }
}


/************************************************************
 * Constructor/Destructor                                    *
 ************************************************************/
SiloIO::SiloIO() : AMP::Utilities::Writer()
{
    d_dim = -1;
#ifdef USE_EXT_SILO
    DBSetAllowEmptyObjects( true );
#endif
}
SiloIO::~SiloIO() = default;


/************************************************************
 * Some basic functions                                      *
 ************************************************************/
Writer::WriterProperties SiloIO::getProperties() const
{
    WriterProperties properties;
    properties.type                   = "Silo";
    properties.extension              = "silo";
    properties.registerMesh           = true;
    properties.registerVector         = false;
    properties.registerVectorWithMesh = true;
    properties.registerMatrix         = false;
    return properties;
}


#ifdef USE_EXT_SILO

// Some internal functions
static void createSiloDirectory( DBfile *FileHandle, const std::string &path );


/************************************************************
 * Function to read a silo file                              *
 ************************************************************/
void SiloIO::readFile( const std::string & ) { AMP_ERROR( "readFile is not implimented yet" ); }


/************************************************************
 * Function to write a silo file                             *
 * Note: it appears that only one prcoessor may write to a   *
 * file at a time, and that once a processor closes the file *
 * it cannot reopen it (or at least doing this on the        *
 * processor that created the file creates problems).        *
 ************************************************************/
void SiloIO::writeFile( const std::string &fname_in, size_t cycle, double time )
{
    PROFILE_START( "writeFile" );
    // Create the directory (if needed)
    createDirectories( fname_in );
    // Create the file name
    std::string fname = fname_in + "_" + std::to_string( cycle ) + "." + getExtension();
    // Check that the dimension is matched across all processors
    PROFILE_START( "sync dim", 1 );
    int dim2 = d_comm.maxReduce( d_dim );
    if ( d_dim == -1 )
        d_dim = dim2;
    else
        AMP_ASSERT( dim2 == d_dim );
    d_comm.barrier();
    PROFILE_STOP( "sync dim", 1 );
// Syncronize all vectors
#ifdef USE_AMP_VECTORS
    PROFILE_START( "makeConsistent", 1 );
    for ( auto &elem : d_vectors ) {
        auto localState = elem->getUpdateStatus();
        if ( localState == AMP::LinearAlgebra::VectorData::UpdateState::ADDING )
            elem->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_ADD );
        else
            elem->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    }
    PROFILE_STOP( "makeConsistent", 1 );
#endif
    // Write the data for each base mesh
    if ( d_decomposition == 1 ) {
        // Write all mesh data to the main file
        for ( int i = 0; i < d_comm.getSize(); ++i ) {
            if ( d_comm.getRank() == i ) {
                // Open the file
                DBfile *FileHandle;
                if ( d_comm.getRank() == 0 ) {
                    FileHandle = DBCreate( fname.c_str(), DB_CLOBBER, DB_LOCAL, nullptr, DB_HDF5 );
                } else {
                    FileHandle = DBOpen( fname.c_str(), DB_HDF5, DB_APPEND );
                }
                // Write the base meshes
                for ( auto &baseMesh : d_baseMeshes ) {
                    auto &data = baseMesh.second;
                    data.file  = fname.c_str();
                    AMP_ASSERT( data.id == baseMesh.first );
                    writeMesh( FileHandle, baseMesh.second, cycle, time );
                }
                // Close the file
                DBClose( FileHandle );
            }
            d_comm.barrier();
        }
    } else if ( d_decomposition == 2 ) {
        // Every rank will write a seperate file
        if ( d_comm.getRank() == 0 )
            Utilities::recursiveMkdir( fname_in + "_silo", ( S_IRUSR | S_IWUSR | S_IXUSR ), false );
        d_comm.barrier();
        auto fname_rank = fname_in + "_silo/" + std::to_string( cycle ) + "." +
                          std::to_string( d_comm.getRank() + 1 ) + "." + getExtension();
        DBfile *FileHandle = DBCreate( fname_rank.c_str(), DB_CLOBBER, DB_LOCAL, nullptr, DB_HDF5 );
        // Write the base meshes
        for ( auto &baseMesh : d_baseMeshes ) {
            auto &data = baseMesh.second;
            data.file  = fname_rank.c_str();
            AMP_ASSERT( data.id == baseMesh.first );
            writeMesh( FileHandle, baseMesh.second, cycle, time );
        }
        // Close the file
        DBClose( FileHandle );
    } else {
        AMP_ERROR( "Unknown file decomposition" );
    }
    // Write the summary results (multimeshes, multivariables, etc.)
    if ( d_decomposition != 1 ) {
        if ( d_comm.getRank() == 0 ) {
            DBfile *FileHandle = DBCreate( fname.c_str(), DB_CLOBBER, DB_LOCAL, nullptr, DB_HDF5 );
            DBClose( FileHandle );
        }
        d_comm.barrier();
    }
    writeSummary( fname, cycle, time );
    PROFILE_STOP( "writeFile" );
}


/************************************************************
 * Function to register a mesh with silo                     *
 ************************************************************/
void SiloIO::registerMesh( AMP::Mesh::Mesh::shared_ptr mesh, int level, const std::string &path )
{
    if ( mesh == nullptr )
        return;
    if ( d_dim == -1 )
        d_dim = mesh->getDim();
    else
        AMP_INSIST( d_dim == mesh->getDim(),
                    "All meshes must have the same number of physical dimensions" );
    AMP_INSIST( level >= 0 && level <= 3, "Invalid value for level" );
    auto multimesh = std::dynamic_pointer_cast<AMP::Mesh::MultiMesh>( mesh );
    if ( multimesh.get() == nullptr ) {
        // We are dealing with a single mesh
        siloBaseMeshData data;
        data.id        = mesh->meshID();
        data.mesh      = mesh;
        data.rank      = mesh->getComm().getRank() + 1;
        data.ownerRank = d_comm.getRank();
        data.meshName  = "rank_" + std::to_string( data.rank );
        data.path      = path + mesh->getName() + "_/";
        if ( d_baseMeshes.find( mesh->meshID() ) == d_baseMeshes.end() )
            d_baseMeshes.insert( std::make_pair( mesh->meshID(), data ) );
        // Create and register a multimesh for the current mesh
        if ( level > 0 ) {
            siloMultiMeshData data2;
            data2.id       = mesh->meshID();
            data2.mesh     = mesh;
            data2.name     = path + mesh->getName();
            data.ownerRank = mesh->getComm().bcast( d_comm.getRank(), 0 );
            d_multiMeshes.insert( std::make_pair( data2.id, data2 ) );
        }
        // Create and register a multimesh for the rank
        if ( level == 3 ) {
            // Create a unique id for each rank
            uint64_t tmp_id = mesh->meshID().getData();
            uint64_t root2  = d_comm.getRank() + 1;
            tmp_id          = ( root2 << 48 ) + tmp_id;
            siloMultiMeshData data2;
            data2.id       = AMP::Mesh::MeshID( tmp_id );
            data2.mesh     = mesh;
            data2.name     = path + mesh->getName() + "_/rank_" + std::to_string( data.rank );
            data.ownerRank = d_comm.getRank();
            d_multiMeshes.insert( std::make_pair( data2.id, data2 ) );
        }
    } else {
        // We are dealing with a multimesh, register the current mesh and sub meshes
        int level2 = level;
        if ( level == 1 )
            level2 = 0;
        auto new_path  = path + mesh->getName() + "_/";
        auto submeshes = multimesh->getMeshes();
        for ( auto &submeshe : submeshes )
            registerMesh( submeshe, level2, new_path );
        if ( level > 0 ) {
            siloMultiMeshData data;
            data.id        = mesh->meshID();
            data.mesh      = mesh;
            data.name      = path + mesh->getName();
            data.ownerRank = mesh->getComm().bcast( d_comm.getRank(), 0 );
            d_multiMeshes.insert( std::make_pair( mesh->meshID(), data ) );
        }
    }
}

/************************************************************
 * Function to get the mesh ids to use for registering       *
 ************************************************************/
std::vector<AMP::Mesh::MeshID> SiloIO::getMeshIDs( AMP::Mesh::Mesh::shared_ptr mesh )
{
    std::vector<AMP::Mesh::MeshID> ids;
    auto multimesh = std::dynamic_pointer_cast<AMP::Mesh::MultiMesh>( mesh );
    if ( multimesh.get() == nullptr ) {
        // We are dealing with a single mesh
        ids = std::vector<AMP::Mesh::MeshID>( 1, mesh->meshID() );
    } else {
        // We are dealining with a multimesh
        auto meshes = multimesh->getMeshes();
        for ( auto &meshe : meshes ) {
            auto ids2 = getMeshIDs( meshe );
            for ( auto &elem : ids2 )
                ids.push_back( elem );
        }
    }
    return ids;
}


/************************************************************
 * Function to register a vector with silo                   *
 ************************************************************/
#ifdef USE_AMP_VECTORS
void SiloIO::registerVector( AMP::LinearAlgebra::Vector::shared_ptr vec,
                             AMP::Mesh::Mesh::shared_ptr mesh,
                             AMP::Mesh::GeomType type,
                             const std::string &name_in )
{
    // Return if the vector or mesh is empty
    if ( !vec || !mesh )
        return;
    // Make sure the mesh has been registered
    registerMesh( mesh );
    // Perform some error checking
    auto DOFs = vec->getDOFManager();
    if ( !DOFs )
        return;
    auto it1 = mesh->getIterator( type, 0 );
    auto it2 = DOFs->getIterator();
    auto it3 = AMP::Mesh::Mesh::getIterator( AMP::Mesh::SetOP::Intersection, it1, it2 );
    if ( it1.size() == 0 || it2.size() == 0 )
        return;
    if ( it1.size() != it3.size() )
        AMP_WARNING( "vector does not cover the entire mesh for the given entity type" );
    std::vector<size_t> dofs;
    DOFs->getDOFs( it1->globalID(), dofs );
    int DOFsPerPoint = dofs.size();
    if ( type == AMP::Mesh::GeomType::Vertex )
        it1 = mesh->getIterator( type, 1 );
    for ( const auto &elem : DOFs->getIterator() ) {
        DOFs->getDOFs( elem.globalID(), dofs );
        AMP_ASSERT( (int) dofs.size() == DOFsPerPoint );
    }
    // Register the vector with the appropriate base meshes
    std::string name = vec->type();
    if ( !name.empty() ) {
        name = name_in;
    }
    auto ids = getMeshIDs( mesh );
    for ( auto id : ids ) {
        const auto &it = d_baseMeshes.find( id );
        if ( it == d_baseMeshes.end() )
            continue;
        it->second.varName.push_back( name );
        it->second.varType.push_back( type );
        it->second.varSize.push_back( DOFsPerPoint );
        it->second.vec.push_back( vec );
    }
    // Register the vector with the appropriate multi-meshes
    auto it = d_multiMeshes.find( mesh->meshID() );
    AMP_ASSERT( it != d_multiMeshes.end() );
    it->second.varName.push_back( name );
    // Add the vector to the list of vectors so we can perform makeConsistent
    d_vectors.push_back( vec );
    // Add the variable name to the list of variables
    d_varNames.insert( name );
}
void SiloIO::registerVector( AMP::LinearAlgebra::Vector::shared_ptr, const std::string & )
{
    AMP_ERROR( "SiloIO currently requires a mesh to register a vector with" );
}
#endif
#ifdef USE_AMP_MATRICES
void SiloIO::registerMatrix( AMP::LinearAlgebra::Matrix::shared_ptr, const std::string & )
{
    AMP_ERROR( "SiloIO does not yet support matrices" );
}
#endif


/************************************************************
 * Function to write a mesh                                  *
 ************************************************************/
void SiloIO::writeMesh( DBfile *FileHandle, const siloBaseMeshData &data, int cycle, double time )
{
    PROFILE_START( "writeMesh", 1 );
    auto mesh = data.mesh;
    // Get the zone (element) lists
    PROFILE_START( "writeMesh - get-elements", 2 );
    auto elem_iterator = mesh->getIterator( mesh->getGeomType(), 0 );
    AMP_ASSERT( elem_iterator.size() > 0 );
    auto type     = elem_iterator->globalID().type();
    auto nodes    = elem_iterator->getElements( AMP::Mesh::GeomType::Vertex );
    int shapesize = nodes.size();
    int shapetype;
    if ( shapesize == 8 && type == AMP::Mesh::GeomType::Volume )
        shapetype = DB_ZONETYPE_HEX;
    else if ( shapesize == 4 && type == AMP::Mesh::GeomType::Volume )
        shapetype = DB_ZONETYPE_TET;
    else if ( shapesize == 4 && type == AMP::Mesh::GeomType::Face )
        shapetype = DB_ZONETYPE_QUAD;
    else if ( shapesize == 3 && type == AMP::Mesh::GeomType::Face )
        shapetype = DB_ZONETYPE_TRIANGLE;
    else if ( shapesize == 2 && type == AMP::Mesh::GeomType::Edge )
        shapetype = DB_ZONETYPE_BEAM;
    else
        AMP_ERROR( "Unknown element type" );
    int shapecnt = elem_iterator.size();
    PROFILE_STOP( "writeMesh - get-elements", 2 );
    // Get the node list (unique integer for each node) and coordinates
    PROFILE_START( "writeMesh - get-nodelist", 2 );
    PROFILE_START( "writeMesh - get-nodelist-1", 3 );
    auto node_iterator = mesh->getIterator( AMP::Mesh::GeomType::Vertex, 1 );
    std::vector<AMP::Mesh::MeshElementID> nodelist_ids( node_iterator.size() );
    for ( size_t i = 0; i < node_iterator.size(); ++i, ++node_iterator )
        nodelist_ids[i] = node_iterator->globalID();
    AMP::Utilities::quicksort( nodelist_ids );
    PROFILE_STOP( "writeMesh - get-nodelist-1", 3 );
    PROFILE_START( "writeMesh - get-nodelist-2", 3 );
    double *coord[3];
    for ( int i = 0; i < d_dim; ++i )
        coord[i] = new double[node_iterator.size()];
    node_iterator = mesh->getIterator( AMP::Mesh::GeomType::Vertex, 1 );
    for ( size_t i = 0; i < node_iterator.size(); ++i ) {
        size_t index = AMP::Utilities::findfirst( nodelist_ids, node_iterator->globalID() );
        AMP_ASSERT( nodelist_ids[index] == node_iterator->globalID() );
        auto elem_coord = node_iterator->coord();
        for ( int j = 0; j < d_dim; ++j )
            coord[j][index] = elem_coord[j];
        ++node_iterator;
    }
    PROFILE_STOP( "writeMesh - get-nodelist-2", 3 );
    PROFILE_START( "writeMesh - get-nodelist-3", 3 );
    elem_iterator = mesh->getIterator( mesh->getGeomType(), 0 );
    std::vector<int> nodelist;
    nodelist.reserve( shapesize * elem_iterator.size() );
    std::vector<AMP::Mesh::MeshElementID> nodeids;
    for ( const auto &elem : elem_iterator ) {
        elem.getElementsID( AMP::Mesh::GeomType::Vertex, nodeids );
        AMP_INSIST( (int) nodeids.size() == shapesize,
                    "Mixed element types is currently not supported" );
        for ( auto &nodeid : nodeids ) {
            size_t index = AMP::Utilities::findfirst( nodelist_ids, nodeid );
            AMP_ASSERT( nodelist_ids[index] == nodeid );
            nodelist.push_back( (int) index );
        }
        ++elem_iterator;
    }
    PROFILE_STOP( "writeMesh - get-nodelist-3", 3 );
    PROFILE_STOP( "writeMesh - get-nodelist", 2 );
    // Create the directory for the mesh
    PROFILE_START( "writeMesh - directory", 2 );
    std::string tmp_path = data.path;
    while ( tmp_path.size() > 0 ) {
        if ( tmp_path[0] == '/' ) {
            tmp_path.erase( 0, 1 );
            continue;
        }
        size_t pos = tmp_path.find_first_of( '/' );
        if ( pos == std::string::npos ) {
            pos = tmp_path.size();
        }
        auto subdir       = tmp_path.substr( 0, pos );
        DBtoc *toc        = DBGetToc( FileHandle );
        bool subdir_found = false;
        for ( int i = 0; i < toc->ndir; ++i ) {
            if ( subdir.compare( toc->dir_names[i] ) == 0 )
                subdir_found = true;
        }
        if ( !subdir_found )
            DBMkDir( FileHandle, subdir.c_str() );
        DBSetDir( FileHandle, subdir.c_str() );
        tmp_path.erase( 0, pos );
    }
    DBSetDir( FileHandle, "/" );
    DBSetDir( FileHandle, data.path.c_str() );
    PROFILE_STOP( "writeMesh - directory", 2 );
    // Write the elements (connectivity)
    PROFILE_START( "writeMesh - elements", 2 );
    std::string meshName  = data.meshName;
    std::string zoneName  = "zone_" + std::to_string( data.rank );
    auto element_iterator = mesh->getIterator( mesh->getGeomType(), 0 );
    auto num_elems        = (int) element_iterator.size();
    DBPutZonelist2( FileHandle,
                    zoneName.c_str(),
                    num_elems,
                    d_dim,
                    &nodelist[0],
                    nodelist.size(),
                    0,
                    0,
                    0,
                    &shapetype,
                    &shapesize,
                    &shapecnt,
                    1,
                    nullptr );
    PROFILE_STOP( "writeMesh - elements", 2 );
    // Write the mesh
    PROFILE_START( "writeMesh - mesh", 2 );
    DBPutUcdmesh( FileHandle,
                  meshName.c_str(),
                  d_dim,
                  nullptr,
                  coord,
                  node_iterator.size(),
                  nodelist.size(),
                  zoneName.c_str(),
                  nullptr,
                  DB_DOUBLE,
                  nullptr );
    for ( int i = 0; i < d_dim; ++i )
        delete[] coord[i];
    PROFILE_STOP( "writeMesh - mesh", 2 );
    // Write the variables
    PROFILE_START( "writeMesh - variables", 2 );
#ifdef USE_AMP_VECTORS
    float ftime        = time;
    DBoptlist *optlist = DBMakeOptlist( 10 );
    DBAddOption( optlist, DBOPT_CYCLE, &cycle );
    DBAddOption( optlist, DBOPT_TIME, &ftime );
    DBAddOption( optlist, DBOPT_DTIME, &time );
    // DBAddOption(optlist, DBOPT_UNITS, (void *)units);
    for ( size_t i = 0; i < data.varName.size(); ++i ) {
        auto DOFs     = data.vec[i]->getDOFManager();
        int nvar      = 0;
        int centering = 0;
        auto var      = new double *[data.varSize[i]];
        for ( int j = 0; j < data.varSize[i]; ++j )
            var[j] = nullptr;
        const char *varnames[] = { "1", "2", "3" };
        if ( data.varType[i] > mesh->getGeomType() ) {
            // We have a mixed mesh type and there will be no data of the given type for this mesh
        } else if ( data.varType[i] == AMP::Mesh::GeomType::Vertex ) {
            // We are saving node-centered data
            centering = DB_NODECENT;
            nvar      = (int) nodelist_ids.size();
            for ( int j = 0; j < data.varSize[i]; ++j )
                var[j] = new double[nvar];
            std::vector<size_t> dofs( data.varSize[i] );
            std::vector<double> vals( data.varSize[i] );
            for ( int j = 0; j < nvar; ++j ) {
                DOFs->getDOFs( nodelist_ids[j], dofs );
                AMP_ASSERT( (int) dofs.size() == data.varSize[i] );
                data.vec[i]->getValuesByGlobalID( data.varSize[i], &dofs[0], &vals[0] );
                for ( int k = 0; k < data.varSize[i]; ++k )
                    var[k][j] = vals[k];
            }
        } else if ( data.varType[i] == mesh->getGeomType() ) {
            // We are saving cell-centered data
            centering = DB_ZONECENT;
            nvar      = (int) num_elems;
            for ( int j = 0; j < data.varSize[i]; ++j )
                var[j] = new double[nvar];
            std::vector<size_t> dofs( data.varSize[i] );
            std::vector<double> vals( data.varSize[i] );
            auto it = element_iterator.begin();
            for ( int j = 0; j < nvar; ++j, ++it ) {
                DOFs->getDOFs( it->globalID(), dofs );
                data.vec[i]->getValuesByGlobalID( data.varSize[i], &dofs[0], &vals[0] );
                for ( int k = 0; k < data.varSize[i]; ++k )
                    var[k][j] = vals[k];
            }
        } else {
            // We are storing edge or face data
            AMP_ERROR( "The silo writer currently only supports GeomType::Vertex and Cell data" );
        }
        std::string varNameRank = data.varName[i] + "P" + std::to_string( data.rank );
        if ( data.varSize[i] == 1 || data.varSize[i] == d_dim ||
             data.varSize[i] == d_dim * d_dim ) {
            // We are writing a scalar, vector, or tensor variable
            DBPutUcdvar( FileHandle,
                         varNameRank.c_str(),
                         meshName.c_str(),
                         data.varSize[i],
                         (char **) varnames,
                         var,
                         nvar,
                         nullptr,
                         0,
                         DB_DOUBLE,
                         centering,
                         optlist );
        } else {
            // Write each component
            for ( int j = 0; j < data.varSize[i]; ++j ) {
                auto vname = varNameRank + "_" + std::to_string( j );
                DBPutUcdvar( FileHandle,
                             vname.c_str(),
                             meshName.c_str(),
                             1,
                             (char **) varnames,
                             &var[j],
                             nvar,
                             nullptr,
                             0,
                             DB_DOUBLE,
                             centering,
                             optlist );
            }
        }
        for ( int j = 0; j < data.varSize[i]; ++j ) {
            if ( var[j] != nullptr )
                delete[] var[j];
        }
        delete[] var;
    }
    DBFreeOptlist( optlist );
    PROFILE_STOP( "writeMesh - variables", 2 );
#endif
    // Change the directory back to root
    DBSetDir( FileHandle, "/" );
    PROFILE_STOP( "writeMesh", 1 );
}


/************************************************************
 * Function to syncronize the multimesh data                 *
 * If root==-1, the data will be synced across all procs     *
 ************************************************************/
void SiloIO::syncMultiMeshData( std::map<AMP::Mesh::MeshID, siloMultiMeshData> &data,
                                int root ) const
{
    PROFILE_START( "syncMultiMeshData", 1 );
    // Convert the data to vectors
    std::vector<AMP::Mesh::MeshID> ids;
    std::vector<siloMultiMeshData> meshdata;
    int myRank = d_comm.getRank();
    for ( const auto &it : data ) {
        // Only send the base meshes that I own
        auto tmp = it.second;
        tmp.meshes.resize( 0 );
        for ( auto &elem : it.second.meshes ) {
            if ( elem.ownerRank == myRank )
                tmp.meshes.push_back( elem );
        }
        // Only the owner rank will send the variable list
        if ( tmp.ownerRank != myRank )
            tmp.varName.resize( 0 );
        // Only send the multimesh if there are base meshes that need to be sent or I own the mesh
        if ( !tmp.meshes.empty() || tmp.ownerRank == myRank ) {
            ids.push_back( it.first );
            meshdata.push_back( it.second );
        }
    }
    // Create buffers to store the data
    size_t send_size = 0;
    for ( size_t i = 0; i < meshdata.size(); ++i ) {
        AMP_ASSERT( ids[i] == meshdata[i].id );
        send_size += meshdata[i].size();
    }
    auto send_buf = new char[send_size];
    char *ptr     = send_buf;
    for ( auto &elem : meshdata ) {
        elem.pack( ptr );
        ptr = &ptr[elem.size()];
    }
    // Send the data and unpack the buffer to a vector
    size_t tot_num = d_comm.sumReduce( meshdata.size() );
    if ( root == -1 ) {
        // Everybody gets a copy
        size_t tot_size = d_comm.sumReduce( send_size );
        auto recv_buf   = new char[tot_size];
        meshdata.resize( tot_num );
        d_comm.allGather( send_buf, send_size, recv_buf );
        ptr = recv_buf;
        for ( size_t i = 0; i < tot_num; ++i ) {
            meshdata[i] = siloMultiMeshData::unpack( ptr );
            ptr         = &ptr[meshdata[i].size()];
        }
        delete[] recv_buf;
    } else {
        AMP_ASSERT( root >= 0 && root < d_comm.getSize() );
        // Only the root gets a copy
        // Note: the root already has his own data
        size_t max_size = d_comm.maxReduce( send_size );
        std::vector<int> recv_num( d_comm.getSize() );
        d_comm.allGather( (int) meshdata.size(), &recv_num[0] );
        if ( root == d_comm.getRank() ) {
            // Recieve all data
            meshdata.resize( 0 );
            meshdata.reserve( tot_num );
            auto recv_buf = new char[max_size];
            for ( int i = 0; i < d_comm.getSize(); ++i ) {
                if ( i == root )
                    continue;
                int recv_size = d_comm.probe( i, 24987 );
                AMP_ASSERT( recv_size <= (int) max_size );
                d_comm.recv( recv_buf, recv_size, i, false, 24987 );
                char *cptr = recv_buf;
                for ( int j = 0; j < recv_num[i]; ++j ) {
                    auto tmp = siloMultiMeshData::unpack( cptr );
                    cptr     = &cptr[tmp.size()];
                    meshdata.push_back( tmp );
                }
            }
            delete[] recv_buf;
        } else {
            // Send my data
            d_comm.send( send_buf, send_size, root, 24987 );
        }
    }
    delete[] send_buf;
    // Add the meshes from other processors (keeping the existing meshes)
    for ( auto &elem : meshdata ) {
        auto iterator = data.find( elem.id );
        if ( iterator == data.end() ) {
            // Add the multimesh
            data.insert( std::make_pair( elem.id, elem ) );
        } else {
            // Add the submeshes
            for ( auto &meshe : elem.meshes ) {
                bool found = false;
                for ( auto &_k : iterator->second.meshes ) {
                    if ( meshe.id == _k.id && meshe.meshName == _k.meshName &&
                         meshe.path == _k.path && meshe.path == _k.file )
                        found = true;
                }
                if ( !found )
                    iterator->second.meshes.push_back( meshe );
            }
            // Add the variables if we don't have them yet
            if ( elem.varName.size() > 0 ) {
                if ( !iterator->second.varName.empty() )
                    AMP_ASSERT( iterator->second.varName.size() == elem.varName.size() );
                iterator->second.varName = elem.varName;
            }
        }
    }
    PROFILE_STOP( "syncMultiMeshData", 1 );
}


/************************************************************
 * Function to syncronize a variable list                    *
 * If root==-1, the data will be synced across all procs     *
 ************************************************************/
void SiloIO::syncVariableList( std::set<std::string> &data_set, int root ) const
{
    PROFILE_START( "syncVariableList", 1 );
    std::vector<std::string> data( data_set.begin(), data_set.end() );
    size_t N_local  = data.size();
    size_t N_global = d_comm.sumReduce( N_local );
    auto size_local = new size_t[N_local];
    for ( size_t i = 0; i < N_local; ++i )
        size_local[i] = data[i].size();
    auto size_global = new size_t[N_global];
    d_comm.allGather( size_local, N_local, size_global );
    size_t tot_size_local = 0;
    for ( size_t i = 0; i < N_local; ++i )
        tot_size_local += size_local[i];
    size_t tot_size_global = 0;
    for ( size_t i = 0; i < N_global; ++i )
        tot_size_global += size_global[i];
    auto send_buf = new char[tot_size_local];
    auto recv_buf = new char[tot_size_global];
    size_t k      = 0;
    for ( size_t i = 0; i < N_local; ++i ) {
        data[i].copy( &send_buf[k], data[i].size(), 0 );
        k += size_local[i];
    }
    if ( root == -1 ) {
        // Everybody gets a copy
        d_comm.allGather( send_buf, tot_size_local, recv_buf );
        k = 0;
        for ( size_t i = 0; i < N_global; ++i ) {
            std::string tmp( &recv_buf[k], size_global[i] );
            data_set.insert( tmp );
            k += size_global[i];
        }
    } else {
        // Only the root gets a copy
        // Note: the root already has his own data
        AMP_ASSERT( root >= 0 && root < d_comm.getSize() );
        std::vector<int> recv_num( d_comm.getSize() );
        d_comm.allGather( (int) N_local, &recv_num[0] );
        if ( root == d_comm.getRank() ) {
            // Recieve all data
            int index = 0;
            for ( int i = 0; i < d_comm.getSize(); ++i ) {
                if ( i == root ) {
                    index += recv_num[i];
                    continue;
                }
                int recv_size = d_comm.probe( i, 24987 );
                d_comm.recv( recv_buf, recv_size, i, false, 24987 );
                k = 0;
                for ( int j = 0; j < recv_num[i]; ++j ) {
                    std::string tmp( &recv_buf[k], size_global[index] );
                    data_set.insert( tmp );
                    k += size_global[index];
                    index++;
                }
                AMP_ASSERT( (int) k == recv_size );
            }
        } else {
            // Send my data
            d_comm.send( send_buf, tot_size_local, root, 24987 );
        }
    }
    delete[] send_buf;
    delete[] recv_buf;
    delete[] size_local;
    delete[] size_global;
    PROFILE_STOP( "syncVariableList", 1 );
}


/************************************************************
 * Function to write the summary data                        *
 ************************************************************/
static inline std::string getFile( const std::string &file, const std::string &root )
{
    AMP_ASSERT( !file.empty() );
    if ( file.compare( 0, root.size(), root ) == 0 )
        return file.substr( root.size() );
    return file;
}
void SiloIO::writeSummary( std::string filename, int cycle, double time )
{
    PROFILE_START( "writeSummary", 1 );
    AMP_ASSERT( !filename.empty() );
    // Add the siloBaseMeshData to the multimeshes
    auto multiMeshes = d_multiMeshes;
    for ( auto &tmp : multiMeshes ) {
        auto mesh     = tmp.second.mesh;
        auto base_ids = getMeshIDs( mesh );
        for ( auto id : base_ids ) {
            const auto &it = d_baseMeshes.find( id );
            if ( it != d_baseMeshes.end() ) {
                siloBaseMeshData data = it->second;
                AMP_ASSERT( it->first == data.id );
                tmp.second.meshes.push_back( data );
            }
        }
    }
    // Add the whole mesh
    /*if ( multiMeshes.size()==0 ) {
        siloMultiMeshData wholemesh;
        wholemesh.id = AMP::Mesh::MeshID((unsigned int)-1,0);
        wholemesh.name = "whole_mesh";
        for (auto it=d_baseMeshes.begin(); it!=d_baseMeshes.end(); ++it) {
            siloBaseMeshData data = it->second;
            AMP_ASSERT(it->first==data.id);
            wholemesh.meshes.push_back(data);
        }
        wholemesh.owner_rank = 0;
        multimeshes.insert( std::make_pair(wholemesh.id,wholemesh)
    );
    }*/
    // Gather the results
    // Note: we only need to guarantee that rank 0 has all the data
    syncMultiMeshData( multiMeshes, 0 );
    syncVariableList( d_varNames, 0 );
    // Write the multimeshes and multivariables
    std::string base_path;
    if ( find_slash( filename ) != std::string::npos )
        base_path = filename.substr( 0, find_slash( filename ) + 1 );
    if ( d_comm.getRank() == 0 ) {
        DBfile *FileHandle = DBOpen( filename.c_str(), DB_HDF5, DB_APPEND );
        // Create the subdirectories
        PROFILE_START( "create directories", 2 );
        std::set<std::string> subdirs;
        for ( const auto &tmp : multiMeshes ) {
            auto data  = tmp.second;
            auto file  = getFile( data.name, base_path );
            size_t pos = find_slash( file );
            if ( pos != std::string::npos )
                subdirs.insert( file.substr( 0, pos ) );
        }
        for ( const auto &subdir : subdirs )
            createSiloDirectory( FileHandle, subdir );
        PROFILE_STOP( "create directories", 2 );
        // Create the multimeshes
        PROFILE_START( "write multimeshes", 2 );
        for ( const auto &tmp : multiMeshes ) {
            const auto &data = tmp.second;
            size_t N         = data.meshes.size();
            std::vector<std::string> meshNames( N );
            for ( size_t i = 0; i < N; ++i ) {
                auto file    = getFile( data.meshes[i].file, base_path );
                meshNames[i] = file + ":" + data.meshes[i].path + "/" + data.meshes[i].meshName;
                strrep( meshNames[i], "//", "/" );
            }
            auto meshnames = new char *[N];
            auto meshtypes = new int[N];
            for ( size_t i = 0; i < N; ++i ) {
                meshnames[i] = (char *) meshNames[i].c_str();
                meshtypes[i] = DB_UCDMESH;
            }
            std::string tree_name = data.name + "_tree";
            DBoptlist *optList    = DBMakeOptlist( 10 );
            DBAddOption( optList, DBOPT_MRGTREE_NAME, (char *) tree_name.c_str() );
            DBPutMultimesh(
                FileHandle, data.name.c_str(), meshNames.size(), meshnames, meshtypes, nullptr );
            DBFreeOptlist( optList );
            delete[] meshnames;
            delete[] meshtypes;
        }
        PROFILE_STOP( "write multimeshes", 2 );
        // Generate the multi-variables
        PROFILE_START( "write multivariables", 2 );
        for ( const auto &tmp : multiMeshes ) {
            const auto &data = tmp.second;
            size_t N         = data.meshes.size();
            // std::cout << data.name << std::endl;
            for ( const auto &varName : data.varName ) {
                std::vector<std::string> varNames( N );
                auto varnames = new char *[N];
                auto vartypes = new int[N];
                for ( size_t i = 0; i < N; ++i ) {
                    std::string rankStr = std::to_string( data.meshes[i].rank );
                    auto file           = getFile( data.meshes[i].file, base_path );
                    varNames[i] = file + ":" + data.meshes[i].path + "/" + varName + "P" + rankStr;
                    strrep( varNames[i], "//", "/" );
                    varnames[i] = (char *) varNames[i].c_str();
                    vartypes[i] = DB_UCDVAR;
                }
                int varSize = 0;
                for ( size_t i = 0; i < data.meshes[0].varName.size(); ++i ) {
                    if ( data.meshes[0].varName[i] == varName ) {
                        varSize = data.meshes[0].varSize[i];
                        break;
                    }
                }
                auto multiMeshName = data.name;
                auto visitVarName  = multiMeshName + "_" + varName;
                float ftime        = time;
                DBoptlist *opts    = DBMakeOptlist( 10 );
                DBAddOption( opts, DBOPT_CYCLE, &cycle );
                DBAddOption( opts, DBOPT_TIME, &ftime );
                DBAddOption( opts, DBOPT_DTIME, &time );
                // DBAddOption( opts, DBOPT_MMESH_NAME, (char*) multiMeshName.c_str() );
                if ( varSize == 1 || varSize == d_dim || varSize == d_dim * d_dim ) {
                    // We are writing a scalar, vector, or tensor variable
                    DBPutMultivar( FileHandle,
                                   visitVarName.c_str(),
                                   varNames.size(),
                                   varnames,
                                   vartypes,
                                   opts );
                } else {
                    // Write each component
                    for ( int j = 0; j < varSize; ++j ) {
                        std::string postfix = "_" + std::to_string( j );
                        std::vector<std::string> varNames2( data.meshes.size() );
                        for ( size_t k = 0; k < data.meshes.size(); ++k ) {
                            varNames2[k] = varNames[k] + postfix;
                            varnames[k]  = (char *) varNames2[k].c_str();
                        }
                        DBPutMultivar( FileHandle,
                                       ( visitVarName + postfix ).c_str(),
                                       varNames.size(),
                                       varnames,
                                       vartypes,
                                       opts );
                    }
                }
                DBFreeOptlist( opts );
                delete[] varnames;
                delete[] vartypes;
            }
        }
        PROFILE_STOP( "write multivariables", 2 );
        DBClose( FileHandle );
    }
    PROFILE_STOP( "writeSummary", 1 );
}


/************************************************************
 * Functions to pack/unpack data to a char array             *
 ************************************************************/
template<class TYPE>
static inline void packData( char *ptr, size_t &pos, const TYPE &data );
template<class TYPE>
static inline TYPE unpackData( const char *ptr, size_t &pos );
template<>
inline void packData<std::string>( char *ptr, size_t &pos, const std::string &data )
{
    int N = data.size();
    memcpy( &ptr[pos], data.c_str(), N + 1 );
    pos += N + 1;
}
template<>
inline std::string unpackData<std::string>( const char *ptr, size_t &pos )
{
    std::string data( &ptr[pos] );
    pos += data.size() + 1;
    return data;
}
template<class TYPE>
static inline void packData( char *ptr, size_t &pos, const TYPE &data )
{
    memcpy( &ptr[pos], &data, sizeof( TYPE ) );
    pos += sizeof( TYPE );
}
template<class TYPE>
static inline TYPE unpackData( const char *ptr, size_t &pos )
{
    TYPE data;
    memcpy( &data, &ptr[pos], sizeof( TYPE ) );
    pos += sizeof( TYPE );
    return data;
}


/************************************************************
 * Functions for siloBaseMeshData                            *
 ************************************************************/
size_t SiloIO::siloBaseMeshData::size() const
{
    size_t N_bytes = sizeof( AMP::Mesh::MeshID ); // Store the mesh id
    N_bytes += sizeof( int );                     // Store the processor rank
    N_bytes += sizeof( int );                     // Store the owner rank
    N_bytes += meshName.size() + 1;               // Store the mesh name
    N_bytes += path.size() + 1;                   // Store the mesh path
    N_bytes += file.size() + 1;                   // Store the mesh file
    N_bytes += sizeof( int );                     // Store the number of variables
    for ( auto &elem : varName ) {
        N_bytes += elem.size() + 1;               // Store the variable name
        N_bytes += sizeof( AMP::Mesh::GeomType ); // Store the variable type
        N_bytes += sizeof( int );                 // Store the number of unknowns per point
    }
    return N_bytes;
}
void SiloIO::siloBaseMeshData::pack( char *ptr ) const
{
    size_t pos = 0;
    packData<AMP::Mesh::MeshID>( ptr, pos, id );
    packData<int>( ptr, pos, rank );
    packData<int>( ptr, pos, ownerRank );
    packData<std::string>( ptr, pos, meshName );
    packData<std::string>( ptr, pos, path );
    packData<std::string>( ptr, pos, file );
    packData<int>( ptr, pos, varName.size() );
    for ( size_t i = 0; i < varName.size(); ++i ) {
        packData<std::string>( ptr, pos, varName[i] );
        packData<AMP::Mesh::GeomType>( ptr, pos, varType[i] );
        packData<int>( ptr, pos, varSize[i] );
    }
    AMP_ASSERT( pos == size() );
}
SiloIO::siloBaseMeshData SiloIO::siloBaseMeshData::unpack( const char *ptr )
{
    siloBaseMeshData data;
    size_t pos     = 0;
    data.id        = unpackData<AMP::Mesh::MeshID>( ptr, pos );
    data.rank      = unpackData<int>( ptr, pos );
    data.ownerRank = unpackData<int>( ptr, pos );
    data.meshName  = unpackData<std::string>( ptr, pos );
    data.path      = unpackData<std::string>( ptr, pos );
    data.file      = unpackData<std::string>( ptr, pos );
    // Store the variables
    size_t N_var = unpackData<int>( ptr, pos );
    data.varName.resize( N_var );
    data.varType.resize( N_var );
    data.varSize.resize( N_var );
#ifdef USE_AMP_VECTORS
    data.vec.resize( N_var );
#endif
    for ( size_t i = 0; i < N_var; ++i ) {
        data.varName[i] = unpackData<std::string>( ptr, pos );
        data.varType[i] = unpackData<AMP::Mesh::GeomType>( ptr, pos );
        data.varSize[i] = unpackData<int>( ptr, pos );
    }
    AMP_ASSERT( pos == data.size() );
    return data;
}


/************************************************************
 * Functions for siloMultiMeshData                           *
 ************************************************************/
SiloIO::siloMultiMeshData::siloMultiMeshData( const SiloIO::siloMultiMeshData &rhs )
    : id( rhs.id ),
      mesh( rhs.mesh ),
      ownerRank( rhs.ownerRank ),
      name( rhs.name ),
      meshes( rhs.meshes ),
      varName( rhs.varName )
{
}
SiloIO::siloMultiMeshData &
SiloIO::siloMultiMeshData::operator=( const SiloIO::siloMultiMeshData &rhs )
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return *this;
    this->id        = rhs.id;
    this->mesh      = rhs.mesh;
    this->ownerRank = rhs.ownerRank;
    this->name      = rhs.name;
    this->meshes    = rhs.meshes;
    this->varName   = rhs.varName;
    return *this;
}
size_t SiloIO::siloMultiMeshData::size() const
{
    size_t N_bytes = sizeof( AMP::Mesh::MeshID ); // Store the mesh id
    N_bytes += sizeof( int );                     // Store the owner rank
    N_bytes += name.size() + 1;                   // Store the mesh name
    N_bytes += sizeof( int );                     // Store the number of sub meshes
    for ( const auto &mesh : meshes )
        N_bytes += mesh.size(); // Store the sub meshes
    N_bytes += sizeof( int );   // Store the number of variables
    for ( const auto &name : varName )
        N_bytes += name.size() + 1; // Store the variable name
    return N_bytes;
}
void SiloIO::siloMultiMeshData::pack( char *ptr ) const
{
    size_t pos = 0;
    packData<AMP::Mesh::MeshID>( ptr, pos, id );
    packData<int>( ptr, pos, ownerRank );
    packData<std::string>( ptr, pos, name );
    // Store the base meshes
    packData<int>( ptr, pos, meshes.size() );
    for ( const auto &mesh : meshes ) {
        mesh.pack( &ptr[pos] );
        pos += mesh.size();
    }
    // Store the variables
    packData<int>( ptr, pos, varName.size() );
    for ( const auto &name : varName )
        packData<std::string>( ptr, pos, name );
    AMP_ASSERT( pos == size() );
}
SiloIO::siloMultiMeshData SiloIO::siloMultiMeshData::unpack( const char *ptr )
{
    size_t pos = 0;
    siloMultiMeshData data;
    data.id        = unpackData<AMP::Mesh::MeshID>( ptr, pos );
    data.ownerRank = unpackData<int>( ptr, pos );
    data.name      = unpackData<std::string>( ptr, pos );
    // Store the base meshes
    int N_meshes = unpackData<int>( ptr, pos );
    data.meshes.resize( N_meshes );
    for ( int i = 0; i < N_meshes; ++i ) {
        data.meshes[i] = siloBaseMeshData::unpack( &ptr[pos] );
        pos += data.meshes[i].size();
    }
    // Store the variables
    int N_var    = unpackData<int>( ptr, pos );
    data.varName = std::vector<std::string>( N_var );
    for ( auto &name : data.varName )
        name = unpackData<std::string>( ptr, pos );
    AMP_ASSERT( pos == data.size() );
    return data;
}


/************************************************************
 * Some utilit functions                                     *
 ************************************************************/
void createSiloDirectory( DBfile *FileHandle, const std::string &path )
{
    // Create a subdirectory tree from the current working path if it does not exist
    char current_dir[256];
    DBGetDir( FileHandle, current_dir );
    // Get the list of directories that may need to be created
    std::vector<std::string> subdirs;
    std::string path2 = path + "/";
    while ( !path2.empty() ) {
        size_t pos = path2.find( "/" );
        if ( pos > 0 ) {
            subdirs.push_back( path2.substr( 0, pos ) );
        }
        path2.erase( 0, pos + 1 );
    }
    // Create the directories as necessary
    for ( auto &subdir : subdirs ) {
        DBtoc *toc  = DBGetToc( FileHandle );
        bool exists = false;
        for ( int j = 0; j < toc->ndir; ++j ) {
            if ( subdir.compare( toc->dir_names[j] ) == 0 )
                exists = true;
        }
        if ( !exists )
            DBMkDir( FileHandle, subdir.c_str() );
        DBSetDir( FileHandle, subdir.c_str() );
    }
    // Return back to the original working directory
    DBSetDir( FileHandle, current_dir );
}


#else
void SiloIO::readFile( const std::string & ) { AMP_ERROR( "SILO not configured" ); }
void SiloIO::writeFile( const std::string &, size_t, double )
{
    AMP_ERROR( "SILO not configured" );
}
void SiloIO::registerMesh( AMP::Mesh::Mesh::shared_ptr, int, const std::string & )
{
    AMP_ERROR( "SILO not configured" );
}
#ifdef USE_AMP_VECTORS
void SiloIO::registerVector( AMP::LinearAlgebra::Vector::shared_ptr,
                             AMP::Mesh::Mesh::shared_ptr,
                             AMP::Mesh::GeomType,
                             const std::string & )
{
    AMP_ERROR( "SILO not configured" );
}
void SiloIO::registerVector( AMP::LinearAlgebra::Vector::shared_ptr, const std::string & )
{
    AMP_ERROR( "SILO not configured" );
}
#endif
#ifdef USE_AMP_MATRICES
void SiloIO::registerMatrix( AMP::LinearAlgebra::Matrix::shared_ptr, const std::string & )
{
    AMP_ERROR( "SILO not configured" );
}
#endif


#endif


} // namespace AMP::Utilities

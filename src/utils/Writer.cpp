#include "AMP/utils/Writer.h"
#include "AMP/utils/Utilities.h"

#include "AMP/utils/AMP_MPI.I"
#include "AMP/utils/AsciiWriter.h"
#include "AMP/utils/HDF5writer.h"
#include "AMP/utils/NullWriter.h"
#include "AMP/utils/SiloWriter.h"

#ifdef USE_AMP_MESH
#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MultiMesh.h"
#endif
#ifdef USE_AMP_VECTORS
#include "AMP/vectors/Vector.h"
#endif
#ifdef USE_AMP_MATRICES
#include "AMP/matrices/Matrix.h"
#endif


namespace AMP::Utilities {


/************************************************************
 * Functions to pack/unpack data to a char array             *
 ************************************************************/
template<class TYPE>
static inline void packData( char *ptr, size_t &pos, const TYPE &data )
{
    if constexpr ( std::is_trivially_copyable<TYPE>::value ) {
        memcpy( &ptr[pos], &data, sizeof( TYPE ) );
        pos += sizeof( TYPE );
    } else if constexpr ( std::is_same<TYPE, std::string>::value ) {
        int N = data.size();
        memcpy( &ptr[pos], data.c_str(), N + 1 );
        pos += N + 1;
    } else {
        throw std::logic_error( "Error packing data" );
    }
}
template<class TYPE>
static inline TYPE unpackData( const char *ptr, size_t &pos )
{
    if constexpr ( std::is_trivially_copyable<TYPE>::value ) {
        TYPE data;
        memcpy( &data, &ptr[pos], sizeof( TYPE ) );
        pos += sizeof( TYPE );
        return data;
    } else if constexpr ( std::is_same<TYPE, std::string>::value ) {
        std::string data( &ptr[pos] );
        pos += data.size() + 1;
        return data;
    } else {
        throw std::logic_error( "Error packing data" );
    }
}


/************************************************************
 * Helper function to get a unique id for each vector        *
 ************************************************************/
Writer::GlobalID Writer::getID( const AMP_MPI &comm ) const
{
    GlobalID id;
    if ( comm.getRank() == 0 ) {
        static uint64_t localID = 0;
        id.objID                = localID++;
        id.ownerRank            = d_comm.getRank();
    }
    return comm.bcast( id, 0 );
}


/************************************************************
 * WriterProperties                                          *
 ************************************************************/
Writer::WriterProperties::WriterProperties()
    : registerMesh( false ),
      registerVector( false ),
      registerVectorWithMesh( false ),
      registerMatrix( false )
{
}


/************************************************************
 * Writer::VectorData                                        *
 ************************************************************/
Writer::VectorData::VectorData( std::shared_ptr<AMP::LinearAlgebra::Vector> vec_,
                                const std::string &name_ )
    : name( name_ ), numDOFs( 0 ), vec( vec_ ), type( static_cast<AMP::Mesh::GeomType>( 0xFF ) )
{
    if ( !vec )
        return;
    if ( name.empty() )
        name = vec->getVariable()->getName();
    if ( name.empty() )
        name = vec->type();
}


/************************************************************
 * Writer::MatrixData                                        *
 ************************************************************/
Writer::MatrixData::MatrixData( std::shared_ptr<AMP::LinearAlgebra::Matrix> mat_,
                                const std::string &name_ )
    : name( name_ ), mat( mat_ )
{
    if ( !mat )
        return;
    if ( name.empty() )
        name = mat->getLeftVector()->getVariable()->getName() + " - " +
               mat->getRightVector()->getVariable()->getName();
}


/************************************************************
 * Builder                                                   *
 ************************************************************/
std::shared_ptr<AMP::Utilities::Writer> Writer::buildWriter( std::string type, AMP_MPI comm )
{
    std::for_each( type.begin(), type.end(), []( char &c ) { c = ::tolower( c ); } );
    std::shared_ptr<AMP::Utilities::Writer> writer;
    if ( type == "none" ) {
        writer.reset( new AMP::Utilities::NullWriter() );
    } else if ( type == "silo" ) {
        writer.reset( new AMP::Utilities::SiloIO() );
    } else if ( type == "hdf5" ) {
        writer.reset( new AMP::Utilities::HDF5writer() );
    } else if ( type == "ascii" ) {
        writer.reset( new AMP::Utilities::AsciiWriter() );
    } else {
        AMP_ERROR( "Unknown writer: " + type );
    }
    writer->d_comm = std::move( comm );
    return writer;
}
std::shared_ptr<AMP::Utilities::Writer> Writer::buildWriter( std::shared_ptr<AMP::Database> db )
{
    auto type   = db->getString( "Name" );
    auto writer = Writer::buildWriter( type );
    if ( db->keyExists( "Decomposition" ) )
        writer->setDecomposition( db->getScalar<int>( "Decomposition" ) );
    return writer;
}


/************************************************************
 * Constructor/Destructor                                    *
 ************************************************************/
Writer::Writer() : d_comm( AMP_COMM_WORLD ) { d_decomposition = 2; }
Writer::~Writer() = default;


/************************************************************
 * Some basic functions                                      *
 ************************************************************/
std::string Writer::getExtension() const { return getProperties().extension; }
void Writer::setDecomposition( int d )
{
    AMP_INSIST( d == 1 || d == 2, "decomposition must be 1 or 2" );
    d_decomposition = d;
}
void Writer::createDirectories( const std::string &filename )
{
    size_t i = filename.rfind( '/' );
    if ( i != std::string::npos && d_comm.getRank() == 0 )
        AMP::Utilities::recursiveMkdir(
            filename.substr( 0, i ), ( S_IRUSR | S_IWUSR | S_IXUSR ), false );
    d_comm.barrier();
}


/************************************************************
 * Register a mesh                                           *
 ************************************************************/
void Writer::registerMesh( AMP::Mesh::Mesh::shared_ptr mesh, int level, const std::string &path )
{
    AMP_INSIST( level >= 0 && level <= 3, "Invalid value for level" );
    std::set<GlobalID> base_ids;
    registerMesh2( mesh, level, path, base_ids );
}
void Writer::registerMesh2( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                            int level,
                            const std::string &path,
                            std::set<GlobalID> &base_ids )
{
    if ( !mesh )
        return;
    if ( !getProperties().registerMesh )
        AMP_ERROR( "registerMesh is not supported for " + getProperties().type );
    auto multimesh = std::dynamic_pointer_cast<AMP::Mesh::MultiMesh>( mesh );
    if ( !multimesh ) {
        // Create a unique id for each rank
        GlobalID id( mesh->meshID().getData(), d_comm.getRank() );
        base_ids.insert( id );
        if ( d_baseMeshes.find( id ) != d_baseMeshes.end() )
            return;
        // We are dealing with a single mesh
        int rank = d_comm.getRank();
        baseMeshData data;
        data.id        = id;
        data.mesh      = mesh;
        data.rank      = rank;
        data.ownerRank = rank;
        data.meshName  = "rank_" + std::to_string( data.rank );
        data.path      = path + mesh->getName() + "_/";
        if ( d_baseMeshes.find( id ) == d_baseMeshes.end() )
            d_baseMeshes.insert( std::make_pair( id, data ) );
        // Create and register a multimesh for the current mesh
        if ( level > -1 ) {
            multiMeshData data2;
            data2.id        = GlobalID( mesh->meshID().getData(), 0 );
            data2.mesh      = mesh;
            data2.name      = path + mesh->getName();
            data2.meshes    = { id };
            data2.ownerRank = mesh->getComm().bcast<int>( rank, 0 );
            d_multiMeshes.insert( std::make_pair( data2.id, data2 ) );
        }
        // Create and register a multimesh for the rank
        if ( level == 3 ) {
            // Create a unique id for each rank
            multiMeshData data2;
            data2.id        = id;
            data2.mesh      = mesh;
            data2.name      = path + mesh->getName() + "_/rank_" + std::to_string( data.rank );
            data2.meshes    = { id };
            data2.ownerRank = rank;
            d_multiMeshes.insert( std::make_pair( data2.id, data2 ) );
        }
    } else {
        // We are dealing with a multimesh, register the current mesh and sub meshes
        GlobalID id( mesh->meshID().getData(), 0 );
        int level2 = level;
        if ( level == 1 )
            level2 = -1;
        auto new_path  = path + mesh->getName() + "_/";
        auto submeshes = multimesh->getMeshes();
        std::set<GlobalID> ids;
        for ( auto &submesh : submeshes )
            registerMesh2( submesh, level2, new_path, ids );
        base_ids.insert( ids.begin(), ids.end() );
        if ( level > 0 ) {
            multiMeshData data;
            data.id        = id;
            data.mesh      = mesh;
            data.name      = path + mesh->getName();
            data.ownerRank = mesh->getComm().bcast<int>( d_comm.getRank(), 0 );
            data.meshes.insert( data.meshes.begin(), ids.begin(), ids.end() );
            d_multiMeshes.insert( std::make_pair( data.id, data ) );
        }
    }
}


/************************************************************
 * Register a vector without a mesh                          *
 ************************************************************/
void Writer::registerVector( std::shared_ptr<AMP::LinearAlgebra::Vector> vec,
                             const std::string &name )
{
#ifdef USE_AMP_VECTORS
    auto id = getID( vec->getComm() );
    VectorData data( vec, name );
    d_vectors.insert( std::make_pair( id, std::move( data ) ) );
#else
    NULL_USE( vec );
    NULL_USE( name );
#endif
}


/************************************************************
 * Register a vector with a mesh                             *
 ************************************************************/
#ifdef USE_AMP_VECTORS
void Writer::registerVector( AMP::LinearAlgebra::Vector::shared_ptr vec,
                             AMP::Mesh::Mesh::shared_ptr mesh,
                             AMP::Mesh::GeomType type,
                             const std::string &name_in )
{
    // Return if the vector or mesh is empty
    if ( !vec || !mesh )
        return;
    // Make sure the mesh has been registered
    std::set<GlobalID> base_ids;
    registerMesh2( mesh, -1, "", base_ids );
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
    VectorData data( vec, name_in );
    data.type    = type;
    data.numDOFs = DOFsPerPoint;
    auto ids     = getMeshIDs( mesh );
    for ( auto id : ids ) {
        for ( auto &[id0, mesh2] : d_baseMeshes ) {
            if ( id0.objID == id.getData() )
                mesh2.vectors.push_back( data );
        }
    }
    // Register the vector with the appropriate multi-meshes
    for ( auto &[id0, mesh2] : d_multiMeshes ) {
        if ( id0.objID == mesh->meshID().getData() )
            mesh2.varName.push_back( data.name );
    }
    // Add the vector to the list of vectors so we can perform makeConsistent
    d_vectorsMesh.push_back( vec );
}
#endif


/************************************************************
 * Register a matrix                                         *
 ************************************************************/
void Writer::registerMatrix( std::shared_ptr<AMP::LinearAlgebra::Matrix> mat,
                             const std::string &name )
{
#ifdef USE_AMP_MATRICES
    auto id = getID( mat->getLeftVector()->getComm() );
    MatrixData data( mat, name );
    d_matrices.insert( std::make_pair( id, std::move( data ) ) );
#else
    NULL_USE( mat );
    NULL_USE( name );
#endif
}


/************************************************************
 * Function to get the mesh ids to use for registering       *
 ************************************************************/
std::vector<AMP::Mesh::MeshID> Writer::getMeshIDs( std::shared_ptr<AMP::Mesh::Mesh> mesh )
{
#ifdef USE_AMP_MESH
    std::vector<AMP::Mesh::MeshID> ids;
    auto multimesh = std::dynamic_pointer_cast<AMP::Mesh::MultiMesh>( mesh );
    if ( !multimesh ) {
        // We are dealing with a single mesh
        ids = std::vector<AMP::Mesh::MeshID>( 1, mesh->meshID() );
    } else {
        // We are dealing with a multimesh
        auto meshes = multimesh->getMeshes();
        for ( auto &meshe : meshes ) {
            auto ids2 = getMeshIDs( meshe );
            for ( auto &elem : ids2 )
                ids.push_back( elem );
        }
    }
    return ids;
#else
    return std::vector<AMP::Mesh::MeshID>();
#endif
}


/************************************************************
 * Functions for baseMeshData                                *
 ************************************************************/
size_t Writer::baseMeshData::size() const
{
    size_t N_bytes = sizeof( id );  // Store the mesh id
    N_bytes += sizeof( int );       // Store the processor rank
    N_bytes += sizeof( int );       // Store the owner rank
    N_bytes += meshName.size() + 1; // Store the mesh name
    N_bytes += path.size() + 1;     // Store the mesh path
    N_bytes += file.size() + 1;     // Store the mesh file
    N_bytes += sizeof( int );       // Store the number of variables
    for ( auto &vec : vectors ) {
        N_bytes += vec.name.size() + 1;           // Store the variable name
        N_bytes += sizeof( AMP::Mesh::GeomType ); // Store the variable type
        N_bytes += sizeof( int );                 // Store the number of unknowns per point
    }
    return N_bytes;
}
void Writer::baseMeshData::pack( char *ptr ) const
{
    size_t pos = 0;
    packData<GlobalID>( ptr, pos, id );
    packData<int>( ptr, pos, rank );
    packData<int>( ptr, pos, ownerRank );
    packData<std::string>( ptr, pos, meshName );
    packData<std::string>( ptr, pos, path );
    packData<std::string>( ptr, pos, file );
    packData<int>( ptr, pos, vectors.size() );
    for ( auto &vec : vectors ) {
        packData<std::string>( ptr, pos, vec.name );
        packData<AMP::Mesh::GeomType>( ptr, pos, vec.type );
        packData<int>( ptr, pos, vec.numDOFs );
    }
    AMP_ASSERT( pos == size() );
}
Writer::baseMeshData Writer::baseMeshData::unpack( const char *ptr )
{
    baseMeshData data;
    size_t pos     = 0;
    data.id        = unpackData<GlobalID>( ptr, pos );
    data.rank      = unpackData<int>( ptr, pos );
    data.ownerRank = unpackData<int>( ptr, pos );
    data.meshName  = unpackData<std::string>( ptr, pos );
    data.path      = unpackData<std::string>( ptr, pos );
    data.file      = unpackData<std::string>( ptr, pos );
    // Store the variables
    size_t N_var = unpackData<int>( ptr, pos );
    data.vectors.resize( N_var );
    for ( size_t i = 0; i < N_var; ++i ) {
        data.vectors[i].name    = unpackData<std::string>( ptr, pos );
        data.vectors[i].type    = unpackData<AMP::Mesh::GeomType>( ptr, pos );
        data.vectors[i].numDOFs = unpackData<int>( ptr, pos );
    }
    AMP_ASSERT( pos == data.size() );
    return data;
}


/************************************************************
 * Functions for multiMeshData                               *
 ************************************************************/
Writer::multiMeshData::multiMeshData( const Writer::multiMeshData &rhs )
    : id( rhs.id ),
      mesh( rhs.mesh ),
      ownerRank( rhs.ownerRank ),
      name( rhs.name ),
      meshes( rhs.meshes ),
      varName( rhs.varName )
{
}
Writer::multiMeshData &Writer::multiMeshData::operator=( const Writer::multiMeshData &rhs )
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
size_t Writer::multiMeshData::size() const
{
    size_t N_bytes = sizeof( uint64_t );           // Store the mesh id
    N_bytes += sizeof( int );                      // Store the owner rank
    N_bytes += name.size() + 1;                    // Store the mesh name
    N_bytes += sizeof( int );                      // Store the number of sub meshes
    N_bytes += meshes.size() * sizeof( GlobalID ); // Store the sub meshes
    N_bytes += sizeof( int );                      // Store the number of variables
    for ( const auto &name : varName )
        N_bytes += name.size() + 1; // Store the variable name
    return N_bytes;
}
void Writer::multiMeshData::pack( char *ptr ) const
{
    /*size_t pos = 0;
    packData<GlobalID>( ptr, pos, id );
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
    AMP_ASSERT( pos == size() );*/
}
Writer::multiMeshData Writer::multiMeshData::unpack( const char *ptr )
{
    /*size_t pos = 0;
    multiMeshData data;
    data.id        = unpackData<GlobalID>( ptr, pos );
    data.ownerRank = unpackData<int>( ptr, pos );
    data.name      = unpackData<std::string>( ptr, pos );
    // Store the base meshes
    int N_meshes = unpackData<int>( ptr, pos );
    data.meshes.resize( N_meshes );
    for ( int i = 0; i < N_meshes; ++i ) {
        data.meshes[i] = baseMeshData::unpack( &ptr[pos] );
        pos += data.meshes[i].size();
    }
    // Store the variables
    int N_var    = unpackData<int>( ptr, pos );
    data.varName = std::vector<std::string>( N_var );
    for ( auto &name : data.varName )
        name = unpackData<std::string>( ptr, pos );
    AMP_ASSERT( pos == data.size() );
    return data;*/
}


/************************************************************
 * Function to synchronize the multimesh data                *
 * If root==-1, the data will be synced across all procs     *
 ************************************************************/
void Writer::syncMultiMeshData( std::map<GlobalID, multiMeshData> &data, int root ) const
{
    if ( d_comm.getSize() == 1 )
        return;
    // Convert the data to vectors
    std::vector<GlobalID> ids;
    std::vector<multiMeshData> meshdata;
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
            meshdata[i] = multiMeshData::unpack( ptr );
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
                    auto tmp = multiMeshData::unpack( cptr );
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
    /*    for ( auto &multimesh : meshdata ) {
            auto iterator = data.find( multimesh.id );
            if ( iterator == data.end() ) {
                // Add the multimesh
                data.insert( std::make_pair( multimesh.id, multimesh ) );
            } else {
                // Add the submeshes
                for ( auto &id : multimesh.meshes ) {
                    bool found = false;
                    for ( auto &id2 : iterator->second.meshes ) {
                        if ( id == _k.id && meshe.meshName == _k.meshName &&
                             meshe.path == _k.path && meshe.path == _k.file )
                            found = true;
                    }
                    if ( !found )
                        iterator->second.meshes.push_back( meshe );
                }
                // Add the variables if we don't have them yet
                if ( multimesh.varName.size() > 0 ) {
                    if ( !iterator->second.varName.empty() )
                        AMP_ASSERT( iterator->second.varName.size() == multimesh.varName.size() );
                    iterator->second.varName = multimesh.varName;
                }
            }
        }*/
    AMP_ERROR( "Not finished" );
}


/************************************************************
 * Helper function to get node and element lists for a mesh  *
 ************************************************************/
void Writer::getNodeElemList( std::shared_ptr<const AMP::Mesh::Mesh> mesh,
                              const AMP::Mesh::MeshIterator &elements,
                              AMP::Array<double> *x,
                              AMP::Array<int> &nodelist,
                              std::vector<AMP::Mesh::MeshElementID> &nodelist_ids )
{
    AMP_ASSERT( elements.size() > 0 );
    int ndim = mesh->getDim();
    // Get the element list
    auto elem_iterator = elements.begin();
    auto nodes         = elem_iterator->getElements( AMP::Mesh::GeomType::Vertex );
    int shapesize      = nodes.size();
    // Get the node list (unique integer for each node) and coordinates
    auto node_iterator = mesh->getIterator( AMP::Mesh::GeomType::Vertex, 1 );
    nodelist_ids.resize( node_iterator.size() );
    for ( size_t i = 0; i < node_iterator.size(); ++i, ++node_iterator )
        nodelist_ids[i] = node_iterator->globalID();
    AMP::Utilities::quicksort( nodelist_ids );
    double *coord[3] = { nullptr };
    for ( int d = 0; d < ndim; ++d ) {
        x[d].resize( node_iterator.size() );
        coord[d] = x[d].data();
    }
    node_iterator = mesh->getIterator( AMP::Mesh::GeomType::Vertex, 1 );
    for ( size_t i = 0; i < node_iterator.size(); ++i ) {
        size_t index = AMP::Utilities::findfirst( nodelist_ids, node_iterator->globalID() );
        AMP_ASSERT( nodelist_ids[index] == node_iterator->globalID() );
        auto elem_coord = node_iterator->coord();
        for ( int d = 0; d < ndim; ++d )
            coord[d][index] = elem_coord[d];
        ++node_iterator;
    }
    elem_iterator = elements.begin();
    nodelist.resize( shapesize, elem_iterator.size() );
    std::vector<AMP::Mesh::MeshElementID> nodeids;
    size_t i = 0;
    for ( const auto &elem : elem_iterator ) {
        elem.getElementsID( AMP::Mesh::GeomType::Vertex, nodeids );
        AMP_INSIST( (int) nodeids.size() == shapesize,
                    "Mixed element types is currently not supported" );
        for ( auto &nodeid : nodeids ) {
            int index = AMP::Utilities::findfirst( nodelist_ids, nodeid );
            AMP_ASSERT( nodelist_ids[index] == nodeid );
            nodelist( i++ ) = index;
        }
        ++elem_iterator;
    }
    AMP_ASSERT( i == nodelist.length() );
}


} // namespace AMP::Utilities

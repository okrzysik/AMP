#include "AMP/utils/Writer.h"
#include "AMP/utils/Utilities.h"

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
 * Helper function to get a unique id for each vector        *
 ************************************************************/
static uint32_t localID = 0;
uint64_t getID( const AMP_MPI &local_comm, const AMP_MPI &global_comm )
{
    uint64_t id = 0;
    if ( local_comm.getRank() == 0 ) {
        uint64_t id1 = static_cast<uint64_t>( global_comm.getRank() ) << 32;
        uint64_t id2 = localID++;
        id           = id1 + id2;
    }
    return local_comm.bcast( id, 0 );
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
    if ( !mesh )
        return;
    if ( !getProperties().registerMesh )
        AMP_ERROR( "registerMesh is not supported for " + getProperties().type );
    AMP_INSIST( level >= 0 && level <= 3, "Invalid value for level" );
    auto multimesh = std::dynamic_pointer_cast<AMP::Mesh::MultiMesh>( mesh );
    if ( !multimesh ) {
        uint64_t id = mesh->meshID().getData();
        // We are dealing with a single mesh
        baseMeshData data;
        data.id        = id;
        data.mesh      = mesh;
        data.rank      = mesh->getComm().getRank() + 1;
        data.ownerRank = d_comm.getRank();
        data.meshName  = "rank_" + std::to_string( data.rank );
        data.path      = path + mesh->getName() + "_/";
        if ( d_baseMeshes.find( id ) == d_baseMeshes.end() )
            d_baseMeshes.insert( std::make_pair( id, data ) );
        // Create and register a multimesh for the current mesh
        if ( level > 0 ) {
            multiMeshData data2;
            data2.id       = id;
            data2.mesh     = mesh;
            data2.name     = path + mesh->getName();
            data.ownerRank = mesh->getComm().bcast( d_comm.getRank(), 0 );
            d_multiMeshes.insert( std::make_pair( data2.id, data2 ) );
        }
        // Create and register a multimesh for the rank
        if ( level == 3 ) {
            // Create a unique id for each rank
            uint64_t tmp_id = id;
            uint64_t root2  = d_comm.getRank() + 1;
            tmp_id          = ( root2 << 48 ) + tmp_id;
            multiMeshData data2;
            data2.id       = tmp_id;
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
        for ( auto &submesh : submeshes )
            registerMesh( submesh, level2, new_path );
        if ( level > 0 ) {
            multiMeshData data;
            data.id        = mesh->meshID().getData();
            data.mesh      = mesh;
            data.name      = path + mesh->getName();
            data.ownerRank = mesh->getComm().bcast( d_comm.getRank(), 0 );
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
    auto id = getID( vec->getComm(), d_comm );
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
    VectorData data( vec, name_in );
    data.type    = type;
    data.numDOFs = DOFsPerPoint;
    auto ids     = getMeshIDs( mesh );
    for ( auto id : ids ) {
        const auto &it = d_baseMeshes.find( id.getData() );
        if ( it == d_baseMeshes.end() )
            continue;
        it->second.vectors.push_back( data );
    }
    // Register the vector with the appropriate multi-meshes
    auto it = d_multiMeshes.find( mesh->meshID().getData() );
    AMP_ASSERT( it != d_multiMeshes.end() );
    it->second.varName.push_back( data.name );
    // Add the vector to the list of vectors so we can perform makeConsistent
    d_vectorsMesh.push_back( vec );
    // Add the variable name to the list of variables
    d_varNames.insert( data.name );
}
#endif


/************************************************************
 * Register a matrix                                         *
 ************************************************************/
void Writer::registerMatrix( std::shared_ptr<AMP::LinearAlgebra::Matrix> mat,
                             const std::string &name )
{
#ifdef USE_AMP_MATRICES
    auto id = getID( mat->getLeftVector()->getComm(), d_comm );
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
        // We are dealining with a multimesh
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
    size_t N_bytes = sizeof( uint64_t ); // Store the mesh id
    N_bytes += sizeof( int );            // Store the processor rank
    N_bytes += sizeof( int );            // Store the owner rank
    N_bytes += meshName.size() + 1;      // Store the mesh name
    N_bytes += path.size() + 1;          // Store the mesh path
    N_bytes += file.size() + 1;          // Store the mesh file
    N_bytes += sizeof( int );            // Store the number of variables
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
    packData<uint64_t>( ptr, pos, id );
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
    data.id        = unpackData<uint64_t>( ptr, pos );
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
 * Functions for multiMeshData                           *
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
    size_t N_bytes = sizeof( uint64_t ); // Store the mesh id
    N_bytes += sizeof( int );            // Store the owner rank
    N_bytes += name.size() + 1;          // Store the mesh name
    N_bytes += sizeof( int );            // Store the number of sub meshes
    for ( const auto &mesh : meshes )
        N_bytes += mesh.size(); // Store the sub meshes
    N_bytes += sizeof( int );   // Store the number of variables
    for ( const auto &name : varName )
        N_bytes += name.size() + 1; // Store the variable name
    return N_bytes;
}
void Writer::multiMeshData::pack( char *ptr ) const
{
    size_t pos = 0;
    packData<uint64_t>( ptr, pos, id );
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
Writer::multiMeshData Writer::multiMeshData::unpack( const char *ptr )
{
    size_t pos = 0;
    multiMeshData data;
    data.id        = unpackData<uint64_t>( ptr, pos );
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
    return data;
}


} // namespace AMP::Utilities

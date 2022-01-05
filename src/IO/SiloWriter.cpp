#include "AMP/IO/SiloWriter.h"
#include "AMP/utils/Utilities.h"

#include "AMP/matrices/Matrix.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MultiMesh.h"
#include "AMP/vectors/Vector.h"

#include "ProfilerApp.h"

#include <chrono>


namespace AMP::IO {


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
SiloIO::SiloIO() : AMP::IO::Writer()
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
    properties.registerVectorWithMesh = true;
    return properties;
}


#if defined( USE_EXT_SILO )

// Some internal functions
static void createSiloDirectory( DBfile *FileHandle, const std::string &path );


/************************************************************
 * Function to read a silo file                              *
 ************************************************************/
void SiloIO::readFile( const std::string & ) { AMP_ERROR( "readFile is not implemented yet" ); }


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
    d_dim = -1;
    for ( auto tmp : d_baseMeshes ) {
        int dim = tmp.second.mesh->getDim();
        if ( d_dim == -1 )
            d_dim = tmp.second.mesh->getDim();
        AMP_INSIST( d_dim == dim, "All meshes must have the same number of physical dimensions" );
    }
    int dim = d_comm.maxReduce( d_dim );
    if ( d_dim == -1 )
        d_dim = dim;
    AMP_INSIST( d_dim == dim, "All meshes must have the same number of physical dimensions" );
    d_comm.barrier();
    PROFILE_STOP( "sync dim", 1 );
    // Synchronize the vectors
    syncVectors();
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
 * Function to write a mesh                                  *
 ************************************************************/
void SiloIO::writeMesh( DBfile *FileHandle, const baseMeshData &data, int cycle, double time )
{
    NULL_USE( cycle );
    NULL_USE( time );
    PROFILE_START( "writeMesh", 1 );
    auto mesh = data.mesh;
    // Get the zone (element) lists
    const auto type     = mesh->getGeomType();
    const auto elements = mesh->getIterator( type, 0 );
    AMP::Array<double> x[3];
    AMP::Array<int> nodelist;
    std::vector<AMP::Mesh::MeshElementID> nodelist_ids;
    getNodeElemList( mesh, elements, x, nodelist, nodelist_ids );
    // Get the shape type
    int shapecnt  = elements.size();
    int shapesize = nodelist.size( 0 );
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
    auto element_iterator = elements.begin();
    int num_elems         = elements.size();
    DBPutZonelist2( FileHandle,
                    zoneName.c_str(),
                    num_elems,
                    d_dim,
                    nodelist.data(),
                    nodelist.length(),
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
    double *coord[3] = { x[0].data(), x[1].data(), x[2].data() };
    DBPutUcdmesh( FileHandle,
                  meshName.c_str(),
                  d_dim,
                  nullptr,
                  coord,
                  x[0].length(),
                  nodelist.length(),
                  zoneName.c_str(),
                  nullptr,
                  DB_DOUBLE,
                  nullptr );
    PROFILE_STOP( "writeMesh - mesh", 2 );
    // Write the variables
    PROFILE_START( "writeMesh - variables", 2 );
    float ftime        = time;
    DBoptlist *optlist = DBMakeOptlist( 10 );
    DBAddOption( optlist, DBOPT_CYCLE, &cycle );
    DBAddOption( optlist, DBOPT_TIME, &ftime );
    DBAddOption( optlist, DBOPT_DTIME, &time );
    // DBAddOption(optlist, DBOPT_UNITS, (void *)units);
    for ( const auto &vector : data.vectors ) {
        int varSize   = vector.numDOFs;
        auto varType  = vector.type;
        auto DOFs     = vector.vec->getDOFManager();
        int nvar      = 0;
        int centering = 0;
        auto var      = new double *[varSize];
        for ( int j = 0; j < varSize; ++j )
            var[j] = nullptr;
        const char *varnames[] = { "1", "2", "3" };
        if ( varType > mesh->getGeomType() ) {
            // We have a mixed mesh type and there will be no data of the given type for this mesh
            continue;
        } else if ( varType == AMP::Mesh::GeomType::Vertex ) {
            // We are saving node-centered data
            centering = DB_NODECENT;
            nvar      = (int) nodelist_ids.size();
            for ( int j = 0; j < varSize; ++j )
                var[j] = new double[nvar];
            std::vector<size_t> dofs( varSize );
            std::vector<double> vals( varSize );
            for ( int j = 0; j < nvar; ++j ) {
                DOFs->getDOFs( nodelist_ids[j], dofs );
                AMP_ASSERT( (int) dofs.size() == varSize );
                vector.vec->getValuesByGlobalID( varSize, &dofs[0], &vals[0] );
                for ( int k = 0; k < varSize; ++k )
                    var[k][j] = vals[k];
            }
        } else if ( varType == mesh->getGeomType() ) {
            // We are saving cell-centered data
            centering = DB_ZONECENT;
            nvar      = (int) num_elems;
            for ( int j = 0; j < varSize; ++j )
                var[j] = new double[nvar];
            std::vector<size_t> dofs( varSize );
            std::vector<double> vals( varSize );
            auto it = element_iterator.begin();
            for ( int j = 0; j < nvar; ++j, ++it ) {
                DOFs->getDOFs( it->globalID(), dofs );
                vector.vec->getValuesByGlobalID( varSize, &dofs[0], &vals[0] );
                for ( int k = 0; k < varSize; ++k )
                    var[k][j] = vals[k];
            }
        } else {
            // We are storing edge or face data
            AMP_ERROR( "The silo writer currently only supports GeomType::Vertex and Cell data" );
        }
        std::string varNameRank = vector.name + "P" + std::to_string( data.rank );
        if ( varSize == 1 || varSize == d_dim || varSize == d_dim * d_dim ) {
            // We are writing a scalar, vector, or tensor variable
            DBPutUcdvar( FileHandle,
                         varNameRank.c_str(),
                         meshName.c_str(),
                         varSize,
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
            for ( int j = 0; j < varSize; ++j ) {
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
        for ( int j = 0; j < varSize; ++j ) {
            if ( var[j] )
                delete[] var[j];
        }
        delete[] var;
    }
    DBFreeOptlist( optlist );
    PROFILE_STOP( "writeMesh - variables", 2 );
    // Change the directory back to root
    DBSetDir( FileHandle, "/" );
    PROFILE_STOP( "writeMesh", 1 );
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
    // Gather the results
    // Note: we only need to guarantee that rank 0 has all the data
    const auto [multiMeshes, baseMeshes] = syncMultiMeshData( 0 );
    // Write the multimeshes and multivariables
    std::string base_path;
    if ( find_slash( filename ) != std::string::npos )
        base_path = filename.substr( 0, find_slash( filename ) + 1 );
    if ( d_comm.getRank() == 0 ) {
        DBfile *FileHandle = DBOpen( filename.c_str(), DB_HDF5, DB_APPEND );
        // Create the subdirectories
        PROFILE_START( "create directories", 2 );
        std::set<std::string> subdirs;
        for ( const auto &data : multiMeshes ) {
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
        for ( const auto &data : multiMeshes ) {
            size_t N = data.meshes.size();
            std::vector<std::string> meshNames( N );
            for ( size_t i = 0; i < N; ++i ) {
                auto it = baseMeshes.find( data.meshes[i] );
                AMP_ASSERT( it != baseMeshes.end() );
                const auto &base = it->second;
                auto file        = getFile( base.file, base_path );
                meshNames[i]     = file + ":" + base.path + "/" + base.meshName;
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
        for ( const auto &data : multiMeshes ) {
            size_t N = data.meshes.size();
            // std::cout << data.name << std::endl;
            for ( const auto &varName : data.varName ) {
                std::vector<std::string> varNames( N );
                auto varnames = new char *[N];
                auto vartypes = new int[N];
                for ( size_t i = 0; i < N; ++i ) {
                    auto it = baseMeshes.find( data.meshes[i] );
                    AMP_ASSERT( it != baseMeshes.end() );
                    auto base    = it->second;
                    auto rankStr = std::to_string( base.rank );
                    auto file    = getFile( base.file, base_path );
                    varNames[i]  = file + ":" + base.path + "/" + varName + "P" + rankStr;
                    strrep( varNames[i], "//", "/" );
                    varnames[i] = (char *) varNames[i].c_str();
                    vartypes[i] = DB_UCDVAR;
                }
                auto it = baseMeshes.find( data.meshes[0] );
                AMP_ASSERT( it != baseMeshes.end() );
                auto base   = it->second;
                int varSize = 0;
                for ( size_t i = 0; i < base.vectors.size(); ++i ) {
                    if ( base.vectors[i].name == varName ) {
                        varSize = base.vectors[i].numDOFs;
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
 * Some utility functions                                    *
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
void SiloIO::readFile( const std::string & ) {}
void SiloIO::writeFile( const std::string &, size_t, double ) {}
#endif


} // namespace AMP::IO

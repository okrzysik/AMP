#include "AMP/utils/HDF5writer.h"
#include "AMP/utils/HDF5_IO.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Xdmf.h"

#ifdef USE_AMP_MESH
#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshIterator.h"
#include "AMP/ampmesh/MultiMesh.h"
#include "AMP/ampmesh/structured/BoxMesh.h"
#else
namespace AMP::Mesh {
enum class GeomType : uint8_t { Vertex = 0, Edge = 1, Face = 2, Volume = 3, null = 0xFF };
}
#endif
#ifdef USE_AMP_VECTORS
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/data/ArrayVectorData.h"
#endif
#ifdef USE_AMP_MATRICES
#include "AMP/matrices/Matrix.h"
#endif

#include "ProfilerApp.h"

#include <chrono>


namespace AMP::Utilities {


/************************************************************
 * Helper functions                                          *
 ************************************************************/
#ifdef USE_AMP_VECTORS
AMP::Array<double> getArrayData( std::shared_ptr<const AMP::LinearAlgebra::Vector> vec )
{
    auto multivec = std::dynamic_pointer_cast<const AMP::LinearAlgebra::MultiVector>( vec );
    if ( multivec ) {
        AMP_ASSERT( multivec->getNumberOfSubvectors() == 1 );
        vec = multivec->getVector( 0 );
    }
    auto data = vec->getVectorData();
    auto arrayData =
        std::dynamic_pointer_cast<const AMP::LinearAlgebra::ArrayVectorData<double>>( data );
    if ( arrayData )
        return arrayData->getArray();
    size_t N = data->getLocalSize();
    AMP::Array<double> data2( N );
    auto it = data->constBegin();
    for ( size_t i = 0; i < N; i++, ++it )
        data2( i ) = *it;
    return data2;
}
#endif


/************************************************************
 * Constructor/Destructor                                    *
 ************************************************************/
HDF5writer::HDF5writer() : AMP::Utilities::Writer() {}
HDF5writer::~HDF5writer() = default;


/************************************************************
 * Some basic functions                                      *
 ************************************************************/
Writer::WriterProperties HDF5writer::getProperties() const
{
    WriterProperties properties;
    properties.type           = "HDF5";
    properties.extension      = "hdf5";
    properties.registerVector = true;
    properties.registerMesh   = true;
    // properties.registerVectorWithMesh = true;
    // properties.registerMatrix = true;
    return properties;
}


/************************************************************
 * Function to read a silo file                              *
 ************************************************************/
void HDF5writer::readFile( const std::string & ) { AMP_ERROR( "readFile is not implimented yet" ); }


/************************************************************
 * Function to write the hdf5 file                           *
 * Note: it appears that only one processor may write to a   *
 * file at a time, and that once a processor closes the file *
 * it cannot reopen it (or at least doing this on the        *
 * processor that created the file creates problems).        *
 ************************************************************/
void HDF5writer::writeFile( const std::string &fname_in, size_t cycle, double time )
{
#ifdef USE_EXT_HDF5
    Xdmf xmf;
    AMP_ASSERT( d_comm.getSize() == 1 );
    // Create the file
    auto filename = fname_in + "_" + std::to_string( cycle ) + ".hdf5";
    auto fid      = openHDF5( filename, "w", Compression::GZIP );
    writeHDF5( fid, "time", time );
    // Add the mesh
#ifdef USE_AMP_MESH
    std::map<GlobalID, Xdmf::MeshData> baseMeshData;
    auto gid = createGroup( fid, "meshes" );
    for ( const auto &[id, mesh] : d_baseMeshes )
        baseMeshData[id] = writeMesh( gid, mesh, Utilities::filename( filename ) + ":/meshes" );
    for ( const auto &[id, mesh] : d_multiMeshes ) {
        NULL_USE( id );
        std::vector<Xdmf::MeshData> data;
        for ( const auto &id2 : mesh.meshes ) {
            auto it = baseMeshData.find( id2 );
            if ( it != baseMeshData.end() )
                data.push_back( it->second );
        }
        xmf.addMultiMesh( mesh.name, data );
    }
    closeGroup( gid );
#endif
    // Add the vectors
#ifdef USE_AMP_VECTORS
    for ( const auto &[id, data] : d_vectors ) {
        NULL_USE( id );
        auto data2 = getArrayData( data.vec );
        writeHDF5( fid, data.name, data2 );
    }
#endif
    // Add the matricies
#ifdef USE_AMP_MATRICES
    for ( size_t i = 0; i < d_matrices.size(); i++ ) {
        AMP_ERROR( "Not finished" );
    }
#endif
    // Close the file
    closeHDF5( fid );
    // Write the Xdmf file
    xmf.gather( d_comm );
    if ( !xmf.empty() ) {
        auto fname = fname_in + "_" + std::to_string( cycle ) + ".xmf";
        xmf.write( fname );
    }
#else
    // No HDF5
    NULL_USE( fname_in );
    NULL_USE( cycle );
    NULL_USE( time );
#endif
}


/************************************************************
 * Function to write a base mesh                             *
 ************************************************************/
Xdmf::MeshData HDF5writer::writeDefaultMesh( hid_t fid,
                                             std::shared_ptr<const AMP::Mesh::Mesh> mesh,
                                             const std::string &name,
                                             const std::string &path )
{
    // Treat the mesh as an unstructured mesh
    const int ndim      = mesh->getDim();
    const auto type     = mesh->getGeomType();
    const auto elements = mesh->getIterator( type, 0 );
    AMP::Array<double> x[3];
    AMP::Array<int> nodelist;
    std::vector<AMP::Mesh::MeshElementID> nodelist_ids;
    getNodeElemList( mesh, elements, x, nodelist, nodelist_ids );
    auto shapetype = AMP::Xdmf::TopologyType::Null;
    int shapesize  = nodelist.size( 0 );
    if ( shapesize == 8 && type == AMP::Mesh::GeomType::Volume )
        shapetype = AMP::Xdmf::TopologyType::Hexahedron;
    else if ( shapesize == 4 && type == AMP::Mesh::GeomType::Volume )
        shapetype = AMP::Xdmf::TopologyType::Tetrahedron;
    else if ( shapesize == 4 && type == AMP::Mesh::GeomType::Face )
        shapetype = AMP::Xdmf::TopologyType::Quadrilateral;
    else if ( shapesize == 3 && type == AMP::Mesh::GeomType::Face )
        shapetype = AMP::Xdmf::TopologyType::Triangle;
    else if ( shapesize == 2 && type == AMP::Mesh::GeomType::Edge )
        shapetype = AMP::Xdmf::TopologyType::Polyline;
    else
        AMP_ERROR( "Unknown element type" );
    const char *x_names[3] = { "x", "y", "z" };
    std::string x_path[3];
    for ( int d = 0; d < ndim; d++ ) {
        writeHDF5( fid, x_names[d], x[d] );
        x_path[d] = path + "/" + x_names[d];
    }
    writeHDF5( fid, "type", static_cast<int>( type ) );
    writeHDF5( fid, "elements", nodelist );
    return AMP::Xdmf::createUnstructuredMesh( name,
                                              ndim,
                                              shapetype,
                                              elements.size(),
                                              path + "/elements",
                                              x[0].length(),
                                              x_path[0],
                                              x_path[1],
                                              x_path[2] );
}
Xdmf::MeshData HDF5writer::writeBoxMesh( hid_t fid,
                                         std::shared_ptr<const AMP::Mesh::BoxMesh> mesh,
                                         const std::string &name,
                                         const std::string &path )
{
    Xdmf::MeshData XdmfData;
    AMP::ArraySize size( mesh->size() );
    if ( size.ndim() != mesh->getDim() )
        return writeDefaultMesh( fid, mesh, name, path ); // We have issues with surface meshes
    auto size2 = size;
    for ( int d = 0; d < size.ndim(); d++ )
        size2.resize( d, size[d] + 1 );
    AMP::Array<double> x( size2 ), y( size2 ), z( size2 );
    for ( size_t k = 0; k < size2[2]; k++ ) {
        for ( size_t j = 0; j < size2[1]; j++ ) {
            for ( size_t i = 0; i < size2[0]; i++ ) {
                double pos[3] = { 0 };
                auto index =
                    AMP::Mesh::BoxMesh::MeshElementIndex( AMP::Mesh::GeomType::Vertex, 0, i, j, k );
                mesh->coord( index, pos );
                x( i, j, k ) = pos[0];
                y( i, j, k ) = pos[1];
                z( i, j, k ) = pos[2];
            }
        }
    }
    if ( mesh->getDim() == 1 ) {
        writeHDF5( fid, "x", x );
        XdmfData = AMP::Xdmf::createCurvilinearMesh( name, size, path + "/x" );
    } else if ( mesh->getDim() == 2 ) {
        writeHDF5( fid, "x", x );
        writeHDF5( fid, "y", y );
        XdmfData = AMP::Xdmf::createCurvilinearMesh( name, size, path + "/x", path + "/y" );
    } else if ( mesh->getDim() == 3 ) {
        writeHDF5( fid, "x", x );
        writeHDF5( fid, "y", y );
        writeHDF5( fid, "z", z );
        XdmfData =
            AMP::Xdmf::createCurvilinearMesh( name, size, path + "/x", path + "/y", path + "/z" );
    } else {
        AMP_ERROR( "Not finished" );
    }
    return XdmfData;
}
Xdmf::MeshData HDF5writer::writeMesh( hid_t fid, const baseMeshData &mesh, const std::string &path )
{
    Xdmf::MeshData XdmfData;
    if ( !mesh.mesh )
        return XdmfData;
    // Write the mesh
    auto name  = mesh.mesh->getName();
    auto gid   = createGroup( fid, name );
    auto path2 = path + "/" + name;
    auto name2 = name + "_" + mesh.meshName;
    writeHDF5( gid, "ndim", (int) mesh.mesh->getDim() );
    writeHDF5( gid, "meshClass", mesh.mesh->meshClass() );
    if ( std::dynamic_pointer_cast<AMP::Mesh::BoxMesh>( mesh.mesh ) ) {
        // We are dealing with a box mesh
        auto mesh2 = std::dynamic_pointer_cast<AMP::Mesh::BoxMesh>( mesh.mesh );
        XdmfData   = writeBoxMesh( gid, mesh2, name2, path2 );
    } else {
        XdmfData = writeDefaultMesh( gid, mesh.mesh, name2, path2 );
    }

    // Write the geometry (if it exists)

    closeGroup( gid );
    return XdmfData;
}


} // namespace AMP::Utilities

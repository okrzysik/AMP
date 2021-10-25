#include "AMP/utils/HDF5writer.h"
#include "AMP/utils/HDF5_IO.h"
#include "AMP/utils/Utilities.h"

#ifdef USE_AMP_MESH
#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MultiMesh.h"
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
auto getArrayData( std::shared_ptr<const AMP::LinearAlgebra::Vector> vec )
{
    auto multivec = std::dynamic_pointer_cast<const AMP::LinearAlgebra::MultiVector>( vec );
    if ( multivec ) {
        AMP_ASSERT( multivec->getNumberOfSubvectors() == 1 );
        vec = multivec->getVector( 0 );
    }
    return std::dynamic_pointer_cast<const AMP::LinearAlgebra::ArrayVectorData<double>>(
        vec->getVectorData() );
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
    properties.type      = "HDF5";
    properties.extension = "hdf5";
#ifdef USE_EXT_HDF5
#ifdef USE_AMP_VECTORS
    properties.registerVector = true;
#endif
#ifdef USE_AMP_MESH
    properties.registerMesh = false;
#endif
#endif
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
    AMP_ASSERT( d_comm.getSize() == 1 );
    // Create the file
    auto filename = fname_in + "_" + std::to_string( cycle ) + ".hdf5";
    auto fid      = openHDF5( filename, "w", Compression::GZIP );
    writeHDF5( fid, "time", time );
    // Add the mesh
#ifdef USE_AMP_MESH
    if ( !d_baseMeshes.empty() || !d_multiMeshes.empty() )
        AMP_ERROR( "Not finished" );
#endif
        // Add the vectors
#ifdef USE_AMP_VECTORS
    for ( auto vec : d_vec ) {
        auto arrayData = getArrayData( vec.vec );
        if ( arrayData ) {
            writeHDF5( fid, vec.name, arrayData->getArray() );
        } else {
            AMP_ERROR( "Not finished" );
        }
    }
#endif
    // Add the matricies
#ifdef USE_AMP_MESH
    for ( size_t i = 0; i < d_mat.size(); i++ ) {
        AMP_ERROR( "Not finished" );
    }
#endif
    // Close the file
    closeHDF5( fid );
#else
    // No HDF5
    NULL_USE( fname_in );
    NULL_USE( cycle );
    NULL_USE( time );
#endif
}
void HDF5writer::writeMesh( MeshData data )
{
    auto multimesh = std::dynamic_pointer_cast<AMP::Mesh::MultiMesh>( data );
    if ( multimesh ) {
        for ( auto mesh : multimesh->getMeshes() )
            writeMesh( mesh );
        AMP_ERROR( "Not finished" );
    } else {
        AMP_ERROR( "Not finished" );
    }
}


} // namespace AMP::Utilities

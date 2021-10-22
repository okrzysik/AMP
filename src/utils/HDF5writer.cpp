#include "AMP/utils/HDF5writer.h"
#include "AMP/utils/HDF5_IO.h"
#include "AMP/utils/Utilities.h"

#ifdef USE_AMP_MESH
#include "AMP/ampmesh/Mesh.h"
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
#endif
    return properties;
}


/************************************************************
 * Function to read a silo file                              *
 ************************************************************/
void HDF5writer::readFile( const std::string & ) { AMP_ERROR( "readFile is not implimented yet" ); }


/************************************************************
 * Function to write the hdf5 file                           *
 * Note: it appears that only one prcoessor may write to a   *
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
    for ( size_t i = 0; i < d_mesh.size(); i++ ) {
        AMP_ERROR( "Not finished" );
    }
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


/************************************************************
 * Function to register a mesh with silo                     *
 ************************************************************/
void HDF5writer::registerMesh( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                               int level,
                               const std::string &path )
{
    /*    d_mesh
          \param level How many sub meshes do we want?
                       0: Only register the local base meshes (advanced users only)
                       1: Register current mesh only (default)
                       2: Register all meshes (do not seperate for the ranks)
                       3: Register all mesh pieces including the individual ranks
          \param path  The directory path for the mesh.  Default is an empty string.
        void registerMesh( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                           int level               = 1,
                           const std::string &path = std::string() ) override;*/
}


/************************************************************
 * Function to register a vector with silo                   *
 ************************************************************/
void HDF5writer::registerVector( std::shared_ptr<AMP::LinearAlgebra::Vector>,
                                 std::shared_ptr<AMP::Mesh::Mesh>,
                                 AMP::Mesh::GeomType,
                                 const std::string & )
{
    AMP_ERROR( "Meshes are not supported yet" );
}
void HDF5writer::registerVector( std::shared_ptr<AMP::LinearAlgebra::Vector> vec,
                                 const std::string &name )
{
    VectorData data;
    data.name = name;
    data.vec  = vec;
    data.type = AMP::Mesh::GeomType::null;
    data.mesh = nullptr;
    d_vec.push_back( data );
}
void HDF5writer::registerMatrix( std::shared_ptr<AMP::LinearAlgebra::Matrix> mat,
                                 const std::string &name )
{
    MatrixData data;
    data.name = name;
    data.mat  = mat;
    d_mat.push_back( data );
}


} // namespace AMP::Utilities

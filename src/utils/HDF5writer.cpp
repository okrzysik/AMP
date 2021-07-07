#include "AMP/utils/HDF5writer.h"
#include "AMP/ampmesh/MultiMesh.h"
#include "AMP/utils/HDF5_IO.h"
#include "AMP/utils/Utilities.h"

#ifdef USE_AMP_MESH
#include "AMP/ampmesh/Mesh.h"
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
    properties.type                   = "HDF5";
    properties.extension              = "hdf";
    properties.registerMesh           = false;
    properties.registerVector         = true;
    properties.registerVectorWithMesh = false;
    properties.registerMatrix         = false;
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
    NULL_USE( fname_in );
    NULL_USE( cycle );
    NULL_USE( time );
#ifdef USE_EXT_HDF5
    AMP_ASSERT( d_comm.getSize() == 1 );
    // Create the file
    auto filename = fname_in + "_" + std::to_string( cycle ) + ".hdf5";
    auto fid      = openHDF5( filename, "w", Compression::GZIP );
    writeHDF5( fid, "time", time );
    // Add the vectors
#ifdef USE_AMP_VECTORS
    for ( size_t i = 0; i < d_vec.size(); i++ ) {
        auto arrayData = getArrayData( d_vec[i].vec );
        if ( arrayData ) {
            writeHDF5( fid, d_vec[i].name, arrayData->getArray() );
        } else {
            AMP_ERROR( "Not finished" );
        }
    }
#endif
    // Close the file
    closeHDF5( fid );
#endif
}


/************************************************************
 * Function to register a mesh with silo                     *
 ************************************************************/
void HDF5writer::registerMesh( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                               int level,
                               const std::string &path )
{
    NULL_USE( mesh );
    NULL_USE( level );
    NULL_USE( path );
    AMP_ERROR( "Meshes are not supported yet" );
}


/************************************************************
 * Function to register a vector with silo                   *
 ************************************************************/
void HDF5writer::registerVector( std::shared_ptr<AMP::LinearAlgebra::Vector>,
                                 AMP::Mesh::Mesh::shared_ptr,
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

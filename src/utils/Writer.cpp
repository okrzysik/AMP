#include "AMP/utils/Writer.h"
#include "AMP/utils/Utilities.h"

#include "AMP/utils/AsciiWriter.h"
#include "AMP/utils/HDF5writer.h"
#include "AMP/utils/NullWriter.h"
#include "AMP/utils/SiloWriter.h"


namespace AMP::Utilities {


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
 * Define default registers                                  *
 ************************************************************/
void Writer::registerMesh( std::shared_ptr<AMP::Mesh::Mesh>, int, const std::string & ) {}
void Writer::registerVector( std::shared_ptr<AMP::LinearAlgebra::Vector>,
                             std::shared_ptr<AMP::Mesh::Mesh>,
                             AMP::Mesh::GeomType,
                             const std::string & )
{
}
void Writer::registerVector( std::shared_ptr<AMP::LinearAlgebra::Vector>, const std::string & ) {}
void Writer::registerMatrix( std::shared_ptr<AMP::LinearAlgebra::Matrix>, const std::string & ) {}


} // namespace AMP::Utilities

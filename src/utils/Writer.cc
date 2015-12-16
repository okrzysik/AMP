#include "utils/Writer.h"
#include "utils/Utilities.h"

#include "utils/AsciiWriter.h"
#include "utils/NullWriter.h"
#ifdef USE_AMP_MESH
#include "ampmesh/SiloIO.h"
#endif


namespace AMP {
namespace Utilities {


/************************************************************
* Builder                                                   *
************************************************************/
AMP::shared_ptr<AMP::Utilities::Writer> Writer::buildWriter( std::string type )
{
    AMP::shared_ptr<AMP::Utilities::Writer> writer;
    if ( type == "None" || type == "none" || type == "NONE" ) {
        writer.reset( new AMP::Utilities::NullWriter() );
    } else if ( type == "Silo" || type == "silo" || type == "SILO" ) {
#if defined( USE_AMP_MESH ) && defined( USE_EXT_SILO )
        writer.reset( new AMP::Mesh::SiloIO() );
#else
        writer.reset( new AMP::Utilities::NullWriter() );
#endif
    } else if ( type == "Ascii" || type == "ascii" || type == "ASCII" ) {
        writer.reset( new AMP::Utilities::AsciiWriter() );
    } else {
        AMP_ERROR( "Unknown writer" );
    }
    return writer;
}
AMP::shared_ptr<AMP::Utilities::Writer> Writer::buildWriter( AMP::shared_ptr<AMP::Database> db )
{
    std::string type                               = db->getString( "Name" );
    AMP::shared_ptr<AMP::Utilities::Writer> writer = Writer::buildWriter( type );
    if ( db->keyExists( "Decomposition" ) )
        writer->setDecomposition( db->getInteger( "Decomposition" ) );
    return writer;
}


/************************************************************
* Constructor/Destructor                                    *
************************************************************/
Writer::Writer()
{
    d_comm          = AMP_MPI( AMP_COMM_WORLD );
    d_decomposition = 2;
}
Writer::~Writer() {}


/************************************************************
* Some basic functions                                      *
************************************************************/
void Writer::setDecomposition( int d )
{
    AMP_INSIST( d == 1 || d == 2, "decomposition must be 1 or 2" );
    d_decomposition = d;
}
}
}

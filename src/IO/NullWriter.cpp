#include "AMP/IO/NullWriter.h"


namespace AMP::IO {


/************************************************************
 * Some basic functions                                      *
 ************************************************************/
Writer::WriterProperties NullWriter::getProperties() const
{
    WriterProperties properties;
    properties.type                   = "Null";
    properties.extension              = "";
    properties.registerMesh           = false;
    properties.registerVector         = false;
    properties.registerVectorWithMesh = false;
    properties.registerMatrix         = false;
    return properties;
}


} // namespace AMP::IO

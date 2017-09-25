#include "vectors/VectorEngine.h"


namespace AMP {
namespace LinearAlgebra {


/********************************************************
 * VectorEngineParameters constructors             *
 ********************************************************/
VectorEngineParameters::VectorEngineParameters( size_t local_size, size_t global_size, AMP_MPI c )
    : d_begin( 0 ), d_end( 0 ), d_global( 0 ), d_comm( c.getCommunicator() )
{
    d_global = global_size;
    d_comm.sumScan( &local_size, &d_end, 1 );
    d_begin = d_end - local_size;
}
VectorEngineParameters::~VectorEngineParameters() {}


} // namespace LinearAlgebra
} // namespace AMP

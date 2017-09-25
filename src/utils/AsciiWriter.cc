#include "utils/AsciiWriter.h"
#include "ProfilerApp.h"

#ifdef USE_AMP_VECTORS
#include "vectors/SimpleVector.h"
#endif


namespace AMP {
namespace Utilities {


/************************************************************
 * Helper function to get a unique id for each vector        *
 ************************************************************/
static unsigned int localID = 0;
AsciiWriter::global_id AsciiWriter::getID( AMP_MPI local_comm, AMP_MPI global_comm )
{
    AsciiWriter::global_id id( 0, 0 );
    if ( local_comm.getRank() == 0 ) {
        unsigned int id2 = localID++;
        id               = AsciiWriter::global_id( global_comm.getRank(), id2 );
    }
    return local_comm.bcast<AsciiWriter::global_id>( id, 0 );
}


/************************************************************
 * Constructor/Destructor                                    *
 ************************************************************/
AsciiWriter::AsciiWriter() : AMP::Utilities::Writer() {}
AsciiWriter::~AsciiWriter() = default;


/************************************************************
 * Some basic functions                                      *
 ************************************************************/
std::string AsciiWriter::getExtension() { return "ascii"; }


/************************************************************
 * Function to read a silo file                              *
 ************************************************************/
void AsciiWriter::readFile( const std::string & )
{
    AMP_ERROR( "readFile is not implemented yet" );
}


/************************************************************
 * Function to write an ascii file                           *
 ************************************************************/
template<class TYPE>
std::set<AsciiWriter::global_id> getKeys( const std::map<AsciiWriter::global_id, TYPE> &local_map,
                                          AMP_MPI comm )
{
    std::set<AsciiWriter::global_id> ids;
    typename std::map<AsciiWriter::global_id, TYPE>::const_iterator it;
    for ( it = local_map.begin(); it != local_map.end(); ++it )
        ids.insert( it->first );
    comm.setGather( ids );
    return ids;
}
void AsciiWriter::writeFile( const std::string &fname_in, size_t iteration_count, double time )
{
    NULL_USE( time );
    PROFILE_START( "writeFile" );
    // Open the file for writing
    FILE *fid = nullptr;
    if ( d_comm.getRank() == 0 ) {
        std::stringstream tmp;
        tmp << fname_in << "_" << iteration_count << "." << getExtension();
        std::string fname = tmp.str();
        fid               = fopen( fname.c_str(), "w" );
        AMP_ASSERT( fid != nullptr );
    }
// Get the ids for the vectors and save the data
#ifdef USE_AMP_VECTORS
    std::set<global_id> vec_ids = getKeys( d_vectors, d_comm );
    for ( const auto &vec_id : vec_ids ) {
        // Send the data to rank 0
        d_comm.barrier();
        AMP::LinearAlgebra::Vector::shared_ptr src_vec;
        if ( d_vectors.find( vec_id ) != d_vectors.end() )
            src_vec = d_vectors[vec_id];
        AMP::LinearAlgebra::Vector::const_shared_ptr dst_vec =
            sendVecToRoot( src_vec, vec_id.first, d_comm );
        // Write the data
        if ( d_comm.getRank() == 0 ) {
            fprintf( fid,
                     "Vector: \"%s\" %i\n",
                     dst_vec->getVariable()->getName().c_str(),
                     static_cast<int>( dst_vec->getGlobalSize() ) );
            for ( size_t i = 0; i < dst_vec->getGlobalSize(); i++ )
                fprintf( fid, "   %0.14e\n", dst_vec->getValueByGlobalID( i ) );
            fprintf( fid, "\n\n" );
        }
    }
#endif
// Get the ids for the matricies and save the data
#ifdef USE_AMP_MATRICES
    std::set<global_id> mat_ids = getKeys( d_matrices, d_comm );
    for ( const auto &mat_id : mat_ids ) {
        // Send the header data to rank 0
        d_comm.barrier();
        AMP::LinearAlgebra::Matrix::shared_ptr mat;
        if ( d_matrices.find( mat_id ) != d_matrices.end() )
            mat = d_matrices[mat_id];
        std::string name;
        size_t size[2] = { 0, 0 };
        if ( mat != nullptr ) {
            name = mat->getLeftVector()->getVariable()->getName() + " - " +
                   mat->getRightVector()->getVariable()->getName();
            size[0] = mat->getLeftVector()->getGlobalSize();
            size[1] = mat->getRightVector()->getGlobalSize();
        }
        name    = d_comm.bcast( name, mat_id.first );
        size[0] = d_comm.bcast( size[0], mat_id.first );
        size[1] = d_comm.bcast( size[1], mat_id.first );
        // Write the data
        if ( d_comm.getRank() == 0 ) {
            fprintf( fid,
                     "Matrix: \"%s\" %i %i\n",
                     name.c_str(),
                     static_cast<int>( size[0] ),
                     static_cast<int>( size[1] ) );
        }
        std::vector<size_t> col;
        std::vector<double> data;
        for ( int row = 0; row < static_cast<int>( size[0] ); row++ ) {
            // Get and print the current row
            sendRowToRoot( mat, d_comm, row, col, data );
            if ( d_comm.getRank() == 0 ) {
                for ( size_t i = 0; i < col.size(); i++ )
                    fprintf( fid, "   %4i %4u  %0.14e\n", row, (unsigned int) col[i], data[i] );
            }
        }
        if ( d_comm.getRank() == 0 ) {
            fprintf( fid, "\n\n" );
        }
    }
#endif
    // Close the file
    if ( fid != nullptr )
        fclose( fid );
    PROFILE_STOP( "writeFile" );
}


/************************************************************
 * Function to register a mesh                               *
 ************************************************************/
#ifdef USE_AMP_MESH
void AsciiWriter::registerMesh( AMP::Mesh::Mesh::shared_ptr, int, std::string )
{
    AMP_ERROR( "registerMesh is not implimented yet" );
}
#endif


/************************************************************
 * Function to register a vector                             *
 ************************************************************/
#if defined( USE_AMP_MESH ) && defined( USE_AMP_VECTORS )
void AsciiWriter::registerVector( AMP::LinearAlgebra::Vector::shared_ptr,
                                  AMP::Mesh::Mesh::shared_ptr,
                                  AMP::Mesh::GeomType,
                                  const std::string & )
{
    AMP_ERROR( "Mesh support is not implimented yet" );
}
#endif
#ifdef USE_AMP_VECTORS
void AsciiWriter::registerVector( AMP::LinearAlgebra::Vector::shared_ptr vec )
{
    global_id id = getID( vec->getComm(), d_comm );
    d_vectors.insert( std::pair<global_id, AMP::LinearAlgebra::Vector::shared_ptr>( id, vec ) );
}
#endif
#ifdef USE_AMP_MATRICES
void AsciiWriter::registerMatrix( AMP::LinearAlgebra::Matrix::shared_ptr mat )
{
    global_id id = getID( mat->getLeftVector()->getComm(), d_comm );
    d_matrices.insert( std::pair<global_id, AMP::LinearAlgebra::Matrix::shared_ptr>( id, mat ) );
}
#endif


/************************************************************
 * Function to copy a vector to rank 0                       *
 ************************************************************/
#ifdef USE_AMP_VECTORS
AMP::LinearAlgebra::Vector::const_shared_ptr AsciiWriter::sendVecToRoot(
    AMP::LinearAlgebra::Vector::const_shared_ptr src_vec, int vec_root, AMP_MPI comm )
{
    int rank = comm.getRank();
    // Boradcast the local vector size and name to all processors for simplicity
    std::string name;
    if ( rank == vec_root )
        name = src_vec->getVariable()->getName();
    name              = comm.bcast( name, vec_root );
    size_t local_size = 0;
    if ( src_vec != nullptr )
        local_size = src_vec->getLocalSize();
    std::vector<size_t> size( comm.getSize(), 0 );
    comm.allGather( local_size, &size[0] );
    size_t global_size = 0;
    for ( int i = 0; i < comm.getSize(); i++ )
        global_size += size[i];
    // If we are not rank 0 and do not have a copy of the vector we are done
    if ( src_vec == nullptr && rank != 0 )
        return AMP::LinearAlgebra::Vector::const_shared_ptr();
    // Send the local data to rank 0
    std::vector<MPI_Request> requests;
    std::vector<double> local_data( local_size, 0 );
    if ( local_size > 0 ) {
        for ( size_t i = 0; i < local_size; i++ )
            local_data[i] = src_vec->getValueByLocalID( i );
        requests.push_back( comm.Isend( &local_data[0], local_size, 0, 123 ) );
    }
    // Rank 0 needs to create the vector and recv all data
    AMP::LinearAlgebra::Vector::shared_ptr dst_vec;
    if ( rank == 0 ) {
        AMP::LinearAlgebra::Variable::shared_ptr var( new AMP::LinearAlgebra::Variable( name ) );
        dst_vec = AMP::LinearAlgebra::SimpleVector<double>::create(
            global_size, var, AMP_MPI( AMP_COMM_SELF ) );
        AMP_ASSERT( dst_vec->numberOfDataBlocks() == 1 );
        auto *ptr = dst_vec->getRawDataBlock<double>( 0 );
        size_t i  = 0;
        for ( int j = 0; j < comm.getSize(); j++ ) {
            if ( size[j] > 0 ) {
                requests.push_back( comm.Irecv( &ptr[i], size[j], j, 123 ) );
                i += size[j];
            }
        }
    }
    if ( !requests.empty() )
        comm.waitAll( requests.size(), &requests[0] );
    return dst_vec;
}
#endif


/************************************************************
 * Function to copy a row to rank 0                          *
 ************************************************************/
#ifdef USE_AMP_MATRICES
void AsciiWriter::sendRowToRoot( AMP::LinearAlgebra::Matrix::const_shared_ptr mat,
                                 AMP_MPI comm,
                                 int row,
                                 std::vector<size_t> &cols,
                                 std::vector<double> &data )
{
    int rank = comm.getRank();
    cols.clear();
    data.clear();
    // Determine who "owns" the row
    int own_rank = 0;
    if ( mat != nullptr ) {
        AMP::Discretization::DOFManager::shared_ptr DOF = mat->getLeftDOFManager();
        if ( row >= (int) DOF->beginDOF() && row < (int) DOF->endDOF() )
            own_rank = rank;
    }
    own_rank = comm.maxReduce( own_rank );
    // Send the data
    std::vector<MPI_Request> requests;
    if ( own_rank == rank ) {
        mat->getRowByGlobalID( row, cols, data );
        if ( rank == 0 )
            return;
        size_t size = cols.size();
        requests.push_back( comm.Isend<size_t>( &size, 1, 0, 124 ) );
        if ( size > 0 ) {
            requests.push_back( comm.Isend<size_t>( &cols[0], size, 0, 125 ) );
            requests.push_back( comm.Isend<double>( &data[0], size, 0, 126 ) );
        }
    }
    // Recv the data
    if ( rank == 0 ) {
        size_t size = 0;
        int length  = 1;
        comm.recv<size_t>( &size, length, own_rank, false, 124 );
        cols.resize( size, 0 );
        data.resize( size, 0 );
        if ( size > 0 ) {
            requests.push_back( comm.Irecv<size_t>( &cols[0], size, own_rank, 125 ) );
            requests.push_back( comm.Irecv<double>( &data[0], size, own_rank, 126 ) );
        }
    }
    if ( !requests.empty() )
        comm.waitAll( requests.size(), &requests[0] );
}
#endif


} // namespace Utilities
} // namespace AMP

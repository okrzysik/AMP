#include "AMP/matrices/ManagedMatrixParameters.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Utilities.h"

namespace AMP {
namespace LinearAlgebra {


ManagedMatrixParameters::ManagedMatrixParameters( AMP::Discretization::DOFManager::shared_ptr left,
                                                  AMP::Discretization::DOFManager::shared_ptr right,
                                                  const AMP_MPI &comm )
    : MatrixParameters( left, right, comm )
{
    d_vEntriesPerRow.resize( getLocalNumberOfRows() );
}


void ManagedMatrixParameters::addColumns( int a, int *b )
{
    for ( int i = 0; i != a; i++ )
        d_sColumns.insert( b[i] );
}
void ManagedMatrixParameters::addColumns( const std::set<size_t> &col )
{
    for ( auto b : col )
        d_sColumns.insert( b );
}


const int *ManagedMatrixParameters::entryList() const { return &*d_vEntriesPerRow.begin(); }

int *ManagedMatrixParameters::entryList() { return &( d_vEntriesPerRow[0] ); }

void ManagedMatrixParameters::setEntriesInRow( int row, int entries )
{
    d_vEntriesPerRow[row] = entries;
}

int &ManagedMatrixParameters::entriesInRow( int i ) { return d_vEntriesPerRow[i]; }

int ManagedMatrixParameters::entriesInRow( int i ) const { return d_vEntriesPerRow[i]; }

int ManagedMatrixParameters::maxEntitiesInRow() const
{
    return *std::max_element( d_vEntriesPerRow.begin(), d_vEntriesPerRow.end() );
}

bool ManagedMatrixParameters::isSquare()
{
    AMP_ASSERT( d_DOFManagerLeft.get() != nullptr );
    AMP_ASSERT( d_DOFManagerRight.get() != nullptr );
    return d_DOFManagerLeft->numLocalDOF() == d_DOFManagerRight->numLocalDOF() &&
           d_DOFManagerLeft->numGlobalDOF() == d_DOFManagerRight->numGlobalDOF();
}


} // namespace LinearAlgebra
} // namespace AMP

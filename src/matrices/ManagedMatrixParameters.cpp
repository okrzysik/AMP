#include "AMP/matrices/ManagedMatrixParameters.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Utilities.h"

namespace AMP::LinearAlgebra {


ManagedMatrixParameters::ManagedMatrixParameters(
    std::shared_ptr<AMP::Discretization::DOFManager> left,
    std::shared_ptr<AMP::Discretization::DOFManager> right,
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
    AMP_ASSERT( d_DOFManagerLeft );
    AMP_ASSERT( d_DOFManagerRight );
    return d_DOFManagerLeft->numLocalDOF() == d_DOFManagerRight->numLocalDOF() &&
           d_DOFManagerLeft->numGlobalDOF() == d_DOFManagerRight->numGlobalDOF();
}


} // namespace AMP::LinearAlgebra

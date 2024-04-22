#include "AMP/matrices/MatrixParameters.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Vector.h"


namespace AMP::LinearAlgebra {


MatrixParameters::MatrixParameters( std::shared_ptr<AMP::Discretization::DOFManager> left,
                                    std::shared_ptr<AMP::Discretization::DOFManager> right,
                                    const AMP_MPI &comm )
    : MatrixParametersBase( comm )
{
    AMP_ASSERT( left );
    AMP_ASSERT( right );
    d_DOFManagerLeft  = left;
    d_DOFManagerRight = right;
    d_vEntriesPerRow.resize( getLocalNumberOfRows() );
}


size_t MatrixParameters::getLocalNumberOfRows() const { return d_DOFManagerLeft->numLocalDOF(); }


size_t MatrixParameters::getLocalNumberOfColumns() const
{
    return d_DOFManagerRight->numLocalDOF();
}

size_t MatrixParameters::getGlobalNumberOfRows() const { return d_DOFManagerLeft->numGlobalDOF(); }

size_t MatrixParameters::getGlobalNumberOfColumns() const
{
    return d_DOFManagerRight->numGlobalDOF();
}

std::shared_ptr<AMP::Discretization::DOFManager> MatrixParameters::getLeftDOFManager()
{
    return d_DOFManagerLeft;
}

std::shared_ptr<AMP::Discretization::DOFManager> MatrixParameters::getRightDOFManager()
{
    return d_DOFManagerRight;
}


void MatrixParameters::addColumns( const std::vector<size_t> &col )
{
    d_vColumns.insert( d_vColumns.end(), col.begin(), col.end() );
}

void MatrixParameters::findUniqueColumns()
{
#warning This code is non-optimal and only present to test matvec performance ideas
  
  std::copy(d_vColumns.begin(), d_vColumns.end(), std::back_inserter(d_vUniqueColumns));
  AMP::Utilities::unique( d_vUniqueColumns );
  // Create map from vector of columns into vector of unique columns
  for ( size_t nc=0; nc < d_vColumns.size(); ++nc ) {
    const auto col = d_vColumns[nc];
    auto lower = std::lower_bound(d_vUniqueColumns.begin(), d_vUniqueColumns.end(), col);
    const auto idx = std::distance(d_vUniqueColumns.begin(), lower);
    d_vUniqueColumnsMap.push_back(idx);
  }
  AMP::pout << "matrix has " << d_vUniqueColumns.size() << " unique cols out of " << d_vColumns.size() << std::endl;
}

const int *MatrixParameters::entryList() const { return &*d_vEntriesPerRow.begin(); }

int *MatrixParameters::entryList() { return &( d_vEntriesPerRow[0] ); }

void MatrixParameters::setEntriesInRow( int row, int entries ) { d_vEntriesPerRow[row] = entries; }

int &MatrixParameters::entriesInRow( int i ) { return d_vEntriesPerRow[i]; }

int MatrixParameters::entriesInRow( int i ) const { return d_vEntriesPerRow[i]; }

} // namespace AMP::LinearAlgebra

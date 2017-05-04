#include "VectorData.h"


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
* Set/Get individual values                                     *
****************************************************************/
void VectorData::setValuesByGlobalID( int numVals, size_t *ndx, const double *vals )
{
    AMP_ASSERT( *d_UpdateState != UpdateState::ADDING );
    *d_UpdateState = UpdateState::SETTING;
    for ( int i = 0; i < numVals; i++ ) {
        if ( ( ndx[i] < getLocalStartID() ) ||
             ( ndx[i] >= ( getLocalStartID() + getLocalMaxID() ) ) ) {
            ( *d_Ghosts )[d_CommList->getLocalGhostID( ndx[i] )] = vals[i];
        } else {
            setLocalValuesByGlobalID( 1, ndx + i, vals + i );
        }
    }
}
void VectorData::setGhostValuesByGlobalID( int numVals, size_t *ndx, const double *vals )
{
    AMP_ASSERT( *d_UpdateState != UpdateState::ADDING );
    *d_UpdateState = UpdateState::SETTING;
    for ( int i = 0; i < numVals; i++ ) {
        if ( ( ndx[i] < getLocalStartID() ) ||
             ( ndx[i] >= ( getLocalStartID() + getLocalMaxID() ) ) ) {
            ( *d_Ghosts )[d_CommList->getLocalGhostID( ndx[i] )] = vals[i];
        } else {
            AMP_ERROR( "Non ghost index" );
        }
    }
}
void VectorData::addValuesByGlobalID( int numVals, size_t *ndx, const double *vals )
{
    AMP_ASSERT( *d_UpdateState != UpdateState::SETTING );
    *d_UpdateState = UpdateState::ADDING;
    for ( int i = 0; i < numVals; i++ ) {
        if ( ( ndx[i] < getLocalStartID() ) ||
             ( ndx[i] >= ( getLocalStartID() + getLocalMaxID() ) ) ) {
            ( *d_AddBuffer )[d_CommList->getLocalGhostID( ndx[i] )] += vals[i];
        } else {
            addLocalValuesByGlobalID( 1, ndx + i, vals + i );
        }
    }
}
void VectorData::getValuesByLocalID( int num, size_t *ndx, double *vals ) const
{
    for ( int i = 0; i != num; i++ ) {
        size_t block_number = 0;
        size_t offset       = ndx[i];
        while ( offset >= sizeOfDataBlock( block_number ) ) {
            offset -= sizeOfDataBlock( block_number );
            block_number++;
            if ( block_number >= numberOfDataBlocks() ) {
                AMP_ERROR( "Bad local id!" );
            }
        }
        vals[i] = getRawDataBlock<double>( block_number )[offset];
    }
}
void VectorData::getValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
{
    for ( int i = 0; i < numVals; i++ ) {
        if ( ( ndx[i] < getLocalStartID() ) ||
             ( ndx[i] >= ( getLocalStartID() + getLocalMaxID() ) ) ) {
            getGhostValuesByGlobalID( 1, ndx + i, vals + i );
        } else {
            getLocalValuesByGlobalID( 1, ndx + i, vals + i );
        }
    }
}
void VectorData::getGhostValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
{
    for ( int i = 0; i < numVals; i++ ) {
        if ( ( ndx[i] < getLocalStartID() ) ||
             ( ndx[i] >= ( getLocalStartID() + getLocalMaxID() ) ) ) {
            vals[i] = ( *d_Ghosts )[d_CommList->getLocalGhostID( ndx[i] )] +
                      ( *d_AddBuffer )[d_CommList->getLocalGhostID( ndx[i] )];
        } else {
            AMP_ERROR( "Tried to get a non-ghost ghost value" );
        }
    }
}
void VectorData::getGhostAddValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
{
    for ( int i = 0; i < numVals; i++ ) {
        if ( ( ndx[i] < getLocalStartID() ) ||
             ( ndx[i] >= ( getLocalStartID() + getLocalMaxID() ) ) ) {
            vals[i] = ( *d_AddBuffer )[d_CommList->getLocalGhostID( ndx[i] )];
        } else {
            AMP_ERROR( "Tried to get a non-ghost ghost value" );
        }
    }
}


} // LinearAlgebra namespace
} // AMP namespace


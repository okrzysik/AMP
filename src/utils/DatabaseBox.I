//
// File:	$URL:
// file:///usr/casc/samrai/repository/AMP/tags/v-2-4-4/source/toolbox/database/DatabaseBox.I $
// Package:	AMP toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	A box structure representing a portion of the AMR index space
//

#include "Utilities.h"

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace AMP {


inline DatabaseBox::DatabaseBox()
{
    d_data.d_dimension = 0;
    d_data.d_lo[0] = d_data.d_hi[0] = 0;
    d_data.d_lo[1] = d_data.d_hi[1] = 0;
    d_data.d_lo[2] = d_data.d_hi[2] = 0;
}

inline DatabaseBox::DatabaseBox( const DatabaseBox &box )
{
    d_data.d_dimension = box.d_data.d_dimension;
    d_data.d_lo[0]     = box.d_data.d_lo[0];
    d_data.d_lo[1]     = box.d_data.d_lo[1];
    d_data.d_lo[2]     = box.d_data.d_lo[2];
    d_data.d_hi[0]     = box.d_data.d_hi[0];
    d_data.d_hi[1]     = box.d_data.d_hi[1];
    d_data.d_hi[2]     = box.d_data.d_hi[2];
}

inline DatabaseBox &DatabaseBox::operator=( const DatabaseBox &box )
{
    d_data.d_dimension = box.d_data.d_dimension;
    d_data.d_lo[0]     = box.d_data.d_lo[0];
    d_data.d_lo[1]     = box.d_data.d_lo[1];
    d_data.d_lo[2]     = box.d_data.d_lo[2];
    d_data.d_hi[0]     = box.d_data.d_hi[0];
    d_data.d_hi[1]     = box.d_data.d_hi[1];
    d_data.d_hi[2]     = box.d_data.d_hi[2];
    return ( *this );
}

inline DatabaseBox::~DatabaseBox() {}

inline int DatabaseBox::getDimension() const { return ( d_data.d_dimension ); }

inline void DatabaseBox::setDimension( const int dimension )
{
    AMP_ASSERT( ( dimension >= 0 ) && ( dimension <= DatabaseBox_MAX_DIM ) );
    d_data.d_dimension = dimension;
}

inline int &DatabaseBox::lower( const int i )
{
    AMP_ASSERT( ( i >= 0 ) && ( i < d_data.d_dimension ) );
    return ( d_data.d_lo[i] );
}

inline int &DatabaseBox::upper( const int i )
{
    AMP_ASSERT( ( i >= 0 ) && ( i < d_data.d_dimension ) );
    return ( d_data.d_hi[i] );
}

inline int DatabaseBox::lower( const int i ) const
{
    AMP_ASSERT( ( i >= 0 ) && ( i < d_data.d_dimension ) );
    return ( d_data.d_lo[i] );
}

inline int DatabaseBox::upper( const int i ) const
{
    AMP_ASSERT( ( i >= 0 ) && ( i < d_data.d_dimension ) );
    return ( d_data.d_hi[i] );
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}

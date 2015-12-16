#include "DataChangeFirer.h"
#include "utils/Utilities.h"

namespace AMP {
namespace LinearAlgebra {

inline DataChangeListener::DataChangeListener() {}

inline DataChangeListener::~DataChangeListener()
{
    for ( auto &elem : *this )
        ( elem )->deregisterListener( this );
}

inline void DataChangeListener::registerWithFirer( DataChangeFirer *firer )
{
    AMP_ASSERT( firer );
    AMP_ASSERT( std::find( begin(), end(), firer ) == end() );

    push_back( firer );
}

inline void DataChangeListener::deregisterFromFirer( DataChangeFirer *firer )
{
    AMP_ASSERT( firer );
    AMP_ASSERT( std::find( begin(), end(), firer ) != end() );

    erase( std::find( begin(), end(), firer ) );
}
}
}

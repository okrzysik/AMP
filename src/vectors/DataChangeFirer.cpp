#include "DataChangeFirer.h"
#include "AMP/utils/Utilities.h"
#include "DataChangeListener.h"

#include <algorithm>


namespace AMP {
namespace LinearAlgebra {

DataChangeFirer::DataChangeFirer() {}

DataChangeFirer::~DataChangeFirer()
{
    for ( auto &x : d_listeners )
        x->deregisterFromFirer( this );
}

void DataChangeFirer::registerListener( std::shared_ptr<DataChangeListener> listener )
{
    AMP_ASSERT( listener );
    AMP_ASSERT( std::find( d_listeners.begin(), d_listeners.end(), listener.get() ) ==
                d_listeners.end() );
    listener->registerWithFirer( this );
    d_listeners.push_back( listener.get() );
}

void DataChangeFirer::deregisterListener( DataChangeListener *listener )
{
    AMP_ASSERT( listener );
    d_listeners.erase( std::find( d_listeners.begin(), d_listeners.end(), listener ) );
}


void DataChangeFirer::fireDataChange()
{
    for ( auto &x : d_listeners )
        x->recieveDataChanged();
}
} // namespace LinearAlgebra
} // namespace AMP

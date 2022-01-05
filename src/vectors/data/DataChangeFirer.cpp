#include "AMP/vectors/data/DataChangeFirer.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/data/DataChangeListener.h"

#include <algorithm>


namespace AMP::LinearAlgebra {

DataChangeFirer::DataChangeFirer() = default;

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

void DataChangeFirer::registerListener( DataChangeListener *listener )
{
    AMP_ASSERT( listener );
    AMP_ASSERT( std::find( d_listeners.begin(), d_listeners.end(), listener ) ==
                d_listeners.end() );
    listener->registerWithFirer( this );
    d_listeners.push_back( listener );
}

void DataChangeFirer::deregisterListener( DataChangeListener *listener )
{
    AMP_ASSERT( listener );
    d_listeners.erase( std::find( d_listeners.begin(), d_listeners.end(), listener ) );
}


void DataChangeFirer::fireDataChange()
{
    for ( auto &x : d_listeners )
        x->receiveDataChanged();
}
} // namespace AMP::LinearAlgebra

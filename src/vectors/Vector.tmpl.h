#ifndef included_AMP_Vector_tmpl
#define included_AMP_Vector_tmpl

namespace AMP {
namespace LinearAlgebra {

template <typename RETURN_TYPE>
RETURN_TYPE *Vector::getRawDataBlock( size_t i )
{
    return static_cast<RETURN_TYPE *>( this->getRawDataBlockAsVoid( i ) );
}

template <typename RETURN_TYPE>
const RETURN_TYPE *Vector::getRawDataBlock( size_t i ) const
{
    return static_cast<const RETURN_TYPE *>( this->getRawDataBlockAsVoid( i ) );
}

template <typename VIEW_TYPE>
Vector::shared_ptr Vector::getView() const
{
    for ( size_t i = 0; i != d_Views->size(); i++ ) {
        if ( ( *d_Views )[i].lock() ) {
            if ( ( *d_Views )[i].lock()->isA<VIEW_TYPE>() ) {
                return Vector::shared_ptr( ( *d_Views )[i] );
            }
        }
    }
    return Vector::shared_ptr();
}

template <typename VIEW_TYPE>
bool Vector::hasView() const
{
    for ( size_t i = 0; i != d_Views->size(); i++ ) {
        if ( ( *d_Views )[i].lock() ) {
            if ( ( *d_Views )[i].lock()->isA<VIEW_TYPE>() ) {
                return true;
            }
        }
    }
    return false;
}


} // namespace LinearAlgebra
} // namespace AMP

#endif

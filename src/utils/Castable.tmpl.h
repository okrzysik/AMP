

namespace AMP {


template <typename T>
T &Castable::castTo ()
{
    T *ret_val = dynamic_cast<T *> (this);
    if ( ret_val == 0 ) 
        AMP_ERROR("Invalid cast");
    return *ret_val;
}

template <typename T>
const T &Castable::castTo () const
{
    const T * ret_val = dynamic_cast<const T *> (this);
    if ( ret_val == 0 ) 
        AMP_ERROR("Invalid cast");
    return *ret_val ;
}


template <typename T>
bool Castable::isA () const
{
    T *ret_val = dynamic_cast<T *>(const_cast<Castable *>(this));
    return ret_val != 0;
}


} // AMP


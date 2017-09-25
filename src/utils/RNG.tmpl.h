namespace AMP {
template<typename T>
RandomVariable<T>::RandomVariable( type low, type high, RNG::shared_ptr r )
    : d_Low( low ), d_Size( high - low ), d_RNG( r )
{
    AMP_ASSERT( high > low );
}


template<typename T>
RandomVariable<T>::operator type()
{
    type retVal;
    d_RNG->fillBuffer( static_cast<void *>( &retVal ), sizeof( T ) );
    return ( abs( retVal ) % d_Size ) + d_Low;
}
} // namespace AMP

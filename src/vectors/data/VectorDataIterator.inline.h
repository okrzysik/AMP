#ifndef included_AMP_VectorIterators_inline
#define included_AMP_VectorIterators_inline


namespace AMP::LinearAlgebra {


/****************************************************************
 * Boolean Operators                                             *
 ****************************************************************/
template<typename TYPE>
inline bool VectorDataIterator<TYPE>::operator==( const VectorDataIterator<TYPE> &rhs ) const
{
    return d_hashcode == rhs.d_hashcode && d_pos == rhs.d_pos;
}
template<typename TYPE>
inline bool VectorDataIterator<TYPE>::operator!=( const VectorDataIterator<TYPE> &rhs ) const
{
    return d_hashcode != rhs.d_hashcode || d_pos != rhs.d_pos;
}
template<typename TYPE>
inline bool VectorDataIterator<TYPE>::operator<( const VectorDataIterator<TYPE> &rhs ) const
{
    return d_hashcode == rhs.d_hashcode ? ( d_pos < rhs.d_pos ) : ( d_hashcode < rhs.d_hashcode );
}
template<typename TYPE>
inline bool VectorDataIterator<TYPE>::operator>( const VectorDataIterator<TYPE> &rhs ) const
{
    return d_hashcode == rhs.d_hashcode ? ( d_pos > rhs.d_pos ) : ( d_hashcode > rhs.d_hashcode );
}
template<typename TYPE>
inline bool VectorDataIterator<TYPE>::operator<=( const VectorDataIterator<TYPE> &rhs ) const
{
    return d_hashcode == rhs.d_hashcode ? ( d_pos <= rhs.d_pos ) : ( d_hashcode <= rhs.d_hashcode );
}
template<typename TYPE>
inline bool VectorDataIterator<TYPE>::operator>=( const VectorDataIterator<TYPE> &rhs ) const
{
    return d_hashcode == rhs.d_hashcode ? ( d_pos >= rhs.d_pos ) : ( d_hashcode >= rhs.d_hashcode );
}


} // namespace AMP::LinearAlgebra


#endif

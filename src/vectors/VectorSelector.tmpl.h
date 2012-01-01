namespace AMP {
namespace LinearAlgebra {


template <typename T>
bool   VS_ByVariableType<T>::isSelected ( Vector::const_shared_ptr v ) const
{
    return v->getVariable()->isA<T>();
}


}
}


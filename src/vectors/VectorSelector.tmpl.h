namespace AMP {
namespace LinearAlgebra {


template <typename T>
bool  VS_ByVariableType<T>::isSelected ( Vector::const_shared_ptr v ) const
{
    T *ret_val = dynamic_cast<T *>(const_cast<Variable *>(v->getVariable().get()));
    return ret_val != 0;
}


}
}


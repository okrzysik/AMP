
namespace AMP::LinearAlgebra {


inline std::shared_ptr<Matrix> Matrix::transpose() const
{
    AMP_ERROR( "not implemented" );
    return std::shared_ptr<Matrix>();
}

inline std::ostream &operator<<( std::ostream &out, const std::shared_ptr<Matrix> p )
{
    return operator<<( out, *p );
}
} // namespace AMP::LinearAlgebra

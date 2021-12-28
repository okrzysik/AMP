
namespace AMP::LinearAlgebra {


inline std::shared_ptr<Matrix> Matrix::transpose() const
{
    AMP_ERROR( "not implemented" );
    return std::shared_ptr<Matrix>();
}

inline void Matrix::addValueByGlobalID( size_t row, size_t col, double value )
{
    addValuesByGlobalID( 1u, 1u, &row, &col, &value );
}

inline void Matrix::setValueByGlobalID( size_t row, size_t col, double value )
{
    setValuesByGlobalID( 1u, 1u, &row, &col, &value );
}

inline double Matrix::getValueByGlobalID( size_t row, size_t col ) const
{
    double rtn = 0.0;
    getValuesByGlobalID( 1u, 1u, &row, &col, &rtn );
    return rtn;
}

inline std::ostream &operator<<( std::ostream &out, const std::shared_ptr<Matrix> p )
{
    return operator<<( out, *p );
}
} // namespace AMP::LinearAlgebra

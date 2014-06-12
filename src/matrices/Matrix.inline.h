
namespace AMP {
namespace LinearAlgebra {


inline Matrix::shared_ptr  Matrix::transpose () const 
{ 
    AMP_ERROR( "not implemented" ); 
    return Matrix::shared_ptr();
}

inline Matrix::Matrix ( const Matrix& rhs ):
    d_comm(rhs.d_comm)
{
}

inline Matrix::Matrix () 
{
}

inline Matrix::Matrix ( MatrixParameters::shared_ptr params ):
    d_comm(params->getComm())
{
}

inline Matrix :: ~Matrix () 
{
}

inline void  Matrix::axpy ( double alpha , const Matrix::shared_ptr  &x ) 
{ 
    axpy( alpha , *x ); 
}

inline void Matrix::addValueByGlobalID ( int row , int col , double value )
{  
    addValuesByGlobalID( 1u , 1u , &row , &col , &value ); 
}

inline void Matrix::setValueByGlobalID ( int row , int col , double value )
{  
    setValuesByGlobalID( 1u , 1u , &row , &col , &value ); 
}

inline double Matrix::getValueByGlobalID ( int row , int col ) const
{  
    double rtn = 0.0;
    getValuesByGlobalID( 1u , 1u , &row , &col , &rtn ); 
    return rtn;
}

inline std::ostream &operator << ( std::ostream &out , const Matrix::shared_ptr p )
{
    return operator << ( out , *p );
}


}
}


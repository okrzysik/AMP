#ifndef NEWTON_SOLVER_H
#define NEWTON_SOLVER_H


#include <cmath>
#include <iostream>


enum class status_t { CONV_RTOL, CONV_ATOL, DIV_MAXIT, DIV_TOL };
std::ostream &operator<<( std::ostream &os, status_t const &status )
{
    if ( status == status_t::CONV_RTOL ) {
        os << "CONV_RTOL";
    } else if ( status == status_t::CONV_ATOL ) {
        os << "CONV_ATOL";
    } else if ( status == status_t::DIV_MAXIT ) {
        os << "DIV_MAXIT";
    } else if ( status == status_t::DIV_TOL ) {
        os << "DIV_TOL";
    } else {
        AMP_ASSERT( false );
    } // end if
    return os;
}
struct solve_status_t {
    solve_status_t( status_t s, size_t i, double r )
        : status( s ), n_iterations( i ), final_residual_norm( r )
    {
    }
    status_t status;
    size_t n_iterations;
    double final_residual_norm;
};

template<typename vector_t>
class newton_solver_t
{
public:
    newton_solver_t() : _atol( 1.0e-16 ), _rtol( 1.0e-12 ), _dtol( 1.0e3 ), _maxit( 30 ) {}
    void set( void ( * )( vector_t const &, vector_t &, void * ),
              void ( * )( vector_t const &, vector_t const &, vector_t &, void * ) );
    solve_status_t apply( vector_t &x, void *parameters );
    double compute_norm( vector_t const &x );
    void update_solution( vector_t const &dx, vector_t &x );

private:
    double _atol, _rtol, _dtol;
    size_t _maxit;
    void ( *compute_residual )( vector_t const &solution, vector_t &residual, void *parameters );
    void ( *compute_action_inverse_jacobian_on_minus_residual )( vector_t const &solution,
                                                                 vector_t const &residual,
                                                                 vector_t &update,
                                                                 void *parameters );
};

template<typename vector_t>
static void newton_solver_t<vector_t>::set(
    void ( *my_compute_residual )( vector_t const &, vector_t &, void * ),
    void ( *my_compute_action_inverse_jacobian_on_minus_residual )(
        vector_t const &, vector_t const &, vector_t &, void * ) )
{
    compute_residual = my_compute_residual;
    compute_action_inverse_jacobian_on_minus_residual =
        my_compute_action_inverse_jacobian_on_minus_residual;
}

template<typename vector_t>
solve_status_t newton_solver_t<vector_t>::apply( vector_t &solution, void *parameters )
{
    vector_t residual( solution );
    vector_t update( solution );
    compute_residual( solution, residual, parameters );
    double const initial_residual_norm = compute_norm( residual );
    double norm_residual;
    for ( size_t i = 1; i <= _maxit; ++i ) {
        compute_action_inverse_jacobian_on_minus_residual( solution, residual, update, parameters );
        update_solution( update, solution );
        compute_residual( solution, residual, parameters );
        norm_residual = compute_norm( residual );
        if ( norm_residual < _atol ) {
            return solve_status_t( status_t::CONV_ATOL, i, norm_residual );
        } else if ( norm_residual < _rtol * initial_residual_norm ) {
            return solve_status_t( status_t::CONV_RTOL, i, norm_residual );
        } else if ( norm_residual > _dtol * initial_residual_norm ) {
            return solve_status_t( status_t::DIV_TOL, i, norm_residual );
        } // end if
    } // end for i
    return solve_status_t( DIV_MAXIT, _maxit, norm_residual );
}

template<>
double newton_solver_t<double>::compute_norm( double const &x )
{
    return std::abs( x );
}

template<>
void newton_solver_t<double>::update_solution( double const &dx, double &x )
{
    x += dx;
}

#endif // NEWTON_SOLVER_H

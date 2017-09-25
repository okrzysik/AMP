#include "solvers/libmesh/Flow1DSolver.h"
#include "operators/subchannel/FlowFrapconJacobian.h"
#include "utils/Utilities.h"
#include "vectors/MultiVector.h"

namespace AMP {
namespace Solver {

Flow1DSolver::Flow1DSolver( AMP::shared_ptr<SolverStrategyParameters> parameters )
    : SolverStrategy( parameters )
{

    AMP_ASSERT( parameters.get() != nullptr );

    initialize( parameters );
}

Flow1DSolver::~Flow1DSolver() {}

void Flow1DSolver::setInitialGuess( AMP::shared_ptr<AMP::LinearAlgebra::Vector> ) {}

void Flow1DSolver::initialize( AMP::shared_ptr<SolverStrategyParameters> const parameters )
{
    getFromInput( parameters->d_db );

    if ( d_pOperator.get() != nullptr ) {
        AMP::shared_ptr<AMP::Operator::FlowFrapconJacobian> Operator =
            AMP::dynamic_pointer_cast<AMP::Operator::FlowFrapconJacobian>( d_pOperator );

        AMP_ASSERT( Operator.get() != nullptr );

        d_inpVariable = Operator->d_inpVariable;
        d_outVariable = Operator->d_outVariable;
        d_K           = Operator->d_K;
        d_De          = Operator->d_De;
        d_G           = Operator->d_G;
        d_Re          = Operator->d_Re;
        d_Pr          = Operator->d_Pr;
        d_Cp          = Operator->getCp();
        d_dCp         = Operator->getdCp();
        zPoints       = Operator->getZLocations();
        d_cladVec     = Operator->d_cladVec;
    }
}

void Flow1DSolver::reset( AMP::shared_ptr<SolverStrategyParameters> )
{

    if ( d_pOperator.get() != nullptr ) {
    }
}

void Flow1DSolver::resetOperator( const AMP::shared_ptr<AMP::Operator::OperatorParameters> params )
{
    if ( d_pOperator.get() != nullptr ) {
        d_pOperator->reset( params );
    }
}

void Flow1DSolver::solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                          AMP::shared_ptr<AMP::LinearAlgebra::Vector> u )
{

    AMP::LinearAlgebra::Vector::shared_ptr flowInputVec =
        u->subsetVectorForVariable( d_inpVariable );
    AMP::LinearAlgebra::Vector::const_shared_ptr flowRhsVec =
        f->constSubsetVectorForVariable( d_inpVariable );

    d_numpoints = zPoints.size();

    for ( int i = 1; i < d_numpoints; i++ ) {

        double cur_node, next_node;

        cur_node  = zPoints[i - 1];
        next_node = zPoints[i];

        // double Heff, he_z, T_b_i, T_b_im1, T_c_i;
        // double T_b =0, Tin, frhs;
        double Heff, he_z, T_b_i, T_b_im1;
        double frhs;

        Heff = ( 0.023 * d_K / d_De ) * pow( d_Re, 0.8 ) * pow( d_Pr, 0.4 );
        he_z = next_node - cur_node;

        // T_c_i   = d_cladVec->getValueByLocalID(i);
        T_b_im1 = flowInputVec->getValueByLocalID( i - 1 );
        T_b_i   = flowInputVec->getValueByLocalID( i );
        frhs    = flowRhsVec->getValueByLocalID( i );

        T_b_i =
            ( T_b_im1 + frhs ) * ( 1. / ( 1.0 + ( ( 4 * Heff * he_z ) / ( d_Cp * d_G * d_De ) ) ) );

        flowInputVec->setValueByLocalID( i, T_b_i );

    } // end for i
}


AMP::LinearAlgebra::Variable::shared_ptr Flow1DSolver::getInputVariable( int )
{
    return d_inpVariable;
}
} // namespace Solver
} // namespace AMP

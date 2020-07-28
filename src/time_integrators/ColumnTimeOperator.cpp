#include "ColumnTimeOperator.h"
#include "AMP/operators/ColumnOperatorParameters.h"
#include "AMP/utils/Database.h"
#include "LinearTimeOperator.h"
#include "TimeOperatorParameters.h"

namespace AMP {
namespace TimeIntegrator {

ColumnTimeOperator::ColumnTimeOperator(
    std::shared_ptr<AMP::Operator::OperatorParameters> in_params )
    : ColumnOperator( in_params )
{
    auto params    = std::dynamic_pointer_cast<TimeOperatorParameters>( in_params );
    auto column_db = params->d_db;
    d_Mesh         = params->d_Mesh;
    d_pRhsOperator = std::dynamic_pointer_cast<ColumnOperator>( params->d_pRhsOperator );
    AMP_INSIST( d_pRhsOperator.get() != nullptr,
                "Error: ColumnTimeOperator::ColumnTimeOperator() rhs "
                "operator must be a non-NULL column operator" );

    d_pMassOperator = std::dynamic_pointer_cast<ColumnOperator>( params->d_pMassOperator );
    AMP_INSIST( d_pRhsOperator.get() != nullptr,
                "Error: ColumnTimeOperator::ColumnTimeOperator() "
                "mass operator must be a non-NULL column operator" );

    d_bCreateLinearTimeOperators = column_db->getWithDefault( "CreateLinearTimeOperators", true );

    // for now we only allow one algebraic component
    d_iAlgebraicComponent = column_db->getWithDefault( "algebraicComponent", -1 );

    d_dCurrentDt = column_db->getWithDefault<double>( "CurrentTime", 0.0 );

    if ( d_bCreateLinearTimeOperators ) {
        const int numberOfOperators = d_pRhsOperator->getNumberOfOperators();

        for ( int i = 0; i < numberOfOperators; i++ ) {
            std::shared_ptr<AMP::Database> timeOperator_db(
                new AMP::Database( "TimeOperatorDatabase" ) );
            // we assume for now that either all operators in the column operator are linear or all
            // are nonlinear
            timeOperator_db->putScalar( "CurrentDt",
                                        column_db->getWithDefault<double>( "CurrentDt", 1.0e-08 ) );
            timeOperator_db->putScalar( "CurrentTime",
                                        column_db->getWithDefault<double>( "CurrentTime", 0.0 ) );
            timeOperator_db->putScalar( "name", "TimeOperator" );
            timeOperator_db->putScalar( "bLinearMassOperator",
                                        column_db->getWithDefault( "bLinearMassOperator", true ) );
            timeOperator_db->putScalar( "bLinearRhsOperator",
                                        column_db->getWithDefault( "bLinearRhsOperator", false ) );
            timeOperator_db->putScalar(
                "ScalingFactor", column_db->getWithDefault<double>( "ScalingFactor", 1.0e6 ) );
            std::shared_ptr<AMP::TimeIntegrator::TimeOperatorParameters> timeOperatorParameters(
                new AMP::TimeIntegrator::TimeOperatorParameters( timeOperator_db ) );
            timeOperatorParameters->d_pRhsOperator = d_pRhsOperator->getOperator( i );

            // if there are algebraic components set the mass operator to NULL
            if ( i != d_iAlgebraicComponent ) {
                timeOperatorParameters->d_pMassOperator = d_pMassOperator->getOperator( i );
            } else {
                timeOperator_db->putScalar( "bAlgebraicComponent", true );
            }

            std::shared_ptr<AMP::TimeIntegrator::LinearTimeOperator> op(
                new AMP::TimeIntegrator::LinearTimeOperator( timeOperatorParameters ) );

            d_Operators.push_back( op );
        }
    } else {
        AMP::pout << "Error: ColumnTimeOperator::ColumnTimeOperator() currently implemented for "
                     "column rhs and mass "
                     "operators where all components are LinearOperators"
                  << std::endl;
    }

    reset( in_params );

    d_pSourceTerm           = params->d_pSourceTerm;
    d_pPreviousTimeSolution = params->d_pPreviousTimeSolution;
}

ColumnTimeOperator::~ColumnTimeOperator() = default;

void ColumnTimeOperator::reset(
    const std::shared_ptr<AMP::Operator::OperatorParameters> &in_params )
{
    std::shared_ptr<TimeOperatorParameters> params =
        std::dynamic_pointer_cast<TimeOperatorParameters>( in_params );
    std::shared_ptr<AMP::Database> column_db = params->d_db;
    std::shared_ptr<AMP::Operator::ColumnOperatorParameters> pRhsParameters =
        std::dynamic_pointer_cast<AMP::Operator::ColumnOperatorParameters>(
            params->d_pRhsOperatorParameters );
    std::shared_ptr<AMP::Operator::ColumnOperatorParameters> pMassParameters =
        std::dynamic_pointer_cast<AMP::Operator::ColumnOperatorParameters>(
            params->d_pMassOperatorParameters );

    AMP_INSIST( params.get() != nullptr, "Error: NULL TimeOperatorParameters object" );

    getFromInput( params->d_db );

    const int numberOfOperators = d_pRhsOperator->getNumberOfOperators();

    for ( int i = 0; i < numberOfOperators; i++ ) {
        std::shared_ptr<AMP::Database> timeOperator_db(
            new AMP::Database( "TimeOperatorDatabase" ) );
        // we assume for now that either all operators in the column operator are linear or all are
        // nonlinear
        timeOperator_db->putScalar( "CurrentDt",
                                    column_db->getWithDefault<double>( "CurrentDt", 1.0e-08 ) );
        timeOperator_db->putScalar( "CurrentTime",
                                    column_db->getWithDefault<double>( "CurrentTime", 0.0 ) );
        timeOperator_db->putScalar( "name", "TimeOperator" );
        timeOperator_db->putScalar( "bLinearMassOperator",
                                    column_db->getWithDefault( "bLinearMassOperator", true ) );
        timeOperator_db->putScalar( "bLinearRhsOperator",
                                    column_db->getWithDefault( "bLinearRhsOperator", false ) );
        timeOperator_db->putScalar( "ScalingFactor",
                                    column_db->getWithDefault<double>( "ScalingFactor", 1.0e6 ) );
        std::shared_ptr<AMP::TimeIntegrator::TimeOperatorParameters> timeOperatorParameters(
            new AMP::TimeIntegrator::TimeOperatorParameters( timeOperator_db ) );
        if ( pRhsParameters.get() != nullptr ) {
            timeOperatorParameters->d_pRhsOperatorParameters =
                ( pRhsParameters->d_OperatorParameters )[i];
        }

        // if there are algebraic components set the mass operator to NULL
        if ( i != d_iAlgebraicComponent ) {
            if ( pMassParameters.get() != nullptr ) {
                timeOperatorParameters->d_pMassOperatorParameters =
                    ( pMassParameters->d_OperatorParameters )[i];
            }
        } else {
            timeOperator_db->putScalar( "bAlgebraicComponent", true );
        }

        d_Operators[i]->reset( timeOperatorParameters );
    }
}

void ColumnTimeOperator::getFromInput( std::shared_ptr<AMP::Database> db )
{
    AMP_INSIST( db->keyExists( "CurrentDt" ), "key CurrentDt missing in input" );

    d_dCurrentDt = db->getScalar<double>( "CurrentDt" );
}

void ColumnTimeOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr /* u */,
                                AMP::LinearAlgebra::Vector::shared_ptr /* f */ )
{
    AMP_ERROR( "Not Finished" );
}

void ColumnTimeOperator::append( std::shared_ptr<Operator> /* op */ )
{
    AMP::pout
        << "Error: ColumnTimeOperator::append(): this routine is disabled for ColumnTimeOperators"
        << std::endl;
}
} // namespace TimeIntegrator
} // namespace AMP

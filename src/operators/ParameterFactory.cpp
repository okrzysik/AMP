#include "AMP/operators/ParameterFactory.h"
#include "AMP/AMP_TPLs.h"

#ifdef AMP_USE_LIBMESH
    #include "AMP/operators/NeutronicsRhsParameters.h"
    #include "AMP/operators/boundary/DirichletMatrixCorrectionParameters.h"
    #include "AMP/operators/mechanics/MechanicsLinearFEOperatorParameters.h"
    #include "AMP/operators/mechanics/MechanicsNonlinearFEOperatorParameters.h"
#endif


#define resetParameters( NAME )                                      \
    do {                                                             \
        if ( name == #NAME )                                         \
            retParameters.reset( new NAME##Parameters( input_db ) ); \
    } while ( 0 )


namespace AMP::Operator {


std::shared_ptr<OperatorParameters>
ParameterFactory::createParameter( std::shared_ptr<AMP::Database> input_db,
                                   std::shared_ptr<AMP::Mesh::Mesh> mesh )
{
    AMP_ASSERT( input_db );
    auto name = input_db->getString( "name" );
    NULL_USE( name );
    std::shared_ptr<OperatorParameters> retParameters;
#ifdef AMP_USE_LIBMESH
    resetParameters( DirichletMatrixCorrection );
    resetParameters( MechanicsLinearFEOperator );
    resetParameters( MechanicsNonlinearFEOperator );
    resetParameters( NeutronicsRhs );
#endif
    if ( retParameters )
        retParameters->d_Mesh = mesh;
    return retParameters;
}


} // namespace AMP::Operator

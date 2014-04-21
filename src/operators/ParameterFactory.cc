#include "operators/ParameterFactory.h"

#ifdef USE_EXT_LIBMESH
    #include "operators/boundary/DirichletMatrixCorrectionParameters.h"
    #include "operators/mechanics/MechanicsLinearFEOperatorParameters.h"
    #include "operators/mechanics/MechanicsNonlinearFEOperatorParameters.h"
    #include "operators/NeutronicsRhsParameters.h"
#endif


#define resetParameters(NAME)                       \
    do {                                            \
        if ( name == #NAME )                        \
            retParameters.reset( new NAME ## Parameters(input_db) ); \
    } while(0)


namespace AMP {
namespace Operator {


boost::shared_ptr<OperatorParameters>
ParameterFactory::createParameter(boost::shared_ptr<AMP::Database>  input_db, AMP::Mesh::Mesh::shared_ptr  mesh)
{
    boost::shared_ptr<OperatorParameters> retParameters;
    std::string name;

    AMP_INSIST(input_db.get()!=NULL, "ParameterFactory::createParameter:: NULL Database object input");
    AMP_INSIST(input_db->keyExists("name"),  "ParameterFactory::createParameter:: key 'name' must be a part of database ");

    name = input_db->getString("name");
  
    #ifdef USE_EXT_LIBMESH
        resetParameters(DirichletMatrixCorrection);
        resetParameters(MechanicsLinearFEOperator);
        resetParameters(MechanicsNonlinearFEOperator);
        resetParameters(NeutronicsRhs);
    #endif

    AMP_ASSERT(retParameters!=NULL);

    retParameters->d_Mesh = mesh;
  
    return retParameters;
}
  

} // namespace Operator
} // namespace AMP


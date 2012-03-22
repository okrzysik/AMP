
#include "operators/Operator.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Operator {


int Operator :: d_iInstance_id=0;


Operator :: Operator(void)
{
    d_iObject_id                 = Operator::d_iInstance_id;

    d_iDebugPrintInfoLevel       = 0;   

    Operator :: d_iInstance_id++;
}
 

Operator :: Operator(const boost::shared_ptr<OperatorParameters> & params)
{
    AMP_INSIST( ((params.get()) != NULL), "NULL parameter" );

    d_iObject_id                 = Operator::d_iInstance_id;

    d_iDebugPrintInfoLevel       = 0;   

    Operator :: d_iInstance_id++;

    d_Mesh = params->d_Mesh;

    // try and keep the next call the last in the function
    // so as not to override any parameters set through it 
    // by accident
    getFromInput(params->d_db);
}
 

void Operator :: reset(const boost::shared_ptr<OperatorParameters>& params)
{
    AMP_INSIST( ((params.get()) != NULL), "NULL parameter" );

    // try and keep the next call the last in the function
    // so as not to override any parameters set through it 
    // by accident
    getFromInput(params->d_db);

}


void Operator :: getFromInput(const boost::shared_ptr<AMP::Database>& db)
{
    AMP_INSIST( ((db.get()) != NULL), "NULL database" );

    d_iDebugPrintInfoLevel = db->getIntegerWithDefault("print_info_level", 0);
}


AMP::LinearAlgebra::Vector::shared_ptr  Operator::subsetOutputVector(AMP::LinearAlgebra::Vector::shared_ptr vec)
{
    AMP::LinearAlgebra::Variable::shared_ptr varSubset = getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr varSubsetVec = vec->subsetVectorForVariable(varSubset);
    if(varSubsetVec == NULL) {
        return varSubsetVec;
    } else if( d_Mesh.get() !=NULL)  {
        AMP::LinearAlgebra::VS_Mesh meshSelector(varSubset->getName(), d_Mesh);
        return varSubsetVec->select(meshSelector, getOutputVariable()->getName());
    } else {
        return varSubsetVec; 
    } 
}


AMP::LinearAlgebra::Vector::shared_ptr  Operator::subsetInputVector(AMP::LinearAlgebra::Vector::shared_ptr vec)
{
    AMP::LinearAlgebra::Variable::shared_ptr varSubset = getInputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr varSubsetVec = vec->subsetVectorForVariable(varSubset);
    if(varSubsetVec == NULL) {
        return varSubsetVec;
    } else if( d_Mesh.get() !=NULL)  {
      AMP_ASSERT(d_Mesh != NULL);
        AMP::LinearAlgebra::VS_Mesh meshSelector(varSubset->getName(), d_Mesh);
        return varSubsetVec->select(meshSelector, getInputVariable()->getName());
    } else {
        return varSubsetVec; 
    }
}


}
}


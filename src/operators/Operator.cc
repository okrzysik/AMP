#include "operators/Operator.h"
#include "utils/Utilities.h"
#include "ProfilerApp.h"


namespace AMP {
namespace Operator {


int Operator :: d_iInstance_id=0;


Operator :: Operator(void)
{
    d_iObject_id                 = Operator::d_iInstance_id;

    d_iDebugPrintInfoLevel       = 0;   

    Operator :: d_iInstance_id++;
}
 

Operator :: Operator(const AMP::shared_ptr<OperatorParameters> & params)
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
 

void Operator :: reset(const AMP::shared_ptr<OperatorParameters>& params)
{
    AMP_INSIST( ((params.get()) != NULL), "NULL parameter" );

    // try and keep the next call the last in the function
    // so as not to override any parameters set through it 
    // by accident
    getFromInput(params->d_db);

}


void Operator :: residual(AMP::LinearAlgebra::Vector::const_shared_ptr f, 
			  AMP::LinearAlgebra::Vector::const_shared_ptr u,
			  AMP::LinearAlgebra::Vector::shared_ptr r)
{
    AMP_INSIST( ((u.get()) != NULL), "NULL Solution Vector" );
    AMP_INSIST( ((r.get()) != NULL), "NULL Residual Vector" );

    apply(u, r);

    AMP::LinearAlgebra::Vector::shared_ptr rInternal = subsetOutputVector(r);
    AMP_INSIST( (rInternal.get() != NULL), "rInternal is NULL" );

    // the rhs can be NULL
    if(f.get() != NULL)  {
      AMP::LinearAlgebra::Vector::const_shared_ptr fInternal = subsetOutputVector(f);
      rInternal->subtract( fInternal, rInternal );
    }
    else {
      rInternal->scale( -1.0 );
    }

    rInternal->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
}

void Operator :: getFromInput(const AMP::shared_ptr<AMP::Database>& db)
{
    AMP_INSIST( ((db.get()) != NULL), "NULL database" );

    d_iDebugPrintInfoLevel = db->getIntegerWithDefault("print_info_level", 0);
}


AMP::LinearAlgebra::Vector::shared_ptr  Operator::subsetOutputVector(AMP::LinearAlgebra::Vector::shared_ptr vec)
{
    PROFILE_START("subsetOutputVector",1);
    AMP::LinearAlgebra::Variable::shared_ptr var = getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr retvec;
    if(d_Mesh.get() != NULL) {
        AMP::LinearAlgebra::VS_Mesh meshSelector(d_Mesh);
        AMP::LinearAlgebra::Vector::shared_ptr meshSubsetVec = vec->select(meshSelector, vec->getVariable()->getName());
        AMP::LinearAlgebra::Vector::shared_ptr varSubsetVec = meshSubsetVec->subsetVectorForVariable(var);
        retvec = varSubsetVec;
    } else {
        retvec = vec->subsetVectorForVariable(var);
    }
    PROFILE_STOP("subsetOutputVector",1);
    return retvec;
}


AMP::LinearAlgebra::Vector::shared_ptr  Operator::subsetInputVector(AMP::LinearAlgebra::Vector::shared_ptr vec)
{
    PROFILE_START("subsetInputVector",1);
    AMP::LinearAlgebra::Variable::shared_ptr var = getInputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr retvec;
    if(d_Mesh.get() != NULL) {
        AMP::LinearAlgebra::VS_Mesh meshSelector(d_Mesh);
        AMP::LinearAlgebra::Vector::shared_ptr meshSubsetVec = vec->select(meshSelector, vec->getVariable()->getName());
        AMP::LinearAlgebra::Vector::shared_ptr varSubsetVec = meshSubsetVec->subsetVectorForVariable(var);
        retvec = varSubsetVec;
    } else {
        retvec = vec->subsetVectorForVariable(var);
    }
    PROFILE_STOP("subsetInputVector",1);
    return retvec;
}


AMP::LinearAlgebra::Vector::const_shared_ptr  Operator::subsetOutputVector(AMP::LinearAlgebra::Vector::const_shared_ptr vec)
{
    PROFILE_START("constSubsetOutputVector",1);
    AMP::LinearAlgebra::Variable::shared_ptr var = getOutputVariable();
    AMP::LinearAlgebra::Vector::const_shared_ptr retvec;
    if(d_Mesh.get() != NULL) {
        AMP::LinearAlgebra::VS_Mesh meshSelector(d_Mesh);
        AMP::LinearAlgebra::Vector::const_shared_ptr meshSubsetVec = vec->constSelect(meshSelector, vec->getVariable()->getName());
        AMP::LinearAlgebra::Vector::const_shared_ptr varSubsetVec = meshSubsetVec->constSubsetVectorForVariable(var);
        retvec = varSubsetVec;
    } else {
        retvec = vec->constSubsetVectorForVariable(var);
    }
    PROFILE_STOP("constSubsetOutputVector",1);
    return retvec;
}


AMP::LinearAlgebra::Vector::const_shared_ptr  Operator::subsetInputVector(AMP::LinearAlgebra::Vector::const_shared_ptr vec)
{
    PROFILE_START("constSubsetInputVector",1);
    AMP::LinearAlgebra::Variable::shared_ptr var = getInputVariable();
    AMP::LinearAlgebra::Vector::const_shared_ptr retvec;
    if(d_Mesh.get() != NULL) {
        AMP::LinearAlgebra::VS_Mesh meshSelector(d_Mesh);
        AMP::LinearAlgebra::Vector::const_shared_ptr meshSubsetVec = vec->constSelect(meshSelector, vec->getVariable()->getName());
        AMP::LinearAlgebra::Vector::const_shared_ptr varSubsetVec = meshSubsetVec->constSubsetVectorForVariable(var);
        retvec = varSubsetVec;
    } else {
        retvec = vec->constSubsetVectorForVariable(var);
    }
    PROFILE_STOP("constSubsetInputVector",1);
    return retvec;
}


}
}


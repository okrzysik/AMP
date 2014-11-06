#include "AsyncMapOperator.h"
#include "AsyncMapOperatorParameters.h"
#include "ampmesh/MultiMesh.h"
#include "utils/ProfilerApp.h"

namespace AMP {
namespace Operator {


AsyncMapOperator::AsyncMapOperator ( const AMP::shared_ptr <OperatorParameters> &p )
    : AsynchronousOperator ( p )
{
    // Fill some basic info
    AMP::shared_ptr<AsyncMapOperatorParameters> params = AMP::dynamic_pointer_cast<AsyncMapOperatorParameters>(p);
    d_MapComm = params->d_MapComm;
    d_mesh1 = params->d_Mesh1;
    d_mesh2 = params->d_Mesh2;
    AMP_INSIST( !d_MapComm.isNull(), "NULL communicator for map is invalid");
    AMP_INSIST( d_MapComm.sumReduce<int>(d_mesh1.get()!=NULL?1:0)>0, "Somebody must own mesh 1");
    AMP_INSIST( d_MapComm.sumReduce<int>(d_mesh2.get()!=NULL?1:0)>0, "Somebody must own mesh 2");
    // Create a multimesh to use for the operator base class for subsetting
    std::vector<AMP::Mesh::Mesh::shared_ptr> meshes;
    if ( d_mesh1.get() != NULL )
        meshes.push_back( d_mesh1 );
    if ( d_mesh2.get() != NULL )
        meshes.push_back( d_mesh2 );
    d_Mesh = AMP::shared_ptr<AMP::Mesh::MultiMesh>(new AMP::Mesh::MultiMesh(d_MapComm,meshes));
    // Get the input variable
    bool var = params->d_db->keyExists("VariableName");
    bool var1 = params->d_db->keyExists("VariableName1");
    bool var2 = params->d_db->keyExists("VariableName2");
    AMP_INSIST(var1||var2||var,"VariableName must exist in database");
    if ( var ) {
        AMP_INSIST(!var1&&!var2,"VariableName is used, VariableName1 and VariableName2cannot be used");
        std::string variableName = params->d_db->getString("VariableName");
        d_inpVariable = AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable(variableName) );
        d_outVariable = AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable(variableName) );
    } else {
        AMP_INSIST(var1&&var2,"Both VariableName1 and VariableName2 must be used");
        std::string variableName1 = params->d_db->getString("VariableName1");
        std::string variableName2 = params->d_db->getString("VariableName2");
        d_inpVariable = AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable(variableName1) );
        d_outVariable = AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable(variableName2) );
    }
}


AsyncMapOperator::~AsyncMapOperator ()
{
}


void AsyncMapOperator::apply(AMP::LinearAlgebra::Vector::const_shared_ptr f,
        AMP::LinearAlgebra::Vector::const_shared_ptr u, 
        AMP::LinearAlgebra::Vector::shared_ptr r,
        const double a, const double b)
{
    PROFILE_START("apply");
    applyStart  ( f , u , r , a , b );
    applyFinish ( f , u , r , a , b );
    if ( requiresMakeConsistentSet() ) {
        AMP_ASSERT(d_OutputVector.get()!=NULL);
        d_OutputVector->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    }
    PROFILE_STOP("apply");
}


bool AsyncMapOperator::requiresMakeConsistentSet()
{ 
    return false;
}

AMP::Mesh::Mesh::shared_ptr AsyncMapOperator::getMesh(int which) 
{
  if(which == 1) {
    return d_mesh1;
  } else if(which == 2) {
    return d_mesh2;
  } else {
    AMP_ERROR("Wrong option!");
    return AMP::Mesh::Mesh::shared_ptr();
  }
}

}
}


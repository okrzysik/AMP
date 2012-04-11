#include "AsyncMapOperator.h"
#include "AsyncMapOperatorParameters.h"
#include "ampmesh/MultiMesh.h"

namespace AMP {
namespace Operator {


AsyncMapOperator::AsyncMapOperator ( const boost::shared_ptr <OperatorParameters> &p )
    : AsynchronousOperator ( p )
{
    // Fill some basic info
    boost::shared_ptr<AsyncMapOperatorParameters> params = boost::dynamic_pointer_cast<AsyncMapOperatorParameters>(p);
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
    d_Mesh = boost::shared_ptr<AMP::Mesh::MultiMesh>(new AMP::Mesh::MultiMesh(d_MapComm,meshes));
    // Get the input variable
    std::string variableName = params->d_db->getString("VariableName");
    d_inpVariable = AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable(variableName) );
}


AsyncMapOperator::~AsyncMapOperator ()
{
}


void AsyncMapOperator::apply(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const  AMP::LinearAlgebra::Vector::shared_ptr &u, 
        AMP::LinearAlgebra::Vector::shared_ptr  &r,
        const double a, const double b)
{
    applyStart  ( f , u , r , a , b );
    applyFinish ( f , u , r , a , b );
    if ( requiresMakeConsistentSet() ) {
        AMP_ASSERT(d_OutputVector.get()!=NULL);
        d_OutputVector->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    }
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


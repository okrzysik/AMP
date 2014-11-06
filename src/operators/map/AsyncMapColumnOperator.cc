#include "AsyncMapOperator.h"
#include "AsyncMapColumnOperator.h"
#include "utils/ProfilerApp.h"

namespace AMP {
namespace Operator {


size_t globalMapTagOffset = 0;      // Initialize the global map tag offset


AsyncMapColumnOperator::AsyncMapColumnOperator ( const AMP::shared_ptr<OperatorParameters> & params )
    : AsynchronousColumnOperator ( params )
{
}


void  AsyncMapColumnOperator::setVector ( AMP::LinearAlgebra::Vector::shared_ptr p )
{
    d_OutputVector = p;
    for (size_t i=0; i<d_Operators.size(); i++)
        AMP::dynamic_pointer_cast<AsyncMapOperator>(d_Operators[i])->setVector ( d_OutputVector );
}


void  AsyncMapColumnOperator::append ( AMP::shared_ptr < Operator > op )
{
    AMP::shared_ptr<AsyncMapColumnOperator>  mapColumn = AMP::dynamic_pointer_cast<AsyncMapColumnOperator> ( op );
    if ( mapColumn )
    {
      std::vector< AMP::shared_ptr < Operator > >::iterator curOp = mapColumn.get()->d_Operators.begin();
      while ( curOp != mapColumn.get()->d_Operators.end() )
      {
        append ( *curOp );
        ++curOp;
      }
    }
    else
    {
      AMP::shared_ptr<AsyncMapOperator>  mapOp = AMP::dynamic_pointer_cast<AsyncMapOperator> ( op );
      AMP_INSIST ( mapOp , "Attempt to add a non-AsyncMapOperator to a AsyncMapColumnOperator" );
      AsynchronousColumnOperator::append ( mapOp );
    }
}


void AsyncMapColumnOperator::apply(AMP::LinearAlgebra::Vector::const_shared_ptr f,
        AMP::LinearAlgebra::Vector::const_shared_ptr u, 
        AMP::LinearAlgebra::Vector::shared_ptr r,
        const double a, const double b)
{
    PROFILE_START("apply");
    this->applyStart  ( f , u , r , a , b );
    this->applyFinish ( f , u , r , a , b );
    if ( requiresMakeConsistentSet() ) {
        AMP_ASSERT(d_OutputVector.get()!=NULL);
        d_OutputVector->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    }
    PROFILE_STOP("apply");
}


bool AsyncMapColumnOperator::requiresMakeConsistentSet()
{ 
    bool test = false;
    for (size_t i=0; i<d_Operators.size(); i++)
        test = test | AMP::dynamic_pointer_cast<AsyncMapOperator>(d_Operators[i])->requiresMakeConsistentSet();
    return test;
}


/********************************************************************
* Function to copy a key from database 1 to database 2              *
* If the key is an array of size N it will only copy the ith value. *
********************************************************************/
static void copyKey(AMP::shared_ptr<AMP::Database> &database1, 
    AMP::shared_ptr<AMP::Database> &database2, std::string key, int N, int i ) 
{
    int size = database1->getArraySize(key);
    AMP::Database::DataType type = database1->getArrayType(key);
    switch (type) {
        case AMP::Database::AMP_INVALID: {
            // We don't know what this is
            AMP_ERROR("Invalid database object");
            } break;
        case AMP::Database::AMP_DATABASE: {
            // Copy the database
            AMP_ERROR("Not programmed for databases yet");
            } break;
        case AMP::Database::AMP_BOOL: {
            // Copy a bool
            std::vector<unsigned char> data = database1->getBoolArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size )
                database2->putBool(key,data[i]);
            else
                database2->putBoolArray(key,data);
            } break;
        case AMP::Database::AMP_CHAR: {
            // Copy a char
            std::vector<char> data = database1->getCharArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size ) {
                database2->putChar(key,data[i]);
            } else {
                // We need to try a search and replace
                database2->putCharArray(key,data);
            }
            } break;
        case AMP::Database::AMP_INT: {
            // Copy an int
            std::vector<int> data = database1->getIntegerArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size )
                database2->putInteger(key,data[i]);
            else
                database2->putIntegerArray(key,data);
            } break;
        case AMP::Database::AMP_COMPLEX: {
            // Copy a complex number
            std::vector< std::complex<double> > data = database1->getComplexArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size )
                database2->putComplex(key,data[i]);
            else
                database2->putComplexArray(key,data);
            } break;
        case AMP::Database::AMP_DOUBLE: {
            // Copy a double number
            std::vector<double> data = database1->getDoubleArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size )
                database2->putDouble(key,data[i]);
            else
                database2->putDoubleArray(key,data);
            } break;
        case AMP::Database::AMP_FLOAT: {
            // Copy a float
            std::vector<float> data = database1->getFloatArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size )
                database2->putFloat(key,data[i]);
            else
                database2->putFloatArray(key,data);
            } break;
        case AMP::Database::AMP_STRING: {
            // Copy a string
            std::vector<std::string> data = database1->getStringArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size )
                database2->putString(key,data[i]);
            else
                database2->putStringArray(key,data);
            } break;
        case AMP::Database::AMP_BOX: {
            // Copy a box
            AMP_ERROR("Not programmed for boxes yet");
            } break;
        default:
            AMP_ERROR("Unknown key type");
    }
}


/************************************************************
* Function to create the databases for the individual maps  *
************************************************************/
std::vector<AMP::shared_ptr<AMP::Database> >  AsyncMapColumnOperator::createDatabases(AMP::shared_ptr<AMP::Database> database1)
{
    AMP_INSIST(database1->keyExists("N_maps"),"N_maps must exist in input database");
    int N_maps = database1->getInteger("N_maps");
    // Create the basic databases for each mesh
    std::vector<AMP::shared_ptr<AMP::Database> > meshDatabases;
    meshDatabases.reserve(N_maps);
    for (int i=0; i<N_maps; i++) {
        // Create a new database from the existing database
        AMP::shared_ptr<AMP::Database> database2( new AMP::MemoryDatabase("MapDatabase") );
        std::vector<std::string> keys = database1->getAllKeys();
        for (size_t k=0; k<keys.size(); k++) {
            if ( keys[k].compare("N_maps")==0 ) {
                // These keys are used by the builder and should not be in the sub database
            } else {
                // We need to copy the key
                copyKey(database1,database2,keys[k],N_maps,i);
            }
        }
        meshDatabases.push_back(database2);
    }
    return meshDatabases;
}


}
}


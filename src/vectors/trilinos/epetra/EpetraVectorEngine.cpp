#include "AMP/vectors/trilinos/epetra/EpetraVectorEngine.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/data/VectorDataCPU.h"


namespace AMP {
namespace LinearAlgebra {


static inline Epetra_Vector &getEpetraVector( Vector &vec )
{
    auto epetra = dynamic_cast<EpetraVectorEngine *>( &vec );
    AMP_INSIST( epetra != nullptr, "Not an EpetraVectorEngine" );
    return epetra->getEpetra_Vector();
}


/********************************************************
 * Constructor                                           *
 ********************************************************/
EpetraVectorEngine::EpetraVectorEngine( std::shared_ptr<EpetraVectorEngineParameters> alias,
                                        std::shared_ptr<VectorData> buf )
    : d_Params( alias )
{
    d_DOFManager = alias->d_DOFManager;
    d_VectorOps  = std::make_shared<EpetraVectorOperations>();
    d_VectorData = EpetraVectorData::create( alias, buf );
}

Vector::shared_ptr EpetraVectorEngine::cloneVector( const Variable::shared_ptr name ) const
{
    auto params = std::dynamic_pointer_cast<EpetraVectorEngineParameters>( d_Params );
    auto buffer = std::make_shared<VectorDataCPU<double>>(
        params->beginDOF(), params->getLocalSize(), params->getGlobalSize() );

    auto retVal = std::make_shared<EpetraVectorEngine>( d_Params, buffer );
    retVal->setVariable( name );
    return retVal;
}

void EpetraVectorEngine::swapVectors( Vector &other )
{
    double *my_pointer;
    double *oth_pointer;
    getEpetra_Vector().ExtractView( &my_pointer );
    getEpetraVector( other ).ExtractView( &oth_pointer );
    getEpetraVector( other ).ResetView( my_pointer );
    getEpetra_Vector().ResetView( oth_pointer );
}

} // namespace LinearAlgebra
} // namespace AMP

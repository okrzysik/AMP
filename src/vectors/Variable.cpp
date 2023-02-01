#include "AMP/vectors/Variable.h"
#include "AMP/IO/RestartManager.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/VectorSelector.h"


namespace AMP::LinearAlgebra {


Variable::Variable( const std::string &name ) : d_VariableName( name ) {}


Variable::~Variable() = default;


std::shared_ptr<Variable> Variable::cloneVariable( const std::string &name ) const
{
    return std::make_shared<Variable>( name );
}


bool Variable::operator==( const Variable &rhs ) const
{
    return d_VariableName == rhs.d_VariableName;
}


bool Variable::operator!=( const Variable &rhs ) const { return !( *this == rhs ); }


void Variable::setUnits( const Units &t ) { d_Units = t; }


const Units &Variable::getUnits() const { return d_Units; }


std::shared_ptr<VectorSelector> Variable::createVectorSelector() const
{
    return std::make_shared<VS_ByVariableName>( d_VariableName );
}


uint64_t Variable::getID() const { return reinterpret_cast<uint64_t>( this ); }


void Variable::writeRestart( int64_t fid ) const
{
    writeHDF5( fid, "type", className() );
    writeHDF5( fid, "name", d_VariableName );
    writeHDF5( fid, "units", d_Units );
}

Variable::Variable( int64_t fid )
{
    readHDF5( fid, "name", d_VariableName );
    readHDF5( fid, "units", d_Units );
}


} // namespace AMP::LinearAlgebra


/********************************************************
 *  Restart operations                                   *
 ********************************************************/
template<>
AMP::IO::RestartManager::DataStoreType<AMP::LinearAlgebra::Variable>::DataStoreType(
    const std::string &name,
    std::shared_ptr<const AMP::LinearAlgebra::Variable> var,
    RestartManager *manager )
    : d_data( var )
{
    d_name        = name;
    d_hash        = var->getID();
    auto multivar = std::dynamic_pointer_cast<const AMP::LinearAlgebra::MultiVariable>( var );
    if ( multivar ) {
        for ( auto var2 : *multivar )
            manager->registerData( var2 );
    }
}
template<>
void AMP::IO::RestartManager::DataStoreType<AMP::LinearAlgebra::Variable>::write(
    hid_t fid, const std::string &name ) const
{
    hid_t gid = createGroup( fid, name );
    AMP_ERROR( "Not finished" );
    closeGroup( gid );
}
template<>
std::shared_ptr<AMP::LinearAlgebra::Variable>
AMP::IO::RestartManager::getData<AMP::LinearAlgebra::Variable>( const std::string &name )
{
    std::shared_ptr<AMP::LinearAlgebra::Variable> var;
    hid_t gid = openGroup( d_fid, name );
    std::string type;
    readHDF5( gid, "type", type );
    if ( type == "Variable" ) {
        var = std::make_shared<AMP::LinearAlgebra::Variable>( gid );
    } else {
        AMP_ERROR( "Unknown variable type" );
    }
    closeGroup( gid );
    return var;
}

#include "AMP/vectors/Variable.h"
#include "AMP/IO/RestartManager.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/VectorSelector.h"


namespace AMP::LinearAlgebra {


Variable::Variable( const std::string &name ) : d_VariableName( name ) {}


Variable::~Variable() = default;


std::shared_ptr<Variable> Variable::clone() const { return clone( d_VariableName ); }

std::shared_ptr<Variable> Variable::clone( const std::string &name ) const
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


/********************************************************
 *  Restart operations                                   *
 ********************************************************/
void Variable::registerChildObjects( AMP::IO::RestartManager * ) const {}
void Variable::writeRestart( int64_t fid ) const
{
    IO::writeHDF5( fid, "type", className() );
    IO::writeHDF5( fid, "name", d_VariableName );
    IO::writeHDF5( fid, "units", d_Units );
}
Variable::Variable( int64_t fid )
{
    IO::readHDF5( fid, "name", d_VariableName );
    IO::readHDF5( fid, "units", d_Units );
}


} // namespace AMP::LinearAlgebra


/********************************************************
 *  Restart operations                                   *
 ********************************************************/
template<>
AMP::IO::RestartManager::DataStoreType<AMP::LinearAlgebra::Variable>::DataStoreType(
    std::shared_ptr<const AMP::LinearAlgebra::Variable> var, RestartManager *manager )
    : d_data( var )
{
    d_hash        = var->getID();
    auto multivar = std::dynamic_pointer_cast<const AMP::LinearAlgebra::MultiVariable>( var );
    if ( multivar ) {
        for ( auto var2 : *multivar )
            manager->registerObject( var2 );
    }
}
template<>
void AMP::IO::RestartManager::DataStoreType<AMP::LinearAlgebra::Variable>::write(
    hid_t fid, const std::string &name ) const
{
    hid_t gid = IO::createGroup( fid, name );
    d_data->writeRestart( gid );
    IO::closeGroup( gid );
}
template<>
std::shared_ptr<AMP::LinearAlgebra::Variable>
AMP::IO::RestartManager::DataStoreType<AMP::LinearAlgebra::Variable>::read(
    hid_t fid, const std::string &name, RestartManager *manager ) const
{
    std::shared_ptr<AMP::LinearAlgebra::Variable> var;
    hid_t gid = IO::openGroup( fid, name );
    std::string type;
    IO::readHDF5( gid, "type", type );
    if ( type == "Variable" ) {
        var = std::make_shared<AMP::LinearAlgebra::Variable>( gid );
    } else if ( type == "MultiVariable" ) {
        var = std::make_shared<AMP::LinearAlgebra::MultiVariable>( gid, manager );
    } else {
        AMP_ERROR( "Unknown variable type: " + type );
    }
    IO::closeGroup( gid );
    return var;
}

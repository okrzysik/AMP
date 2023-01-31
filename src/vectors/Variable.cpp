#include "AMP/vectors/Variable.h"
#include "AMP/IO/RestartManager.h"
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


} // namespace AMP::LinearAlgebra


/********************************************************
 *  Restart operations                                   *
 ********************************************************/
template<>
AMP::IO::RestartManager::DataStoreType<AMP::LinearAlgebra::Variable>::DataStoreType(
    const std::string &name,
    std::shared_ptr<const AMP::LinearAlgebra::Variable> data,
    RestartManager *manager )
    : d_data( data )
{
    d_name = name;
    d_hash = data->getID();
    AMP_ERROR( "Not finished" );
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
    hid_t gid = openGroup( d_fid, name );
    AMP_ERROR( "Not finished" );
    closeGroup( gid );
}

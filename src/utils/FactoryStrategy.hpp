#ifndef AMP_FactoryStrategy_HPP_
#define AMP_FactoryStrategy_HPP_

#include "AMP/utils/Database.h"
#include "AMP/utils/UtilityMacros.h"
#include "AMP/utils/shared_ptr.h"
#include <map>
#include <sstream>
#include <string>

namespace AMP {

// template for factories for operators and solvers

// NOTE: ideally the create function should return a
// unique_ptr but since AMP does not currently have that
// we will return a shared_ptr

template<typename TYPE, typename PARAMETERS>
class FactoryStrategy
{
private:
    /**
     * Private constructor.
     */
    FactoryStrategy();

    using FunctionPtr = AMP::shared_ptr<TYPE> ( * )( AMP::shared_ptr<PARAMETERS> parameters );
    using FunctionMap = std::map<std::string, FunctionPtr>;

    FunctionMap d_factories; //! maps from names to factories

    //! given name and parameter object create object
    AMP::shared_ptr<TYPE> create( std::string name, AMP::shared_ptr<PARAMETERS> parameters );

    //! register a function that creates an object
    void registerFunction( std::string name, FactoryStrategy::FunctionPtr ptr );

public:
    /**
     * Destructor.
     */
    ~FactoryStrategy();

    //! get a singleton instance of the factory
    static FactoryStrategy &getFactory();

    /**
     * Factory method for generating objects with characteristics
     * specified by parameters.
     */
    static AMP::shared_ptr<TYPE> create( AMP::shared_ptr<PARAMETERS> parameters );

    // public interface for registering a function that creates an object
    static void registerFactory( std::string name, FactoryStrategy::FunctionPtr ptr );
};

template<typename TYPE, typename PARAMETERS>
FactoryStrategy<TYPE, PARAMETERS>::FactoryStrategy()
{
}

template<typename TYPE, typename PARAMETERS>
FactoryStrategy<TYPE, PARAMETERS>::~FactoryStrategy()
{
}

template<typename TYPE, typename PARAMETERS>
FactoryStrategy<TYPE, PARAMETERS> &FactoryStrategy<TYPE, PARAMETERS>::getFactory( void )
{
    static FactoryStrategy<TYPE, PARAMETERS> singletonInstance;
    return singletonInstance;
}

template<typename TYPE, typename PARAMETERS>
AMP::shared_ptr<TYPE>
FactoryStrategy<TYPE, PARAMETERS>::create( std::string name,
                                           AMP::shared_ptr<PARAMETERS> parameters )
{
    AMP::shared_ptr<TYPE> obj;
    auto it = d_factories.find( name );
    if ( it != d_factories.end() ) {
        obj = it->second( parameters );
    } else {
        std::stringstream err;
        err << "Unable to create object " << name;
        AMP_ERROR( err.str() );
    }

    return obj;
}

template<typename TYPE, typename PARAMETERS>
AMP::shared_ptr<TYPE>
FactoryStrategy<TYPE, PARAMETERS>::create( AMP::shared_ptr<PARAMETERS> parameters )
{
    AMP_ASSERT( parameters != nullptr );

    std::string objectName = "";

    AMP::shared_ptr<AMP::Database> inputDatabase = parameters->d_db;

    if ( inputDatabase->keyExists( "name" ) ) {
        objectName = inputDatabase->getString( "name" );
    } else {
        AMP_ERROR( "FactoryStrategy"
                   << " -- Required key `name'"
                   << " missing in input." );
    }

    auto &factory = FactoryStrategy<TYPE, PARAMETERS>::getFactory();
    return factory.create( objectName, parameters );
}

template<typename TYPE, typename PARAMETERS>
void FactoryStrategy<TYPE, PARAMETERS>::registerFactory(
    std::string name, FactoryStrategy<TYPE, PARAMETERS>::FunctionPtr ptr )
{
    auto &factory = FactoryStrategy<TYPE, PARAMETERS>::getFactory();
    factory.registerFunction( name, ptr );
}

template<typename TYPE, typename PARAMETERS>
void FactoryStrategy<TYPE, PARAMETERS>::registerFunction(
    std::string name, FactoryStrategy<TYPE, PARAMETERS>::FunctionPtr ptr )
{
    d_factories[name] = ptr;
}
} // namespace AMP
#endif

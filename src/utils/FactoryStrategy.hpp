#ifndef AMP_FactoryStrategy_HPP_
#define AMP_FactoryStrategy_HPP_

#include <memory>
#include <map>
#include <string>
#include <sstream>
#include "utils/InputDatabase.h"
#include "utils/UtilityMacros.h"

namespace AMP{

// template for factories for operators and solvers
template <typename TYPE, typename PARAMETERS>
class FactoryStrategy
{
private:
    /**
     * Private constructor.
     */
    FactoryStrategy();

    using FunctionPtr = std::unique_ptr<TYPE>(*)( std::shared_ptr<PARAMETERS> parameters );
    using FunctionMap = std::map<std::string, FunctionPtr>;

    FunctionMap d_factories; //! maps from names to factories

    //! given name and parameter object create object
    std::unique_ptr<TYPE> create( std::string name, std::shared_ptr<PARAMETERS> parameters );
    
    //! register a function that creates an object
    void registerFunction(std::string name, FactoryStrategy::FunctionPtr ptr);

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
    static std::unique_ptr<TYPE>
    create( std::shared_ptr<PARAMETERS> parameters );

    // public interface for registering a function that creates an object
    static void registerFactory(std::string name, FactoryStrategy::FunctionPtr ptr);

};

template <typename TYPE, typename PARAMETERS>
FactoryStrategy<TYPE,PARAMETERS>::FactoryStrategy()
{
}

template <typename TYPE, typename PARAMETERS>
FactoryStrategy<TYPE,PARAMETERS>::~FactoryStrategy()
{    
}

template <typename TYPE, typename PARAMETERS>
FactoryStrategy<TYPE,PARAMETERS> &FactoryStrategy<TYPE,PARAMETERS>::getFactory(void )
{
    static FactoryStrategy<TYPE,PARAMETERS> singletonInstance;
    return singletonInstance;
}

template <typename TYPE, typename PARAMETERS>
std::unique_ptr<TYPE> FactoryStrategy<TYPE,PARAMETERS>::create(std::string name, 
                                                               std::shared_ptr<PARAMETERS> parameters)
{
    std::unique_ptr<TYPE> obj;
    auto it = d_factories.find(name);
    if( it != d_factories.end() ) {
        obj = it->second(parameters);
    } else {
        std::stringstream err;
        err << "Unable to create object " << name;
        SAMR_ERROR(err.str());
    }

    return obj;
}

template <typename TYPE, typename PARAMETERS>
std::unique_ptr<TYPE>
FactoryStrategy<TYPE, PARAMETERS>::create( std::shared_ptr<PARAMETERS> parameters )
{
    SAMR_ASSERT( parameters != nullptr );

    std::string objectName = "";

    std::shared_ptr<SAMRAI::tbox::Database> inputDatabase = parameters->d_db;

    if ( inputDatabase->keyExists( "name" ) ) {
        objectName = inputDatabase->getString( "name" );
    } else {
        SAMR_ERROR( "FactoryStrategy"
                    << " -- Required key `name'"
                    << " missing in input." );
    }

    auto &factory = FactoryStrategy<TYPE, PARAMETERS>::getFactory();
    return factory.create(objectName, parameters);
}

template <typename TYPE, typename PARAMETERS>
void FactoryStrategy<TYPE, PARAMETERS>::registerFactory(std::string name, 
                                                        FactoryStrategy<TYPE, PARAMETERS>::FunctionPtr ptr)
{
    auto &factory = FactoryStrategy<TYPE, PARAMETERS>::getFactory();
    factory.registerFunction(name, ptr);
}

template <typename TYPE, typename PARAMETERS>
void FactoryStrategy<TYPE, PARAMETERS>::registerFunction(std::string name, 
                                                         FactoryStrategy<TYPE, PARAMETERS>::FunctionPtr ptr)
{
    d_factories[name] = ptr;
}


}
#endif

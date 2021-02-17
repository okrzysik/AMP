#ifndef FACTORY_H
#define FACTORY_H

/**
 *    \brief        This provides a templated method of creating self registering factories
 *                Each of these factories must take the same number/type of constrctor args
 *    \author        Nathan Barnett
 */


#include "AMP/utils/Singleton.h"
#include "AMP/utils/UtilityMacros.h"
#include "AMP/utils/getkeys.h"

#include <map>
#include <memory>
#include <string>


namespace AMP {
namespace voodoo {

/**
 *    \class        RegistrationBase
 *    \brief        We create one registration object per base class. This is a template for the
 * base class of the
 * registration objects
 */
template<class BaseClass, typename Key = std::string, unsigned long n = 0>
struct RegistrationBase {
};

/**
 *    \class        Registration
 *    \brief        Template for the registration objects
 */
template<class BaseClass, class Derived, typename Key = std::string, unsigned long n = 0>
struct Registration : public RegistrationBase<BaseClass, Key, n> {
};

/**
 *    \class        Factory
 *    \brief        factory class is a singleton class
 */
template<class BaseClass, typename Key = std::string, unsigned long n = 0>
class Factory : public Singleton<Factory<BaseClass, Key, n>>
{
};


template<class BaseClass, typename Key>
struct RegistrationBase<BaseClass, Key, 0> {
    virtual std::shared_ptr<BaseClass> create() const = 0;
    virtual ~RegistrationBase() {}
};


template<class BaseClass, typename Key>
class Factory<BaseClass, Key, 0> : public Singleton<Factory<BaseClass, Key, 0>>
{
public:
    friend class Singleton<Factory<BaseClass, Key, 0>>;
    bool Register( const Key &key, std::shared_ptr<RegistrationBase<BaseClass, Key, 0>> reg )
    {
        return regMapPtr->insert( typename RegistrationMap::value_type( key, reg ) ).second;
    }

    bool Register( const Factory<BaseClass, Key, 0> &other )
    {
        std::shared_ptr<RegistrationMap> tmp( new RegistrationMap( *regMapPtr ) );
        RegistrationMapIterator iter( other.regMapPtr->begin() );
        while ( iter != other.regMapPtr->end() ) {
            if ( !tmp->insert( typename RegistrationMap::value_type( *iter++ ) ).second )
                return false;
        }
        std::swap( regMapPtr, tmp );
        return true;
    }

    std::shared_ptr<BaseClass> create( const Key &id ) const
    {
        RegistrationMapIterator iter( regMapPtr->find( id ) );
        if ( iter == regMapPtr->end() )
            AMP_ERROR( "Unregistered creator" );
        return iter->second->create();
    }

    std::vector<Key> getKeys() const { return AMP::voodoo::getKeys( *regMapPtr ); }

private:
    typedef typename std::map<Key, std::shared_ptr<RegistrationBase<BaseClass, Key, 0>>>
        RegistrationMap;
    typedef typename RegistrationMap::const_iterator RegistrationMapIterator;
    std::shared_ptr<RegistrationMap> regMapPtr;
    Factory() : regMapPtr( new RegistrationMap ) {}
};


template<class BaseClass, class Derived, typename Key>
struct Registration<BaseClass, Derived, Key, 0> : public RegistrationBase<BaseClass, Key, 0> {
    explicit Registration( const Key &key )
    {
        try {
            Factory<BaseClass, Key, 0>::instance().Register(
                key,
                std::shared_ptr<RegistrationBase<BaseClass, Key, 0>>(
                    new Registration<BaseClass, Derived, Key, 0> ) );
        } catch ( ... ) {
        }
    }
    virtual ~Registration() {}

private:
    std::shared_ptr<BaseClass> create() const override
    {
        return std::shared_ptr<BaseClass>( new Derived() );
    }
    Registration() {}
    Registration( const Registration<BaseClass, Derived, Key, 0> & );
    Registration<BaseClass, Derived, Key, 0> &
    operator=( const Registration<BaseClass, Derived, Key, 0> & );
};


} // namespace voodoo
} // namespace AMP


#endif // FACTORY_H

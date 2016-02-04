#ifndef FACTORY_H
#define FACTORY_H

/**
 *    \brief        This provides a templated method of creating self registering factories
 *                Each of these factories must take the same number/type of constrctor args
 *    \author        Nathan Barnett
 */


#include <map>
#include <string>

#include "Singleton.h"
#include "Utilities.h"
#include "getkeys.h"

#include "utils/shared_ptr.h"


namespace AMP {
namespace voodoo {

/**
 *    \class        RegistrationBase
 *    \brief        We create one registration object per base class. This is a template for the
 * base class of the
 * registration objects
 */
template <class BaseClass, typename Key = std::string, unsigned long n = 0>
struct RegistrationBase {
};

/**
 *    \class        Registration
 *    \brief        Template for the registration objects
 */
template <class BaseClass, class Derived, typename Key = std::string, unsigned long n = 0>
struct Registration : public RegistrationBase<BaseClass, Key, n> {
};

/**
 *    \class        Factory
 *    \brief        factory class is a singleton class
 */
template <class BaseClass, typename Key = std::string, unsigned long n = 0>
class Factory : public Singleton<Factory<BaseClass, Key, n>>
{
};


template <class BaseClass, typename Key>
struct RegistrationBase<BaseClass, Key, 0> {
    virtual AMP::shared_ptr<BaseClass> create() const = 0;
    virtual ~RegistrationBase() {}
};


template <class BaseClass, typename Key>
class Factory<BaseClass, Key, 0> : public Singleton<Factory<BaseClass, Key, 0>>
{
public:
    friend class Singleton<Factory<BaseClass, Key, 0>>;
    bool Register( const Key &key, AMP::shared_ptr<RegistrationBase<BaseClass, Key, 0>> reg )
    {
        return regMapPtr->insert( typename RegistrationMap::value_type( key, reg ) ).second;
    }

    bool Register( const Factory<BaseClass, Key, 0> &other )
    {
        AMP::shared_ptr<RegistrationMap> tmp( new RegistrationMap( *regMapPtr ) );
        RegistrationMapIterator iter( other.regMapPtr->begin() );
        while ( iter != other.regMapPtr->end() ) {
            if ( !tmp->insert( typename RegistrationMap::value_type( *iter++ ) ).second )
                return false;
        }
        std::swap( regMapPtr, tmp );
        return true;
    }

    AMP::shared_ptr<BaseClass> create( const Key &id ) const
    {
        RegistrationMapIterator iter( regMapPtr->find( id ) );
        if ( iter == regMapPtr->end() )
            AMP_ERROR( "Unregistered creator" );
        return iter->second->create();
    }

    std::vector<Key> getKeys() const { return AMP::voodoo::getKeys( *regMapPtr ); }

private:
    typedef typename std::map<Key, AMP::shared_ptr<RegistrationBase<BaseClass, Key, 0>>>
        RegistrationMap;
    typedef typename RegistrationMap::const_iterator RegistrationMapIterator;
    AMP::shared_ptr<RegistrationMap> regMapPtr;
    Factory() : regMapPtr( new RegistrationMap ) {}
};


template <class BaseClass, class Derived, typename Key>
struct Registration<BaseClass, Derived, Key, 0> : public RegistrationBase<BaseClass, Key, 0> {
    Registration( const Key &key )
    {
        try {
            Factory<BaseClass, Key, 0>::instance().Register(
                key,
                AMP::shared_ptr<RegistrationBase<BaseClass, Key, 0>>(
                    new Registration<BaseClass, Derived, Key, 0> ) );
        } catch ( ... ) {
        }
    }
    virtual ~Registration() {}
private:
    AMP::shared_ptr<BaseClass> create() const
    {
        return AMP::shared_ptr<BaseClass>( new Derived() );
    }
    Registration() {}
    Registration( const Registration<BaseClass, Derived, Key, 0> & );
    Registration<BaseClass, Derived, Key, 0> &
    operator=( const Registration<BaseClass, Derived, Key, 0> & );
};


} // namespace voodoo
} // namespace amp


#endif // FACTORY_H

#ifndef FACTORY_H
#define FACTORY_H

#include "AMP/utils/Singleton.h"
#include "AMP/utils/UtilityMacros.h"
#include "AMP/utils/getkeys.h"

#include <map>
#include <memory>
#include <string>


namespace AMP::voodoo {

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
        return d_map.insert( std::pair( key, reg ) ).second;
    }

    bool Register( const Factory<BaseClass, Key, 0> &other )
    {
        std::map<Key, std::shared_ptr<RegistrationBase<BaseClass, Key, 0>>> tmp;
        auto iter = other.d_map.begin();
        while ( iter != other.d_map.end() ) {
            if ( !tmp.insert( *iter++ ).second )
                return false;
        }
        std::swap( d_map, tmp );
        return true;
    }

    std::shared_ptr<BaseClass> create( const Key &id ) const
    {
        auto iter = d_map.find( id );
        if ( iter == d_map.end() )
            AMP_ERROR( "Unregistered creator" );
        return iter->second->create();
    }

    std::vector<Key> getKeys() const { return AMP::voodoo::getKeys( d_map ); }

    void clear() { return d_map.clear(); }

private:
    std::map<Key, std::shared_ptr<RegistrationBase<BaseClass, Key, 0>>> d_map;
    Factory() {}
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


} // namespace AMP::voodoo


#endif // FACTORY_H

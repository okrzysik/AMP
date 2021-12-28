#ifndef SINGLETON_H
#define SINGLETON_H

/**
 *	\brief		Provides a templated Meyer's singleton class
 *	\author		Nathan Barnett
 */


namespace AMP::voodoo {


template<class Object>
class Singleton
{
public:
    static Object &instance()
    {
        static Object object;
        return object;
    }

protected:
    Singleton() {}
    ~Singleton() {}

private: // emphasize the following members are private
    Singleton( const Singleton & );
    const Singleton &operator=( const Singleton & );
};


} // namespace AMP::voodoo

#endif // SINGLETON_H

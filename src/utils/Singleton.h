#ifndef SINGLETON_H
#define SINGLETON_H

/**
 *	\file 		singleton.h
 *	\brief		Provides a templated Meyer's singleton class
 *	\author		Nathan Barnett
 */



namespace AMP {
namespace voodoo 
{


template<class Object>
struct Singleton
{
public:
    static Object& instance()
    {
        static Object object;
        return object;
    }
protected:
      Singleton() {}
      ~Singleton() {}
private:  // emphasize the following members are private
      Singleton( const Singleton& );
      const Singleton& operator=( const Singleton& );
};


} // namespace VoodooFactory
} // namespace amp

#endif // SINGLETON_H






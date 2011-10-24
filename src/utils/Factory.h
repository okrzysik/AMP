#ifndef FACTORY_H
#define FACTORY_H

/**
 *	\file 		factory.h
 *	\brief		This provides a templated method of creating self registering factories
 *				Each of these factories must take the same number/type of constrctor args
 *	\author		Nathan Barnett
 */


#include <string>
#include <map>

#include "Singleton.h"
#include "getkeys.h"
#include "Utilities.h"

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include <boost/mpl/size.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/at.hpp>

#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>


namespace AMP 
{
	namespace voodoo
	{
		/**
		 *	\class		RegistrationBase
		 *	\brief		We create one registration object per base class. This is a template for the base class of the registration objects
		 */
		template 
		< class BaseClass, typename Key = std::string, unsigned long  = boost::mpl::size<typename BaseClass::factory_constructor_typeList>::value >
		struct RegistrationBase
		{
		};

		/**
		 *	\class		Registration
		 *	\brief		Template for the registration objects
		 */
		template
		< class BaseClass, class Derived, typename Key = std::string, unsigned long n = boost::mpl::size<typename BaseClass::factory_constructor_typeList>::value >
		struct Registration : public RegistrationBase<BaseClass, Key, n>
		{
		};

		/**
		 *	\class		Factory
		 *	\brief		factory class is a singleton class
		 */
		template
		< class BaseClass, typename Key = std::string , unsigned long  n = boost::mpl::size<typename BaseClass::factory_constructor_typeList>::value >
		class Factory : public Singleton< Factory<BaseClass, Key, n> > 
		{
		};

		// Create a typedef for each - named Parameter+unique_number (e.g. Param12)
		#define PARAM(z, TypeListNumber, data)  typename boost::mpl::at_c< typename BaseClass::factory_constructor_typeList, TypeListNumber>::type Param##TypeListNumber

		#define  REGISTRATION_BASE(z, TypeListNumber, data)														\
		template<class BaseClass, typename Key>  																\
		struct RegistrationBase<BaseClass, Key, TypeListNumber> 												\
		{																										\
			virtual boost::shared_ptr<BaseClass> create(BOOST_PP_ENUM(TypeListNumber, PARAM, ~)) const=0;		\
		};

		// Currently, BOOST_PP_LIMIT_REPEAT=256, which is the maximum number of factories that can be created.
		// This could be increased by making a nested call or increasing the number in boost/preprocessor/config/limits.hpp
		#define AMP_FACTORY_REPEAT 1              
		
		// First number is how many times the preprocessor invokes the call
		// The second parameter is the function to call
		// The third is a piece of data to pass (we don't need to pass anything => ~)
		BOOST_PP_REPEAT(AMP_FACTORY_REPEAT, REGISTRATION_BASE, ~)


		#define FACTORY(z, TypeListNumber, data)																								\
		template < class BaseClass, typename Key >																								\
		class Factory<BaseClass, Key, TypeListNumber> : public Singleton< Factory<BaseClass, Key, TypeListNumber> >								\
		{																																		\
		public:																																	\
			friend class Singleton< Factory<BaseClass, Key, TypeListNumber> >;																	\
			bool Register( const Key& key, boost::shared_ptr< RegistrationBase<BaseClass, Key, TypeListNumber> > reg )							\
			{																																	\
					return regMapPtr->insert(typename RegistrationMap::value_type(key, reg)).second;											\
			}																																	\
			        																															\
			bool Register( const Factory<BaseClass, Key, TypeListNumber>& other )																\
			{																																	\
				boost::scoped_ptr<RegistrationMap> tmp(new RegistrationMap(*regMapPtr));														\
				RegistrationMapIterator iter(other.regMapPtr->begin());																			\
				while (iter!=other.regMapPtr->end()) 																							\
				{																																\
					if(!tmp->insert(typename RegistrationMap::value_type(*iter++)).second)														\
						return false;																											\
				}																																\
				std::swap(regMapPtr,tmp);																										\
				return true;																													\
			}																																	\
																																				\
			boost::shared_ptr<BaseClass> create( const Key& id BOOST_PP_COMMA_IF(TypeListNumber) BOOST_PP_ENUM(TypeListNumber,PARAM,~) ) const 	\
			{																																	\
				RegistrationMapIterator iter(regMapPtr->find(id));																				\
				if(iter==regMapPtr->end()) 																										\
					AMP_ERROR("Unregistered creator");																							\
				return iter->second->create( BOOST_PP_ENUM_PARAMS(TypeListNumber,Param) ); 														\
			}																																	\
																																				\
			std::vector<Key> getKeys() const																									\
			{																																	\
				return AMP::voodoo::getKeys(*regMapPtr);																						\
			}																																	\
																																				\
		private:																																\
			typedef typename std::map< Key, boost::shared_ptr< RegistrationBase<BaseClass, Key, TypeListNumber> > > RegistrationMap;			\
			typedef typename RegistrationMap::const_iterator RegistrationMapIterator;															\
			boost::scoped_ptr<RegistrationMap> regMapPtr;																						\
			Factory() : regMapPtr(new RegistrationMap){}																						\
		};

		BOOST_PP_REPEAT(AMP_FACTORY_REPEAT, FACTORY, ~)

		// The macro defintions of the registration Objects

		#define REGISTRATION(z, TypeListNumber, data)																																										\
		template<class BaseClass, class Derived, typename Key>																																								\
		struct Registration<BaseClass, Derived, Key, TypeListNumber> : public RegistrationBase<BaseClass, Key, TypeListNumber> 																								\
		{																																																					\
			Registration(const Key& key)																																													\
			{																																																				\
				try 																																																		\
				{																																																			\
					Factory<BaseClass, Key, TypeListNumber> ::instance().Register( key, boost::shared_ptr< RegistrationBase<BaseClass, Key, TypeListNumber> >(new Registration<BaseClass, Derived, Key, TypeListNumber>) );	\
				}																																																			\
				catch (...){}																																																\
			}																																																				\
			virtual ~Registration(){}																																														\
																																																							\
		private:																																																			\
			boost::shared_ptr<BaseClass> create( BOOST_PP_ENUM(TypeListNumber, PARAM, ~) ) const																															\
			{																																																				\
				return boost::shared_ptr<BaseClass>( new Derived( BOOST_PP_ENUM_PARAMS(TypeListNumber, Param) ) );																											\
			}																																																				\
			Registration(){}																																																\
			Registration( const Registration<BaseClass, Derived, Key, TypeListNumber>& );																																	\
			Registration<BaseClass, Derived, Key, TypeListNumber>& operator=( const Registration<BaseClass, Derived, Key, TypeListNumber>& );																				\
		};

		BOOST_PP_REPEAT(AMP_FACTORY_REPEAT, REGISTRATION, ~)

	} // namespace voodoo
} // namespace amp

#undef PARAM
#undef REGISTRATION
#undef REGISTRATION_BASE
#undef FACTORY

#endif // FACTORY_H




#ifndef DETAILS_H
#define DETAILS_H

/**
 *	\file 		getkeys.h
 *	\brief		Provides a method for collecting keys from a map and shoving into a vector
 *	\author		Nathan Barnett
 */


#include<map>
#include<vector>
#include<iostream>


namespace AMP
{
	namespace voodoo
	{
		/**
		 * \brief	Given a map that is referenced by a key value, takes the keys and
		 *			shoves them into a vector that is returned
		 */
		template<typename key, class T>
		std::vector<key> getKeys( const std::map<key, T>& keyMap ) 
		{
			std::vector<key> keyVector;
			
			// Create an iterator to mark the end, so that we don't have to constantly re-evaluate
			typedef typename std::map<key,T>::const_iterator const_iterator;
			const_iterator	mapEnd( keyMap.end() );

			// Grab the keys from the map and shove into a vector
			for(const_iterator iter(keyMap.begin()); iter!=mapEnd; ++iter) 
			{
				keyVector.push_back(iter->first);
			}

			return keyVector;
		}

	} // namespace VoodooFactory
} // namespace amp


#endif // DETAILS_H


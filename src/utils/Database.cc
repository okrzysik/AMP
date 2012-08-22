//
// File:	$URL: file:///usr/casc/samrai/repository/AMP/tags/v-2-4-4/source/toolbox/database/Database.C $
// Package:	AMP toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2128 $
// Modified:	$LastChangedDate: 2008-04-11 15:29:55 -0700 (Fri, 11 Apr 2008) $
// Description:	An abstract base class for the AMP database objects
//

#include "Database.h"

#ifdef DEBUG_NO_INLINE
#include "Database.I"
#endif

#include "Utilities.h"

namespace AMP {


Database::Database()
{
}


Database::~Database()
{
}



/*  
 * Boolean
 */ 

/*
 * Create a boolean scalar entry in the database with the specified
 * key name.  A scalar entry is an array of one.
 *                                                                      
*/

void Database::putBool(
   const std::string& key, 
   const bool& data)
{
   AMP_ASSERT(!key.empty());
   unsigned char uchar_data = (unsigned char) data;
   
   putBoolArray(key, &uchar_data, 1);
}

/*
 * Create a boolean array entry in the database with the specified       
 * key 
 */

void Database::putBoolArray(
   const std::string& key, 
   const std::vector<unsigned char>& data)
{
   AMP_ASSERT(!key.empty());
   if ( data.size() > 0 ) {
     putBoolArray(key, &(data[0]), data.size());
   } else {
      AMP_ERROR("Database::putBoolArray() error in database "
		 << getName()
		 << "\n    Attempt to put zero-length array with key = "
		 << key << std::endl);
   }
}

/*                                                                      
 * Get boolean scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a boolean type.
 */

bool Database::getBool(const std::string& key)
{
   AMP_ASSERT(!key.empty());
   bool ret_val;
   getBoolArray(key, &ret_val, 1);

   return(ret_val);
}

/*
 * Get boolean scalar entry from the database with the specified key  
 * name. An error message is printed and the program exits if the      
 * specified key does not exist in the database or is not associated   
 * with a boolean type.                                                
 */

bool Database::getBoolWithDefault(
   const std::string& key, 
   const bool& defaultvalue)
{
   AMP_ASSERT(!key.empty());
   if(keyExists(key)) {
      std::vector<unsigned char> local_bool = getBoolArray(key);
      bool *locptr = (bool *)&(local_bool[0]);
      return ( locptr == NULL ? defaultvalue : *locptr);
   } else {
      return defaultvalue;
   }
}

void Database::getBoolArray(
   const std::string& key,
   bool* data, 
   const int nelements)
{
   AMP_ASSERT(!key.empty());
   std::vector<unsigned char> tmp = getBoolArray(key);
   const int tsize = tmp.size();

   if (nelements != tsize) {
      AMP_ERROR("Database::getBoolArray() error in database "
		 << getName()
		 << "\n    Incorrect array size = " << nelements
		 << " given for key = " << key
		 << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }

}

void Database::getScalar(const std::string& key, bool& scalar)
{
   scalar = getBool(key);
}

void Database::putScalar(const std::string& key, const bool scalar)
{
   putBool(key, scalar);
}

void Database::getArray(const std::string& key, std::vector<unsigned char>& array)
{
   array = getBoolArray(key);
}

void Database::putArray(const std::string& key, const std::vector<unsigned char> array)
{
   putBoolArray(key, array);
}

/*
*************************************************************************
*                                                                       *
* Create a box entry in the database with the specified                 *
* key name.  A box entry is an array of one.                            *
*                                                                       *
*************************************************************************
*/

void Database::putDatabaseBox(
   const std::string& key, 
   const DatabaseBox& data)
{
   AMP_ASSERT(!key.empty());
   putDatabaseBoxArray(key, &data, 1);
}

/*
*************************************************************************
*                                                                       *
* Create a box array entry in the database with the specified key name. *
*                                                                       *
*************************************************************************
*/

void Database::putDatabaseBoxArray(
   const std::string& key, 
   const std::vector<DatabaseBox>& data)
{
   AMP_ASSERT(!key.empty());
   if ( data.size() > 0 ) {
     putDatabaseBoxArray(key, &(data[0]), data.size());
   } else {
      AMP_ERROR("Database::putDatabaseBoxArray() error in database "
		 << getName()
		 << "\n    Attempt to put zero-length array with key = "
		 << key << std::endl);
   }
}

/*
************************************************************************
*                                                                      *
* Get box scalar entry from the database with the specified key        *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a box type.                                                     *
*                                                                      *
************************************************************************
*/

DatabaseBox Database::getDatabaseBox(const std::string& key)
{
   AMP_ASSERT(!key.empty());
   DatabaseBox ret_val;
   getDatabaseBoxArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Get box scalar entry from the database with the specified key        *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a box type.                                                     *
*                                                                      *
************************************************************************
*/

DatabaseBox Database::getDatabaseBoxWithDefault(
   const std::string& key,
   const DatabaseBox& defaultvalue)
{
   AMP_ASSERT(!key.empty());
   if(keyExists(key)) {
      std::vector<DatabaseBox> local_box = getDatabaseBoxArray(key);
      DatabaseBox *locptr = &(local_box[0]);
      return ( locptr == NULL ? defaultvalue : *locptr);
   } else {
      return defaultvalue;
   }

}


void Database::getDatabaseBoxArray(
   const std::string& key,
   DatabaseBox* data,
   const int nelements)
{
   AMP_ASSERT(!key.empty());
   std::vector<DatabaseBox> tmp = getDatabaseBoxArray(key);
   const int tsize = tmp.size();

   if (nelements != tsize) {
      AMP_ERROR("Database::getDatabaseBoxArray() error in database "
		 << getName()
		 << "\n    Incorrect array size = " << nelements
		 << " given for key = " << key
		 << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*  
 * Char
 */ 

/*
*************************************************************************
*                                                                       *
* Create a char scalar entry in the database with the specified         *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void Database::putChar(
   const std::string& key, 
   const char& data)
{
   AMP_ASSERT(!key.empty());
   putCharArray(key, &data, 1);

}

/*
*************************************************************************
*                                                                       *
* Create a char array entry in the database with the specified          *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void Database::putCharArray(
   const std::string& key, 
   const std::vector<char>& data)
{
   AMP_ASSERT(!key.empty());
   if ( data.size() > 0 ) {
     putCharArray(key, &(data[0]), data.size());
   } else { 
      AMP_ERROR("Database::putCharArray() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
************************************************************************
*                                                                      *
* Get char scalar entry from the database with the specified key       *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a char type.                                                    *
*                                                                      *
************************************************************************
*/

char Database::getChar(const std::string& key)
{
   AMP_ASSERT(!key.empty());
   char ret_val;
   getCharArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Get char scalar entry from the database with the specified key       *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a char type.                                                    *
*                                                                      *
************************************************************************
*/

char Database::getCharWithDefault(
   const std::string& key,
   const char& defaultvalue)
{
   AMP_ASSERT(!key.empty());
   if(keyExists(key)) {
      std::vector<char> local_char = getCharArray(key);
      char *locptr = &(local_char[0]);
      return ( locptr == NULL ? defaultvalue : *locptr);
   } else {
      return defaultvalue;
   }

}


void Database::getCharArray(
   const std::string& key,
   char* data,
   const int nelements)
{
   AMP_ASSERT(!key.empty());
   std::vector<char> tmp = getCharArray(key);
   const int tsize = tmp.size();

   if (nelements != tsize) {
      AMP_ERROR("Database::getCharArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

void Database::getScalar(const std::string& key, char& scalar)
{
   scalar = getChar(key);
}

void Database::putScalar(const std::string& key, const char scalar)
{
   putChar(key, scalar);
}

void Database::getArray(const std::string& key, std::vector<char>& array)
{
   array = getCharArray(key);
}

void Database::putArray(const std::string& key, const std::vector<char> array)
{
   putCharArray(key, array);
}

/*  
 * Complex
 */ 

/*
*************************************************************************
*                                                                       *
* Create a complex scalar entry in the database with the specified      *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void Database::putComplex(
   const std::string& key,
   const std::complex<double>& data)
{
   AMP_ASSERT(!key.empty());
   putComplexArray(key, &data, 1);
}

/*
*************************************************************************
*                                                                       *
* Create a complex array entry in the database with the specified       *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void Database::putComplexArray(
   const std::string& key,
   const std::vector<std::complex<double> >& data)
{
   AMP_ASSERT(!key.empty());
   if ( data.size() > 0 ) {
     putComplexArray(key, &(data[0]), data.size());
   } else {
      AMP_ERROR("Database::putComplexArray() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
************************************************************************
*                                                                      *
* Get complex scalar entry from the database with the specified key    *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a complex type.                                                 *
*                                                                      *
************************************************************************
*/

std::complex<double> Database::getComplex(const std::string& key)
{
   AMP_ASSERT(!key.empty());
   std::complex<double> ret_val;
   getComplexArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Get complex scalar entry from the database with the specified key    *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a complex type.                                                 *
*                                                                      *
************************************************************************
*/

std::complex<double> Database::getComplexWithDefault(
   const std::string& key,
   const std::complex<double>& defaultvalue)
{
   AMP_ASSERT(!key.empty());
   if(keyExists(key)) {
      std::vector<std::complex<double> > local_complex = getComplexArray(key);
      std::complex<double> *locptr = &(local_complex[0]);
      return ( locptr == NULL ? defaultvalue : *locptr);   
   } else {
      return defaultvalue;
   }

}

void Database::getComplexArray(
   const std::string& key,
   std::complex<double>* data,
   const int nelements)
{
   AMP_ASSERT(!key.empty());
   std::vector<std::complex<double> > tmp = getComplexArray(key);
   const int tsize = tmp.size();

   if (nelements != tsize) {
      AMP_ERROR("Database::getComplexArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}




void Database::getScalar(const std::string& key, std::complex<double>& scalar)
{
   scalar = getComplex(key);
}

void Database::putScalar(const std::string& key, const std::complex<double> scalar)
{
   putComplex(key, scalar);
}

void Database::getArray(const std::string& key, std::vector<std::complex<double> >& array)
{
   array = getComplexArray(key);
}

void Database::putArray(const std::string& key, const std::vector<std::complex<double> > array)
{
   putComplexArray(key, array);
}

/*  
 * Float
 */ 

/*
*************************************************************************
*                                                                       *
* Create a float scalar entry in the database with the specified        *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void Database::putFloat(
   const std::string& key, 
   const float& data)
{
   AMP_ASSERT(!key.empty());
   putFloatArray(key, &data, 1);
}

/*
*************************************************************************
*                                                                       *
* Create a float array entry in the database with the specified         *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void Database::putFloatArray(
   const std::string& key,
   const std::vector<float>& data)
{
   AMP_ASSERT(!key.empty());
   if ( data.size() > 0 ) {
     putFloatArray(key, &(data[0]), data.size());
   } else {
      AMP_ERROR("Database::putFloatArray() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
  
}

/*
************************************************************************
*                                                                      *
* Get float scalar entry from the database with the specified key      *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a float type.                                                   *
*                                                                      *
************************************************************************
*/

float Database::getFloat(const std::string& key)
{
   AMP_ASSERT(!key.empty());
   float ret_val;
   getFloatArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Get float scalar entry from the database with the specified key      *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a float type.                                                   *
*                                                                      *
************************************************************************
*/

float Database::getFloatWithDefault(
   const std::string& key, 
   const float& defaultvalue)
{
   AMP_ASSERT(!key.empty());
   if(keyExists(key)) {
      std::vector<float> local_float = getFloatArray(key);
      float *locptr = &(local_float[0]);
      return ( locptr == NULL ? defaultvalue : *locptr);
   } else {
      return defaultvalue;
   }

}


void Database::getFloatArray(
   const std::string& key,
   float* data,
   const int nelements)
{
   AMP_ASSERT(!key.empty());
   std::vector<float> tmp = getFloatArray(key);
   const int tsize = tmp.size();

   if (nelements != tsize) {
      AMP_ERROR("Database::getFloatArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

void Database::getScalar(const std::string& key, float& scalar)
{
   scalar = getFloat(key);
}

void Database::putScalar(const std::string& key, const float scalar)
{
   putFloat(key, scalar);
}

void Database::getArray(const std::string& key, std::vector<float>& array)
{
   array = getFloatArray(key);
}

void Database::putArray(const std::string& key, const std::vector<float> array)
{
   putFloatArray(key, array);
}

/*  
 * Double
 */ 


/*
*************************************************************************
*                                                                       *
* Create a double scalar entry in the database with the specified       *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void Database::putDouble(
   const std::string& key, 
   const double& data)
{
   AMP_ASSERT(!key.empty());
   putDoubleArray(key, &data, 1);
}

/*
*************************************************************************
*                                                                       *
* Create a double array entry in the database with the specified        *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void Database::putDoubleArray(
   const std::string& key,
   const std::vector<double>& data)
{
   AMP_ASSERT(!key.empty());
   if ( data.size() > 0 ) {
     putDoubleArray(key, &(data[0]), data.size());
   } else {
      AMP_ERROR("Database::putDoubleArray() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
************************************************************************
*                                                                      *
* Get double scalar entry from the database with the specified key     *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a double type.                                                  *
*                                                                      *
************************************************************************
*/

double Database::getDouble(const std::string& key)
{
   AMP_ASSERT(!key.empty());
   double ret_val;
   getDoubleArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Get double scalar entry from the database with the specified key     *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a double type.                                                  *
*                                                                      *
************************************************************************
*/

double Database::getDoubleWithDefault(
   const std::string& key, 
   const double& defaultvalue)
{
   AMP_ASSERT(!key.empty());
   if(keyExists(key)) {
      std::vector<double> local_double = getDoubleArray(key);
      double *locptr = &(local_double[0]);
      return ( locptr == NULL ? defaultvalue : *locptr); 
   } else {
      return defaultvalue;
   }
}

void Database::getDoubleArray(
   const std::string& key, 
   double* data, 
   const int nelements)
{
   AMP_ASSERT(!key.empty());
   std::vector<double> tmp = getDoubleArray(key);
   const int tsize = tmp.size();

   if (nelements != tsize) {
      AMP_ERROR("Database::getDoubleArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}



void Database::getScalar(const std::string& key, double& scalar)
{
   scalar = getDouble(key);
}

void Database::putScalar(const std::string& key, const double scalar)
{
   putDouble(key, scalar);
}

void Database::getArray(const std::string& key, std::vector<double>& array)
{
   array = getDoubleArray(key);
}

void Database::putArray(const std::string& key, const std::vector<double> array)
{
   putDoubleArray(key, array);
}

/*  
 * Integer
 */ 

/*
*************************************************************************
*                                                                       *
* Create a integer scalar entry in the database with the specified      *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void Database::putInteger(
   const std::string& key,
   const int& data)
{
   AMP_ASSERT(!key.empty());
   putIntegerArray(key, &data, 1);
}

/*
*************************************************************************
*                                                                       *
* Create an integer array entry in the database with the specified      *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void Database::putIntegerArray(
   const std::string& key, 
   const std::vector<int>& data)
{
   AMP_ASSERT(!key.empty());
   if ( data.size() > 0 ) {
     putIntegerArray(key, &(data[0]), data.size());
   } else {
      AMP_ERROR("Database::putIntegerArray() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}


/*
************************************************************************
*                                                                      *
* Get integer scalar entry from the database with the specified key    *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a integer type.                                                 *
*                                                                      *
************************************************************************
*/

int Database::getInteger(const std::string& key)
{
   AMP_ASSERT(!key.empty());
   int ret_val;
   getIntegerArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Get integer scalar entry from the database with the specified key    *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a integer type.                                                 *
*                                                                      *
************************************************************************
*/

int Database::getIntegerWithDefault(
   const std::string& key, 
   const int& defaultvalue)
{
   AMP_ASSERT(!key.empty());
   if(keyExists(key)) {
      std::vector<int> local_int = getIntegerArray(key);
      int *locptr = &(local_int[0]);
      return ( locptr == NULL ? defaultvalue : *locptr);
   } else {
      return defaultvalue;
   }

}

void Database::getIntegerArray(
   const std::string& key, 
   int* data, 
   const int nelements)
{
   AMP_ASSERT(!key.empty());
   std::vector<int> tmp = getIntegerArray(key);
   const int tsize = tmp.size();

   if (nelements != tsize) {
      AMP_ERROR("Database::getIntegerArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

void Database::getScalar(const std::string& key, int& scalar)
{
   scalar = getInteger(key);
}

void Database::putScalar(const std::string& key, const int scalar)
{
   putInteger(key, scalar);
}

void Database::getArray(const std::string& key, std::vector<int>& array)
{
   array = getIntegerArray(key);
}

void Database::putArray(const std::string& key, const std::vector<int> array)
{
   putIntegerArray(key, array);
}

/*
 * String
 */ 

/*
*************************************************************************
*                                                                       *
* Create a string scalar entry in the database with the specified       *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void Database::putString(
   const std::string& key, 
   const std::string& data)
{
   AMP_ASSERT(!key.empty());
   putStringArray(key, &data, 1);
}

/*
*************************************************************************
*                                                                       *
* Create a string array entry in the database with the specified        *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void Database::putStringArray(
   const std::string& key, 
   const std::vector<std::string>& data)
{
   AMP_ASSERT(!key.empty());
   if ( data.size() > 0 ) {
     putStringArray(key, &(data[0]), data.size());
   } else {
      AMP_ERROR("Database::putStringArray() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
************************************************************************
*                                                                      *
* Get string scalar entry from the database with the specified key     *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a string type.                                                  *
*                                                                      *
************************************************************************
*/

std::string Database::getString(const std::string& key)
{
   AMP_ASSERT(!key.empty());
   std::string ret_val;
   getStringArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Get string scalar entry from the database with the specified key     *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a string type.                                                  *
*                                                                      *
************************************************************************
*/

std::string Database::getStringWithDefault(
   const std::string& key, 
   const std::string& defaultvalue)
{
   AMP_ASSERT(!key.empty());
   if(keyExists(key)) {
      std::vector<std::string> local_string = getStringArray(key);
      std::string *locptr = &(local_string[0]);
      return ( locptr == NULL ? defaultvalue : *locptr);
   } else {
      return defaultvalue;
   }

}

void Database::getStringArray(
   const std::string& key, 
   std::string* data, 
   const int nelements)
{
   AMP_ASSERT(!key.empty());
   std::vector<std::string> tmp = getStringArray(key);
   const int tsize = tmp.size();

   if (nelements != tsize) {
      AMP_ERROR("Database::getStringArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }

}

}


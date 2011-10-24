#ifndef __included_AMP_Counter
#define __included_AMP_Counter

#include <string>
#include <map>
#include <vector>
#include "utils/AMP_MPI.h"

#include "utils/Utilities.h"

namespace AMP
{

  class  Counter
  {
    private:
      static std::map<std::string , size_t>   d_Data;

    public:
      static void initCounter ();

      static void   increment ( const std::string &key , size_t amt = 1 );
      static size_t total ( const std::string &key , AMP_MPI comm = AMP_MPI(AMP_COMM_WORLD) );
  };

#ifdef  _USE_COUNTER
#define INCREMENT_COUNT(x) Counter::increment(x)
#define ADD_COUNT(x,y) Counter::increment(x,y)
  inline 
  void Counter::increment ( const std::string &key , size_t amt )
  {
    AMP_INSIST ( d_Data.find ( key ) != d_Data.end() , "Uninitialized key" );
    d_Data[key] += amt;
  }

  inline
  size_t Counter::total ( const std::string &key , AMP_MPI comm )
  {
    AMP_INSIST ( d_Data.find ( key ) != d_Data.end() , "Uninitialized key" );
    unsigned long retVal = comm.sumReduce((unsigned long) d_Data[key]);
    return (size_t) retVal;
  }
#else
#define INCREMENT_COUNT(x)
#define ADD_COUNT(x,y)

  inline 
  void Counter::increment ( const std::string & , size_t )
  {
  }

  inline
  size_t Counter::total ( const std::string & , AMP_MPI )
  {
    return 0;
  }
#endif

}


#endif

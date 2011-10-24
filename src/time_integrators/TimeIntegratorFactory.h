/*
Copyright 2005, The Regents of the University 
of California. This software was produced under
a U.S. Government contract (W-7405-ENG-36) 
by Los Alamos National Laboratory, which is
operated by the University of California for the
U.S. Department of Energy. The U.S.
Government is licensed to use, reproduce, and
distribute this software. Permission is granted
to the public to copy and use this software
without charge, provided that this Notice and
any statement of authorship are reproduced on
all copies. Neither the Government nor the
University makes any warranty, express or
implied, or assumes any liability or
responsibility for the use of this software.
*/

#ifndef included_TimeIntegratorFactory
#define included_TimeIntegratorFactory

#ifndef included_AMP_config

#endif

#ifndef included_TimeIntegratorParameters
#include "TimeIntegratorParameters.h"
#endif

#ifndef included_TimeIntegrator
#include "TimeIntegrator.h"
#endif

#include <memory>


namespace AMP{
namespace TimeIntegrator{

/**\class TimeIntegratorFactory
 * 
 * TimeIntegratorFactory is a factory class that creates specific multilevel
 * solver classes.  These are used to provide methods that operate on
 * a SAMR hierarchy
 */
class TimeIntegratorFactory
{
public:
   /**
    * Constructor.
    */
   TimeIntegratorFactory();

   /**
    * Destructor.
    */
   ~TimeIntegratorFactory();

   /**
    * Factory method for generating multilevel solvers with characteristics
    * specified by parameters.
    */
   boost::shared_ptr<TimeIntegrator> createTimeIntegrator( boost::shared_ptr<TimeIntegratorParameters> timeIntegratorParameters);

protected:
private:  
};

}
}

#endif //included_TimeIntegratorFactory


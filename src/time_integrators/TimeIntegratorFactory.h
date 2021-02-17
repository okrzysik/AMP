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

#ifndef included_AMP_TimeIntegratorFactory
#define included_AMP_TimeIntegratorFactory

#include "TimeIntegrator.h"
#include "TimeIntegratorParameters.h"


namespace AMP {
namespace TimeIntegrator {

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
    virtual ~TimeIntegratorFactory();

    /**
     * Factory method for generating multilevel solvers with characteristics
     * specified by parameters.
     */
    std::shared_ptr<TimeIntegrator>
    createTimeIntegrator( std::shared_ptr<TimeIntegratorParameters> timeIntegratorParameters );

protected:
private:
};
} // namespace TimeIntegrator
} // namespace AMP

#endif

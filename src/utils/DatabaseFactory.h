//
// File:	$URL$
// Package:	AMP toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2086 $
// Modified:	$LastChangedDate: 2008-03-28 15:10:07 -0700 (Fri, 28 Mar 2008) $
// Description:	An abstract base class for a DatabaseFactory
//

#ifndef included_DatabaseFactory
#define included_DatabaseFactory


#include "Database.h"

namespace AMP {


/**
 * @brief Abstract base class factory used to build Database objects.
 *
 * Used to build database objects.  For example, RestartManager
 * may use a DatabaseFactory to build databases when creating
 * a restart file.
 */
class DatabaseFactory {
public:
    /*
     * Build a new Database instance.
     */
    virtual AMP::shared_ptr<Database> allocate( const std::string &name ) = 0;
    virtual ~DatabaseFactory() {}
};
}

#endif

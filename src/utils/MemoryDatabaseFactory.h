//
// File:	$URL$
// Package:	AMP toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2122 $
// Modified:	$LastChangedDate: 2008-04-08 15:37:28 -0700 (Tue, 08 Apr 2008) $
// Description:	A factory for building MemoryDatabases
//

#ifndef included_MemoryDatabaseFactory
#define included_MemoryDatabaseFactory


#include "DatabaseFactory.h"

namespace AMP {


/**
 * @brief MemoryDatabase factory.
 *
 * Builds a new MemoryDatabase.
 */
class MemoryDatabaseFactory : public DatabaseFactory
{
    /**
     * Build a new Database object.
     */
    virtual AMP::shared_ptr<Database> allocate( const std::string &name ) override;
};
} // namespace AMP

#endif

//
// File:	$URL$
// Package:	AMP toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2122 $
// Modified:	$LastChangedDate: 2008-04-08 15:37:28 -0700 (Tue, 08 Apr 2008) $
// Description:	A factory for building MemoryDatabases
//

#include "MemoryDatabaseFactory.h"
#include "MemoryDatabase.h"

namespace AMP {
  

/**
 * Build a new MemoryDatabase object.
 */
 boost::shared_ptr<Database> MemoryDatabaseFactory::allocate(const std::string& name) {
  boost::shared_ptr<MemoryDatabase> database(new MemoryDatabase(name));
   return database;
}


}

//
// File:	$URL:
// file:///usr/casc/samrai/repository/AMP/trunk/source/toolbox/inputdb/InputDatabase.h
// $
// Package:	AMP toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2662 $
// Modified:	$LastChangedDate: 2008-11-24 16:44:41 -0800 (Mon, 24 Nov 2008) $
// Description:	An input database structure that stores (key,value) pairs
//
// NOTE: This is a modification of the MemoryDatabase class from SAMRAI
//       We have simply used it with modifications


#ifndef included_AMP_InputDatabase
#define included_AMP_InputDatabase


#include "MemoryDatabase.h"

namespace AMP {


/**
 * @brief Class InputDatabase stores (key,value) pairs in a hierarchical
 * database.
 *
 * This is just another name for the MemoryDatabase. @see MemoryDatabase
 *
 * It is normally filled with data using a Parser (@see
 * Parser) and used to pass user supplied input from input files
 * to constructors for problem setup.
 *
 */
typedef AMP::MemoryDatabase InputDatabase;
} // namespace AMP
#endif

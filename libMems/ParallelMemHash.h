/*******************************************************************************
 * $Id: ParallelMemHash.h,v 1.4 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifndef _ParallelMemHash_h_
#define _ParallelMemHash_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MemHash.h"

namespace mems {

/**
 * The default amount of sequence to preprocess for matches before spawning threads
 * to continue searching.  Preprocessing is done to avoid heavy lock contention among
 * threads.  This should be a number between 0 and 1. 
 */
const float DEFAULT_PREPROCESS_PERCENT = .01; 

/**
 * Implements a parallelized version of exact match extraction
 */
class ParallelMemHash : public MemHash
{
public:
	ParallelMemHash();
	~ParallelMemHash();
	ParallelMemHash(const ParallelMemHash& mh);
	virtual ParallelMemHash* Clone() const;
	virtual void Clear();
	/**
	 * Search for matches using the specified number of threads.
	 * @param thread_count The number of threads to spawn.  
	 * Generally this should be set to the number of processors in the system.
	 */
	virtual void CreateMems(uint32 thread_count);

protected:
	float m_preprocess_percent;
};

}

#endif //_ParallelMemHash_h_

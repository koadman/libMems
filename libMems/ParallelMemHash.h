/*******************************************************************************
 * $Id: ParallelMemHash.h,v 1.23 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _ParallelMemHash_h_
#define _ParallelMemHash_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef _OPENMP

#include <omp.h>
#include "libMems/MemHash.h"

namespace mems {


/**
 * ParallelMemHash implements an algorithm for finding exact matches of a certain minimal
 * length in several sequences.
 */
class ParallelMemHash : public MemHash{

   

public:
	ParallelMemHash();
	ParallelMemHash(const ParallelMemHash& mh);
	ParallelMemHash& operator=( const ParallelMemHash& mh );
	virtual ParallelMemHash* Clone() const;
	
	/**
	 * Finds (in parallel) all matches in the sequences contained by "match_list"
	 * The resulting list of matches is stored within "match_list"
	 */
	virtual void FindMatches( MatchList& match_list );


	/**
	 * Sets the size of the hash table to new_table_size.
	 * @param new_table_size The new hash table size
	 */
	virtual void SetTableSize(uint32 new_table_size);	

protected:
	virtual MatchHashEntry* AddHashEntry(MatchHashEntry& mhe);

	std::vector< omp_lock_t > omp_locks;
};



}

#endif // _OPENMP


}

#endif //_ParallelMemHash_h_

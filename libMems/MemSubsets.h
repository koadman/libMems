/*******************************************************************************
 * $Id: MemSubsets.h,v 1.8 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifndef _MemSubsets_h_
#define _MemSubsets_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/Match.h"
#include <list>
#include <vector>
#include <map>

namespace mems {

/**
 * This virtual class provides a general framework for locating
 * evenly overlapping mems.
 */
class MemSubsets
{
public:
	MemSubsets();
	MemSubsets( const MemSubsets& ms );
//	~MemSubsets(){};
	/**
	 * Call this function to initiate action on a mem list.
	 * This function locates each evenly overlapping mem and
	 * calls HandleOverlap to perform some function on the overlap.
	 */
	virtual void Subsets( std::vector<Match*>& mem_list );
	
	/**
	 * Eliminates all area in subset matches which is covered by a match of
	 * higher multiplicity
	 */
	virtual void EliminateLinkedInclusions( std::vector<Match*>& mem_list );

protected:
	virtual void HandleOverlap( Match* m1, Match* m2 ) = 0;

	static uint32 quadratic_li(uint32 listI){return (listI*(listI+1))/2;}
	
	// set to true to enable debugging messages
	boolean debug_algorithm;

	std::list<Match*>::iterator mhe, mem_iter, prev_iter;
	//define an array of lists
	//first index is the number of sequences in the match
	//second index is the sequence that the list is sorted on.
	std::list<Match*>** mhe_list;

	uint32 list_count;
	Match* new_mie;

	//keep a vector of iterators with the current positions to check in
	//each level of sorted match lists
	std::vector<std::list<Match*>::iterator>* start_iters;
	std::vector<std::list<Match*>::iterator>* cur_starts;

	int32 listI;

	std::list<Match*>::iterator cur_start_iter, next_start_iter, matching_iter;
	boolean erased_start;

};

class SubsetLinker : public MemSubsets
{
protected:
	virtual void HandleOverlap( Match* m1, Match* m2 );
};

class SubsetInclusionRemover : public MemSubsets
{
protected:
	virtual void HandleOverlap( Match* m1, Match* m2 );
};

}

#endif // _MemSubsets_h_


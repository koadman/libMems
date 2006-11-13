#ifndef _MemSubsets_h_
#define _MemSubsets_h_

#include "MemHashEntry.h"
#include <list>
#include <vector>
#include <map>

/**
 * This virtual class provides a general framework for locating
 * evenly overlapping mems.
 */
class MemSubsets
{
public:
	MemSubsets(){};
	MemSubsets( const MemSubsets& ms );
	~MemSubsets(){};
	/**
	 * Call this function to initiate action on a mem list.
	 * This function locates each evenly overlapping mem and
	 * calls HandleOverlap to perform some function on the overlap.
	 */
	virtual void Subsets( list<MemHashEntry*>& mem_list );
	
	/**
	 * Eliminates all area in subset matches which is covered by a match of
	 * higher multiplicity
	 */
	virtual void EliminateLinkedInclusions( list<MemHashEntry*>& mem_list );

	/**
	 * Eliminates all area in subset matches which is covered by a match of
	 * higher multiplicity
	 */
	void EliminateLinkedInclusions( list<Match*>& mem_list );

	void UnlinkMatch( Match* match, map<MatchID_t, Match*> match_id_map );
	void UpdateSupersets( Match* match, map<MatchID_t, Match*> match_id_map );

protected:
	virtual void HandleOverlap( MemHashEntry* m1, MemHashEntry* m2 ) = 0;

	static uint32 quadratic_li(uint32 listI){return (listI*(listI+1))/2;}
	
	// set to true to enable debugging messages
	boolean debug_algorithm;

	list<MemHashEntry*>::iterator mhe, mem_iter, prev_iter;
	//define an array of lists
	//first index is the number of sequences in the match
	//second index is the sequence that the list is sorted on.
	list<MemHashEntry*>** mhe_list;

	uint32 list_count;
	MemHashEntry* new_mie;

	//keep a vector of iterators with the current positions to check in
	//each level of sorted match lists
	vector<list<MemHashEntry*>::iterator>* start_iters;
	vector<list<MemHashEntry*>::iterator>* cur_starts;

	int32 listI;

	list<MemHashEntry*>::iterator cur_start_iter, next_start_iter, matching_iter;
	boolean erased_start;
};

class SubsetLinker : public MemSubsets
{
protected:
	virtual void HandleOverlap( MemHashEntry* m1, MemHashEntry* m2 );
};

class SubsetInclusionRemover : public MemSubsets
{
protected:
	virtual void HandleOverlap( MemHashEntry* m1, MemHashEntry* m2 );
};


#endif // _MemSubsets_h_


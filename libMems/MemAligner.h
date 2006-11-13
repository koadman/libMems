#ifndef _MemAlignerz_h_
#define _MemAlignerz_h_
#include "gn/gnClone.h"
#include <gn/gnSequence.h>
#include <set>
#include <map>
#include <list>

#include "MimHashEntry.h"

/** This class is a toy for experimenting with different alignment ideas
  * Do not use this for real alignments!  It doesn't work.
  */
class MemAligner : public gnClone{

public:
	MemAligner(uint32 seqcount);
	MemAligner( const MemAligner& ma );
	MemAligner* Clone() const { return new MemAligner( *this ); }
	void CreateAlignment(list<MimHashEntry*>* mem_list);
	void CreateOrderedSet();

protected:
	void CreateSet(MimHashEntry* current_mem);
	void ExamineUnderlaps(MimHashEntry* parent);
	void ExamineOverlaps(MimHashEntry* child);

	float64 DistanceScore(MimHashEntry* m1, MimHashEntry* m2);
	float64 DistanceScore(MemHashEntry* m1, MemHashEntry* m2);

	uint32 seq_count;

	//this set of mem pointers is sorted on size
	//so that its easy to always pull the biggest one out.
	multiset<MimHashEntry*, MimLengthCompare> mem_size_set;

	vector<uint32> ordered_set_starts;
	vector<MimHashEntry*> ordered_supersets;
	vector<MimHashEntry*> core_sets;
	//make a set of orthologous alignments
	set<MimHashEntry*> alignment_group;
	vector<uint32> group_mem_count;
	vector<uint32> group_coverage_count;

	//this is an array of vectors
	//each array entry is 
	vector<MimHashEntry*>* current_set;

	uint32 current_alignment;

	//maybe derive gnSequence to create a gap class which
	//doesn't have to store a gap character for every gap.
	//or better yet, create a gnGapSpec class which implements
	//gnContigSpec.
	vector<gnSequence> alignment;
};

#endif  //_MemAligner_h_

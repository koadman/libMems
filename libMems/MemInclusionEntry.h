#ifndef _MemInclusionEntry_h_
#define _MemInclusionEntry_h_
#include "MemHashEntry.h"

class MemInclusionEntry : public MemHashEntry{
public:
	MemInclusionEntry();
	MemInclusionEntry(const MemInclusionEntry& mie);
	MemInclusionEntry& operator=(const MemInclusionEntry& mie);
	MemInclusionEntry(const MemHashEntry& mie);
	MemInclusionEntry& operator=(const MemHashEntry& mie);
	MemInclusionEntry* Clone() const;

	set<MemInclusionEntry*> m_overlaps;
	set<MemInclusionEntry*> m_underlaps;
};

class MieCompare {
public:
	boolean operator()(const MemInclusionEntry* a, const MemInclusionEntry* b){
		return MemHashEntry::start_lessthan_ptr((const MemHashEntry*)a, (const MemHashEntry*)b);
	}
};

class MieLengthCompare {
public:
	/**
	 * Compares two mems to determine which is the longer of the two.
	 * The first criteria is the number of matching sequences
	 * The second criteria is the length
	 * The third criteria for comparison is the start position.
	 */
	boolean operator()(const MemInclusionEntry* a, const MemInclusionEntry* b){
		boolean retval = false;
		if(a->Multiplicity() == b->Multiplicity()){
			retval = a->Length() > b->Length();
			if(a->Length() == b->Length()){
				if(a->SeqCount() == b->SeqCount()){
					for(uint32 seqI = 0; seqI < a->SeqCount(); seqI++){
						if(a->Start(seqI) == b->Start(seqI))
							continue;
						return a->Start(seqI) > b->Start(seqI);
					}
				}
				return a->SeqCount() > b->SeqCount();
			}
			return retval;
		}
		return a->Multiplicity() > b->Multiplicity();
	}
};

#endif // _MemInclusionEntry_h_

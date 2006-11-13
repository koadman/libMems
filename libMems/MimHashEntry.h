#ifndef _MimHashEntry_h_
#define _MimHashEntry_h_

#include "gn/gnClone.h"
#include "MemHashEntry.h"
#include "gn/gnSequence.h"
#include <set>

class MimHashEntry : public gnClone{
public:
	MimHashEntry();
	~MimHashEntry(){};
	MimHashEntry( const MimHashEntry& mh);
	MimHashEntry* Clone() const;
	void AddMem(MemHashEntry* mhe);
	void Merge(MimHashEntry& mhe);
	void DeleteMem(MemHashEntry* mhe);

	int64 start(uint32 startI){return m_start[startI];}
	int64 end(uint32 endI){return m_end[endI];}
	//compare the start of this mim to the end of the given mem
	int start_compare(MemHashEntry& mhe, uint32 tolerance);
	//compare the end of this mim to the start of the given mem
	int end_compare(MemHashEntry& mhe, uint32 tolerance);
	//return the number of mems which makes up this mim
	uint32 MemCount(){return mem_list.size();}
	gnSeqI BaseLength() const;
	
	void PrintMismatches(ostream& os, vector<string>& seq_strings, uint32 context);
	
	/**
	 * Writes the location of this mem to the specified output stream (e.g. cout).
	 */
	friend std::ostream& operator<<(std::ostream& os, const MimHashEntry& mhe); //write to source.
protected:
	set<MemHashEntry*, MheCompare> mem_list;
	vector<int64> m_start;
	vector<int64> m_end;
	gnSeqI m_baseLength;
};

std::ostream& operator<<(std::ostream& os, const MimHashEntry& mhe); //write to source.

class MimLengthCompare {
public:
	/**
	 * Compares two mems to determine which is the longer of the two.
	 */
	boolean operator()(const MimHashEntry* a, const MimHashEntry* b) const{
		return a->BaseLength() < b->BaseLength();
	}
};


#endif
//  _MimHashEntry_h_

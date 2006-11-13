#ifndef _MimHash_h_
#define _MimHash_h_

#include "MemHash.h"
#include "MimHashEntry.h"
#include "MemInclusionEntry.h"
#include <set>

/**
 * The MemCollapseEntry class is used to represent a doubly linked list whose forward links are sorted
 * differently than the previous links.
 */
class MemCollapseEntry{
public:
	MemCollapseEntry();
	MemCollapseEntry(MemHashEntry* mie);
	
	MemHashEntry* m_mie;
	MemCollapseEntry* next_start;
	MemCollapseEntry* prev_start;
	MemCollapseEntry* next_end;
	MemCollapseEntry* prev_end;
	void DeleteSelf();
	static boolean start_lessthan_ptr(const MemCollapseEntry* a, const MemCollapseEntry* b);
	static boolean end_lessthan_ptr(const MemCollapseEntry* a, const MemCollapseEntry* b);
};


class MimHash : public MemHash{
public:
	MimHash(){}
	~MimHash(){}
	MimHash(const MimHash& mh);
	virtual MimHash* Clone() const;
	virtual void Clear(){m_mim_list.clear();}
	virtual void ResolveMismatches( uint32 tolerance = 1);
	virtual void GetMimList( vector<MimHashEntry>& mim_list );
	virtual void LocateSubsetOverlaps();
	virtual void PrintMismatches(ostream& os, vector<gnSequence*>& seqs);
	virtual void EliminateInclusions();

protected:	
	vector<MimHashEntry*> m_mim_list;
};

inline
boolean MemCollapseEntry::start_lessthan_ptr(const MemCollapseEntry* a, const MemCollapseEntry* b){
	return MemHashEntry::start_lessthan_ptr(a->m_mie, b->m_mie);
};

inline
boolean MemCollapseEntry::end_lessthan_ptr(const MemCollapseEntry* a, const MemCollapseEntry* b){
	return MemHashEntry::end_lessthan(*(a->m_mie), *(b->m_mie));
};

#endif
// _MimHash_h_

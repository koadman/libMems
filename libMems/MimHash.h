#ifndef _MimHash_h_
#define _MimHash_h_

#include "MemHash.h"
#include "MimHashEntry.h"

class MemCollapseEntry{
public:
	MemCollapseEntry();
	MemCollapseEntry(MemHashEntry& mhe);
	
	MemHashEntry m_mem;
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
	MimHash(){m_mim_list.clear();}
	MimHash(uint32 my_table_size){table_size = my_table_size;}
	~MimHash(){}
	MimHash(const MimHash& mh);
	virtual MimHash* Clone() const;
	virtual boolean ResolveMismatches( uint32 tolerance = 1);
	virtual boolean GetMimList( vector<MimHashEntry>& mim_list );
	
private:
	vector<MimHashEntry> m_mim_list;
};

inline
boolean MemCollapseEntry::start_lessthan_ptr(const MemCollapseEntry* a, const MemCollapseEntry* b){
	return MemHashEntry::start_lessthan_ptr(&a->m_mem, &b->m_mem);
};

inline
boolean MemCollapseEntry::end_lessthan_ptr(const MemCollapseEntry* a, const MemCollapseEntry* b){
	return MemHashEntry::end_lessthan(a->m_mem, b->m_mem);
};

#endif
// _MimHash_h_

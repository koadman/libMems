#include "MimHashEntry.h"
#include "MemHashEntry.h"

MimHashEntry::MimHashEntry( const MimHashEntry& mh){
	mem_list = mh.mem_list;
	m_start = mh.m_start;
	m_end = mh.m_end;
	m_baseLength = mh.m_baseLength;
	m_gapLength = mh.m_gapLength;
}

MimHashEntry* MimHashEntry::Clone() const{
	return new MimHashEntry(*this);
}

void MimHashEntry::AddMem(MemHashEntry& mhe){
	//add in order
	uint32 memI;
	for(memI = 0; memI < mem_list.size(); memI++)
		if(MemHashEntry::start_lessthan(mhe, mem_list[memI]))
			break;
	mem_list.insert(mem_list.begin()+memI, mhe);

	m_baseLength = 0;
	for(uint32 i=0; i < mem_list.size(); i++)
		m_baseLength += mem_list[i].Length();

	m_gapLength = 0;
	uint32 fs = mem_list[0].FirstStart();
	for(uint32 i=1; i < mem_list.size(); i++)
		m_gapLength += mem_list[i][fs] - mem_list[i-1][fs] - mem_list[i-1].Length();
}

int MimHashEntry::start_compare(MemHashEntry& mhe, uint32 tolerance){
	return MemHashEntry::fuzzy_align(mhe, mem_list[0], tolerance);
}
//returns -1 if the start position of mhe is < the end position of
//this mim.  returns 0 if the end position of this mim is within the
//tolerable range of the start of the next mem.  returns 1 if this
//mim ends before the tolerable range of the mem's start
//otherwise returns -1
int MimHashEntry::end_compare(MemHashEntry& mhe, uint32 tolerance){
	return MemHashEntry::fuzzy_align(mem_list[mem_list.size()-1], mhe, tolerance);
}

gnSeqI MimHashEntry::BaseLength() const{
	return m_baseLength;
}

gnSeqI MimHashEntry::GapLength() const{
	return m_gapLength;
}

ostream& operator<<(ostream& os, const MimHashEntry& mhe){ //write to stream.
	os << mhe.mem_list.size() << '\t';
	os << mhe.BaseLength() << '\t' << mhe.GapLength() << '\t';
	vector<MemHashEntry>::const_iterator mlist = mhe.mem_list.begin();
	uint32 m_seq_count = mlist->SeqCount();
	if(m_seq_count > 0)
		os << (*mlist)[0];
	for(uint32 i=1; i < m_seq_count; i++)
		os << '\t' << (*mlist)[i];
	mlist = mhe.mem_list.end() - 1;
	gnSeqI last_len = mlist->Length();
	for(uint32 i=0; i < m_seq_count; i++)
		os << '\t' << (*mlist)[i] + last_len;
	return os;
}

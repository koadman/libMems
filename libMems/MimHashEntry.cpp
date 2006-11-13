#include "MimHashEntry.h"
#include "MemHashEntry.h"
#include "MimHash.h"
#include "gn/gnFilter.h"
#include <algorithm>

MimHashEntry::MimHashEntry(){
}

MimHashEntry::MimHashEntry( const MimHashEntry& mh){
	mem_list = mh.mem_list;
	m_start = mh.m_start;
	m_end = mh.m_end;
	m_baseLength = mh.m_baseLength;
}

MimHashEntry* MimHashEntry::Clone() const{
	return new MimHashEntry(*this);
}

void MimHashEntry::AddMem(MemHashEntry* mhe){
	//add in order
//	uint32 memI;
//	for(memI = 0; memI < mem_list.size(); memI++)
//		if(MemHashEntry::start_lessthan(*mhe, *(mem_list[memI])))
//			break;

// let the set class do the sorting for us
	mem_list.insert(mhe);
	mhe->SetMim( this );

	set<MemHashEntry*, MheCompare>::iterator mem_iter = mem_list.begin();
	
	//this code is all useless now that mems can overlap each other
	m_baseLength = 0;
	for(; mem_iter != mem_list.end(); mem_iter++)
		m_baseLength += (*mem_iter)->Length();

}

void MimHashEntry::Merge(MimHashEntry& mhe){
	// let the set class do the merge work
	mem_list.insert(mhe.mem_list.begin(), mhe.mem_list.end());
//	set<MemHashEntry*, MheCompare>::iterator old_iter = mhe.mem_list.begin();
//	for(; old_iter != mhe.mem_list.end(); old_iter++)
//		(*old_iter)->m_mim = this;
}

void MimHashEntry::DeleteMem(MemHashEntry* mhe){
	mem_list.erase(mhe);
}

int MimHashEntry::start_compare(MemHashEntry& mhe, uint32 tolerance){
	return MemHashEntry::fuzzy_align(mhe, **(mem_list.begin()), tolerance);
}
//returns -1 if the start position of mhe is < the end position of
//this mim.  returns 0 if the end position of this mim is within the
//tolerable range of the start of the next mem.  returns 1 if this
//mim ends before the tolerable range of the mem's start
//otherwise returns -1
int MimHashEntry::end_compare(MemHashEntry& mhe, uint32 tolerance){
	set<MemHashEntry*, MheCompare>::const_iterator end_mem = mem_list.end();
	return MemHashEntry::fuzzy_align(**(--end_mem), mhe, tolerance);
}

gnSeqI MimHashEntry::BaseLength() const{
	return m_baseLength;
}


void MimHashEntry::PrintMismatches(ostream& os, vector<string>& seq_strings, uint32 context){
	vector<MemHashEntry*> mie_list;
	set<MemHashEntry*, MheCompare>::iterator iter = mem_list.begin();
	const gnFilter* revCompFilter = gnFilter::DNAComplementFilter();
	for(; iter != mem_list.end(); iter++)
		//only do n-way matches
		if((*iter)->SeqCount() == seq_strings.size())
			mie_list.push_back(*iter);
	if(mie_list.size() < 2)
		return; //no sense in printing mismatches in less than 2 mems!
	sort(mie_list.begin(), mie_list.end(), &MemHashEntry::start_lessthan_ptr);
	int64* diff_vector = new int64[seq_strings.size()];
	int32 max_diff = INT32_MIN;
	int32 min_diff = INT32_MAX;
	boolean too_far;

	for(uint32 memI = 0; memI < mie_list.size() - 1; memI++){
		too_far = false;
		max_diff = INT32_MIN;
		min_diff = INT32_MAX;
		for(uint32 seqI = 0; seqI < seq_strings.size(); seqI++){
			try{
				diff_vector[seqI] = mie_list[memI]->GapSize(*mie_list[memI+1], seqI);

				if(abs( (long) diff_vector[seqI] ) > 1){
					too_far = true;
					break;
				}
				if(diff_vector[seqI] > max_diff)
					max_diff = diff_vector[seqI];
				if(diff_vector[seqI] < min_diff)
					min_diff = diff_vector[seqI];
			}catch(gnException& e){
				cout << "EEEEEEEEERRRROOR!\n";
				diff_vector[seqI] = 0;
				too_far = true;
				break;
			}
		}
		if(too_far)
			continue;
		if((max_diff == min_diff) && (min_diff < 0)){
			os << "Bug alert: Even overlap, min diff = " << min_diff << "\n";
			os << "Bad Starting: " << *mie_list[memI] << "\n";
			os << "Bad Ending:   " << *mie_list[memI+1] << "\n";
			continue;
		}
//		os << "Starting: " << *mie_list[memI] << "\n";
//		os << "Ending:   " << *mie_list[memI+1] << "\n";
//		os << "Min diff: " << min_diff << "\n";
//		os << "Max diff: " << max_diff << "\n";
		for(uint32 seqI = 0; seqI < seq_strings.size(); seqI++){
			int64 sec_start = (*mie_list[memI+1])[seqI];
			if(sec_start < 0)
				sec_start = -(*mie_list[memI])[seqI];
			gnSeqI print_offset = sec_start - diff_vector[seqI] - context - 1;
			if(diff_vector[seqI] < max_diff){
				string toPrint = seq_strings[seqI].substr( print_offset, context*2);
				if((*mie_list[memI+1])[seqI] < 0 && context > 0){
					uint32 dumb_size = toPrint.size();
					gnSeqC* dumb_bug = new gnSeqC[dumb_size + 1];
					memcpy(dumb_bug, toPrint.c_str(), dumb_size);
					revCompFilter->ReverseFilter(&dumb_bug, dumb_size);
					dumb_bug[toPrint.size()] = 0;
					toPrint = dumb_bug;
					delete[] dumb_bug;
				}

				toPrint.insert(context, "-");
				os << toPrint;
			}else{
				string toPrint = seq_strings[seqI].substr( print_offset, context*2 + 1);
				if((*mie_list[memI+1])[seqI] < 0){
					uint32 dumb_size = toPrint.size();
					gnSeqC* dumb_bug = new gnSeqC[dumb_size + 1];
					memcpy(dumb_bug, toPrint.c_str(), dumb_size);
					revCompFilter->ReverseFilter(&dumb_bug, dumb_size);
					dumb_bug[toPrint.size()] = 0;
					toPrint = dumb_bug;
					delete[] dumb_bug;
				}

				os << toPrint;
			}
			os << '\t';
		}
		os << *this << '\n';
	}
	delete[] diff_vector;
}

ostream& operator<<(ostream& os, const MimHashEntry& mhe){ //write to stream.
	os << mhe.mem_list.size() << '\t';
	os << mhe.BaseLength() << '\t';
	set<MemHashEntry*, MheCompare>::const_iterator mlist = mhe.mem_list.begin();
	uint32 m_seq_count = (*mlist)->SeqCount();
	if(m_seq_count > 0)
		os << (**mlist)[0];
	for(uint32 i=1; i < m_seq_count; i++)
		os << '\t' << (**mlist)[i];
	mlist = mhe.mem_list.end();
	mlist--;
	gnSeqI last_len = (*mlist)->Length();
	for(uint32 i=0; i < m_seq_count; i++)
		os << '\t' << (**mlist)[i] + last_len;
	return os;
}

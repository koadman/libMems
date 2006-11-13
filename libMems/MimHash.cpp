#include "MimHash.h"
#include "MimHashEntry.h"
#include <algorithm>

MimHash::MimHash(const MimHash& mh){
	m_mim_list = mh.m_mim_list;
}

MimHash* MimHash::Clone() const{
	return new MimHash(*this);
}

MemCollapseEntry::MemCollapseEntry(){
	next_start = NULL;
	prev_start = NULL;
	next_end = NULL;
	prev_end = NULL;
}

MemCollapseEntry::MemCollapseEntry(MemHashEntry& mhe){
	next_start = NULL;
	prev_start = NULL;
	next_end = NULL;
	prev_end = NULL;
	m_mem = mhe;
}

void MemCollapseEntry::DeleteSelf(){
	if(prev_start != NULL)
		prev_start->next_start = next_start;
	if(next_start != NULL)
		next_start->prev_start = prev_start;
	if(prev_end != NULL)
		prev_end->next_end = next_end;
	if(next_end != NULL)
		next_end->prev_end = prev_end;
}
//worst case for a 1bp indel is a distance of seq_count - 1 hash buckets
//scanning range for the adjacent mems can be calculated as
//mismatch_range = tolerance * (seq_count - 1)
//try sliding a window of size mismatch_range centered on the current
//across the hash table bucket across the hash tables.  sort the
//entries in the window on their starting base pair and look for the
//first aligned mem following the current mem
//possibly look for snp's at this stage too..
//try keeping a list iterator for each bucket in the window which
//points to the first possible mem which could follow the current

//simpler approach:  extract all the mems to a single list.  sort the
//list on start positions and proceed sequentially through the list
//extending matches as far as possible in each direction.  For each
//mem, seek past its end in the list and then check each successive
//mem first if its generalized offset falls within mismatch_range, and
//then if each start coordinate is within the tolerance range of
//the original mem's end coordinates.  if it satisfies these two
//conditions then it is added as part of a mim.
//extending forward is easy in this manner, but extending backwards
//is not.  what I need is a doubly linked list which is linked in
//the forward direction based on start positions and in the reverse
//direction based on end positions
//how the fuck do you create one of these?
//put everything into an array and sort it. Set the 'next' pointers
//to their correct values.  now sort the array
//on end positions and link everything up.  the first and last
//elements will be different in each direction.
boolean MimHash::ResolveMismatches( uint32 tolerance ){
	vector<MemCollapseEntry*> mem_del_array;
	MemCollapseEntry* forward_begin, *forward_end, *reverse_begin, *reverse_end;
	list<MemHashEntry>::iterator mhe;
	mem_del_array.reserve(m_mem_count);
	for(uint32 i=0; i < table_size; i++){
		mhe = mem_table[i].begin();
		for(; mhe != mem_table[i].end(); mhe++){
			MemCollapseEntry* mce = new MemCollapseEntry(*mhe);
			mem_del_array.push_back(mce);
		}
	}

	uint32 i;
	MemHashEntry::seq_compare_start = 0;
	sort(mem_del_array.begin(), mem_del_array.end(), &MemCollapseEntry::start_lessthan_ptr);
	mem_del_array[0]->prev_start = NULL;
	for(i = 1; i < mem_del_array.size(); i++){
		mem_del_array[i - 1]->next_start = mem_del_array[i];
		mem_del_array[i]->prev_start = mem_del_array[i - 1];
	}
	mem_del_array[i-1]->next_start = NULL;
	forward_begin = mem_del_array[0];
	forward_end = mem_del_array[i-1];

	sort(mem_del_array.begin(), mem_del_array.end(), &MemCollapseEntry::end_lessthan_ptr);
	mem_del_array[0]->prev_end = NULL;
	for(i = 1; i < mem_del_array.size(); i++){
		mem_del_array[i-1]->next_end = mem_del_array[i];
		mem_del_array[i]->prev_end = mem_del_array[i-1];
	}
	mem_del_array[i-1]->next_end = NULL;
	reverse_begin = mem_del_array[0];
	reverse_end = mem_del_array[i - 1];
	
	MemCollapseEntry *next_mce, *prev_mce, *cur_mce = forward_begin;
	while(cur_mce != NULL){
		MimHashEntry cur_mim;
		cur_mim.AddMem(cur_mce->m_mem);
		next_mce = cur_mce->next_start;

		//check for fuzzy alignment in the forward direction
		while(next_mce != NULL && cur_mim.end_compare(next_mce->m_mem, 0) < 0)
			next_mce = next_mce->next_start;
		int end_comp = 0;
		while(next_mce != NULL && end_comp <= 0){
			end_comp = cur_mim.end_compare(next_mce->m_mem, tolerance);
			if(end_comp == 0){
				cur_mim.AddMem(next_mce->m_mem);
				//delete the mem
				next_mce->DeleteSelf();
			}
			next_mce = next_mce->next_start;
		}
		
		//now check for fuzzy alignment in the reverse direction
		prev_mce = cur_mce->prev_end;
		while(prev_mce != NULL && cur_mim.start_compare(prev_mce->m_mem, 0) < 0)
			prev_mce = prev_mce->prev_end;
		int start_comp = 0;
		while(prev_mce != NULL && start_comp <= 0){
			start_comp = cur_mim.start_compare(prev_mce->m_mem, tolerance);
			if(start_comp == 0){
				cur_mim.AddMem(prev_mce->m_mem);
				//delete the mem
				prev_mce->DeleteSelf();
			}
			prev_mce = prev_mce->prev_end;
		}
		m_mim_list.push_back(cur_mim);
		//delete the mem
		cur_mce->DeleteSelf();
		cur_mce = cur_mce->next_start;
	}
	
	
	return true;
}

boolean MimHash::GetMimList( vector<MimHashEntry>& mim_list ){
	mim_list = m_mim_list;
	return true;
}

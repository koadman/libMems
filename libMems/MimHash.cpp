#include "MimHash.h"
#include "MimHashEntry.h"
#include <algorithm>
#include <fstream>

MemInclusionEntry::MemInclusionEntry(const MemInclusionEntry& mie) : MemHashEntry(mie){
	m_overlaps = mie.m_overlaps;
	m_underlaps = mie.m_underlaps;
}

MemInclusionEntry::MemInclusionEntry(const MemHashEntry& mie) : MemHashEntry(mie){
}

MemInclusionEntry& MemInclusionEntry::operator=(const MemInclusionEntry& mie){
	MemHashEntry::operator=(mie);
	m_overlaps = mie.m_overlaps;
	m_underlaps = mie.m_underlaps;
	return *this;
}

MemInclusionEntry& MemInclusionEntry::operator=(const MemHashEntry& mie){
	MemHashEntry::operator=(mie);
	return *this;
}

MemInclusionEntry* MemInclusionEntry::Clone() const{
	return new MemInclusionEntry(*this);
}

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

MemCollapseEntry::MemCollapseEntry(MemHashEntry* mie){
	next_start = NULL;
	prev_start = NULL;
	next_end = NULL;
	prev_end = NULL;
	m_mie = mie;
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

void MimHash::ResolveMismatches( uint32 tolerance ){
	vector<MemCollapseEntry*> mem_del_array;
	set<MemHashEntry*, MheCompare>::iterator mem_table_iter;
	MemCollapseEntry* forward_begin, *forward_end, *reverse_begin, *reverse_end;
	
	mem_del_array.reserve(m_mem_count);

	for(uint32 bucketI = 0; bucketI < table_size; bucketI++){
		mem_table_iter = mem_table[bucketI].begin();
		for(; mem_table_iter != mem_table[bucketI].end(); mem_table_iter++){
			if((*mem_table_iter)->Multiplicity() == seq_count){
				MemCollapseEntry* mce = new MemCollapseEntry(*mem_table_iter);
				mem_del_array.push_back(mce);
			}
		}
	}

	uint32 i;

	//sort on the start positions of the first sequence
	MemHashEntry::seq_compare_start = 0;
	sort(mem_del_array.begin(), mem_del_array.end(), &MemCollapseEntry::start_lessthan_ptr);

	//link each mem up in forward order
	mem_del_array[0]->prev_start = NULL;
	for(i = 1; i < mem_del_array.size(); i++){
		mem_del_array[i - 1]->next_start = mem_del_array[i];
		mem_del_array[i]->prev_start = mem_del_array[i - 1];
	}
	mem_del_array[i-1]->next_start = NULL;
	forward_begin = mem_del_array[0];
	forward_end = mem_del_array[i-1];

	//sort on the end positions of the first sequence
	sort(mem_del_array.begin(), mem_del_array.end(), &MemCollapseEntry::end_lessthan_ptr);
	//link each mem up in reverse order
	mem_del_array[0]->prev_end = NULL;
	for(i = 1; i < mem_del_array.size(); i++){
		mem_del_array[i-1]->next_end = mem_del_array[i];
		mem_del_array[i]->prev_end = mem_del_array[i-1];
	}
	mem_del_array[i-1]->next_end = NULL;
	reverse_begin = mem_del_array[0];
	reverse_end = mem_del_array[i - 1];
	
	//starting with the first one in the list
	MemCollapseEntry *next_mce, *prev_mce, *cur_mce = forward_begin;
	while(cur_mce != NULL){
		MimHashEntry cur_mim;
		cur_mim.AddMem(cur_mce->m_mie);
		next_mce = cur_mce->next_start;

		//check for fuzzy alignment in the forward direction

		//first scan past the end of the mim to avoid overlaps
		while(next_mce != NULL && cur_mim.end_compare(*(next_mce->m_mie), 0) < 0)
			next_mce = next_mce->next_start;

		int end_comp = 0;
		//search for the nearest mem which starts within the specified tolerance range in each genome
		while(next_mce != NULL && end_comp <= 0){
			end_comp = cur_mim.end_compare(*(next_mce->m_mie), tolerance);
			if(end_comp == 0){
				cur_mim.AddMem(next_mce->m_mie);
				//delete the mem
				next_mce->DeleteSelf();
			}
			next_mce = next_mce->next_start;
		}
		
		//now check for fuzzy alignment in the reverse direction
		prev_mce = cur_mce->prev_end;
		while(prev_mce != NULL && cur_mim.start_compare(*(prev_mce->m_mie), 0) < 0)
			prev_mce = prev_mce->prev_end;
		int start_comp = 0;
		while(prev_mce != NULL && start_comp <= 0){
			start_comp = cur_mim.start_compare(*(prev_mce->m_mie), tolerance);
			if(start_comp == 0){
				cur_mim.AddMem(prev_mce->m_mie);
				//delete the mem
				prev_mce->DeleteSelf();
			}
			prev_mce = prev_mce->prev_end;
		}
		m_mim_list.push_back(cur_mim.Clone());
		//delete the mem
		cur_mce->DeleteSelf();
		cur_mce = cur_mce->next_start;
	}
}


void MimHash::GetMimList( vector<MimHashEntry>& mim_list ){
	for(uint32 i=0; i < m_mim_list.size(); i++)
		mim_list.push_back(*m_mim_list[i]);
}


void MimHash::PrintMismatches(ostream& os, vector<gnSequence*>& seqs){
	vector<string> seq_strings;
	for(uint32 i=0; i < seqs.size(); i++){
		seq_strings.push_back(string());
		seqs[i]->ToString(seq_strings[i]);
	}
	for(uint32 mimI = 0; mimI < m_mim_list.size(); mimI++){
		m_mim_list[mimI]->PrintMismatches(os, seq_strings, 0);
	}
}

void MimHash::EliminateInclusions() {
}


// Locates overlapping subset matches in a set of mems. 
void MimHash::LocateSubsetOverlaps(){
	//can't find subset matches in two genomes!
	if(seq_count <= 2)
		return;
	
	uint32 debug_level = 100;
	boolean debug_algorithm = false;
	
	set<MemHashEntry*, MheCompare>::iterator mem_table_iter;
	list<MemHashEntry*>::iterator mhe, mem_iter, prev_iter;
	//define an array of lists
	//first index is the number of sequences in the match
	//second index is the sequence that the list is sorted on.
	list<MemHashEntry*>** mhe_list;

	uint32 list_count = seq_count - 1;
	mhe_list = new list<MemHashEntry*>*[list_count];
	for(uint32 sdf = 0; sdf < list_count; sdf++)
		mhe_list[sdf] = new list<MemHashEntry*>[list_count];
	MemHashEntry* new_mie;

	//get all the mems out of the hash table and add them to the list which
	//corresponds to their match count in the first sort level.
	for(uint32 bucketI = 0; bucketI < table_size; bucketI++){
		mem_table_iter = mem_table[bucketI].begin();
		for(; mem_table_iter != mem_table[bucketI].end(); mem_table_iter++){
			mhe_list[new_mie->Multiplicity()-2][0].push_back( *mem_table_iter );
		}
		mem_table[bucketI].clear();
	}

	for(uint32 bellI = list_count; bellI > 0; bellI--){
		cout << "There are " << mhe_list[bellI-1][0].size() << " " << bellI+1 << "-way matches\n";
	}

	//keep a vector of iterators with the current positions to check in
	//each level of sorted match lists
	vector<list<MemHashEntry*>::iterator>* start_iters;
	vector<list<MemHashEntry*>::iterator>* cur_starts;
	//first index is the sequence that mems are sorted on
	//second index is a pointer to each successive match level's start level in that list.
	//got it?  good.
	start_iters = new vector<list<MemHashEntry*>::iterator>[list_count];
	cur_starts = new vector<list<MemHashEntry*>::iterator>[list_count];
	
	//sort the highest match level on each sequence.
	mhe_list[list_count - 1][0].sort(&MemHashEntry::start_lessthan_ptr);
	for(uint32 seq_matchI = 0; seq_matchI < list_count; seq_matchI++){
		if(seq_matchI > 0){
			mhe_list[list_count - 1][seq_matchI] = mhe_list[list_count - 1][0];
			MemHashEntry::seq_compare_start = seq_matchI;
			mhe_list[list_count - 1][seq_matchI].sort(&MemHashEntry::start_lessthan_ptr);
		}

		//add the new start_iters for the current level...
		mem_iter = mhe_list[list_count - 1][seq_matchI].begin();
		prev_iter = mem_iter;
		start_iters[seq_matchI].push_back(mem_iter);

		//the last start iter on each level signifies the end of the level
		start_iters[seq_matchI].push_back(mhe_list[list_count - 1][seq_matchI].end());
		cur_starts[seq_matchI] = start_iters[seq_matchI];
	}

	//now for each sequence match level starting from the second highest,
	//look for inclusions in each higher level
	for(int32 listI = list_count - 2; listI >= 0; listI--){
		//sort the current level
		MemHashEntry::seq_compare_start = 0;
		mhe_list[listI][0].sort(&MemHashEntry::start_lessthan_ptr);

		//find inclusions starting with the highest sequence match level
		for(int32 checkI = list_count - 1; checkI > listI; checkI--){
			//the following two vars are the ranges of start iterators to check against in this match level
			uint32 cur_minI = quadratic_li(list_count - checkI - 1) + (list_count - checkI - 1);
			uint32 cur_maxI = quadratic_li(list_count - checkI) + (list_count - checkI);

			boolean erased_start = false;
			mem_iter = mhe_list[listI][0].begin();
			uint32 cur_start_level;
			//check against every mem in the current match level
			for(; mem_iter != mhe_list[listI][0].end(); mem_iter++){
				//if we erased the first mem in the list then we need special handling
				if(erased_start){
					mem_iter--;
					erased_start = false;
				}
				cur_start_level = (*mem_iter)->FirstStart();
				if(cur_start_level >= debug_level)
					debug_algorithm = true;
				//look for each sequence start level in the current sequence match level
				if(debug_algorithm)
					cout << "Checking " << **mem_iter << "\n";
				list<MemHashEntry*> overlap_list;
				//note: we don't always need to scan through this entire loop.  It should be fixed such that
				//each if FirstStart() of the cur_start is > FirstStart() of mem_iter we skip the loop.
				list<MemHashEntry*>::iterator cur_start_iter, next_start_iter, matching_iter;
				for(uint32 cur_startI = cur_minI; cur_startI < cur_maxI - 1; cur_startI++){
					cur_start_iter = cur_starts[cur_start_level][cur_startI];
					next_start_iter = cur_starts[cur_start_level][cur_startI + 1];
					//scan only until the cur_starts start after mem_iter does.
					if(cur_start_level < (*cur_start_iter)->FirstStart())
						break;

					if(debug_algorithm)
						if(cur_start_iter != next_start_iter)
							cout << "Against " << **cur_start_iter << "\n";
					while(cur_start_iter != next_start_iter &&
					      MemHashEntry::start_compare(**cur_start_iter, **mem_iter) < 0){
						cur_start_iter++;
						if(debug_algorithm)
							if(cur_start_iter != next_start_iter)
								cout << "Against " << **cur_start_iter << "\n";
					}
					if(cur_start_iter == next_start_iter)
						continue;

					// loop till we've scanned the entire range of this mem.
					// each time an even overlap is discovered link the two mems.
					// this is the easy way to do it.  there may be trickier ways
					// which get the job done faster.
					matching_iter = cur_start_iter;
					MemHashEntry tmp_entry = **mem_iter;
					tmp_entry.CropStart(tmp_entry.Length()-1);
					while(matching_iter != next_start_iter &&
						MemHashEntry::start_compare(tmp_entry, **matching_iter) >= 0){

						//now check for intersections
						if((*mem_iter)->EvenlyOverlaps(**matching_iter)){
							// they evenly overlap.  link them to each other.
							(*matching_iter)->LinkSubset(*mem_iter);
						}
						matching_iter++;
					}
					//now check for intersections
// don't know what to do with the following shite.
/*
					if((*mem_iter)->EvenlyOverlaps(**cur_start_iter)){
						//add the mem above it to the intersect list
//						intersect_list.push_back(**cur_starts[cur_start_level][cur_startI]);
						
						//delete the inclusion and insert fragments
						MemHashEntry unique_mhe, *new_mie;
						MimHashEntry* mimmie = (*cur_start_iter)->m_mim;
						list<MemHashEntry*>::iterator tmp_iter = mem_iter;
						MemHashEntry bad_mie, debug_mie = **cur_start_iter;
						bad_mie = **mem_iter;
						//insert a unique beginning fragment, if it exists
						if((*mem_iter)->GetUniqueStart(**cur_start_iter, unique_mhe)){
							new_mie = unique_mhe.Clone();
							SetMimPointers( *mem_iter, new_mie, *cur_start_iter);
							// set mem_iter's mim ptr to NULL to avoid confusion if it also has an end segment
							(*mem_iter)->m_mim = NULL;
							if(new_mie->Length() > 140000)
								cout << "weekesas\n";
							mhe_list[listI][0].insert(mem_iter, new_mie);
							if(debug_algorithm)
								cout << "Added unique start " << unique_mhe << "\n";
						}

						debug_mie = **cur_start_iter;
						bad_mie = **mem_iter;
						
						//now insert the end fragment
						if((*mem_iter)->GetUniqueEnd(**cur_start_iter, unique_mhe)){
							new_mie = unique_mhe.Clone();
							SetMimPointers( *mem_iter, new_mie, *cur_start_iter);

							//scan forward to find the insert point for the end mem
							while(tmp_iter != mhe_list[listI][0].end() && 
							      (*tmp_iter)->FirstStart() == unique_mhe.FirstStart() &&
							      MemHashEntry::start_compare(**tmp_iter, unique_mhe) < 0){
								tmp_iter++;
							}
							if(new_mie->Length() > 140000)
								cout << "weekesas\n";
							mhe_list[listI][0].insert(tmp_iter, new_mie);
							if(debug_algorithm)
								cout << "Added unique end " << unique_mhe << "\n";
						}
						
						//delete the middle piece
						tmp_iter = mem_iter++;
						mhe_list[listI][0].erase(tmp_iter);
						if(mem_iter != mhe_list[listI][0].begin())
							mem_iter--;
						else
							erased_start = true;
					}
*/
					//commit any modifications to our array of current starts to the start iter list.
					cur_starts[cur_start_level][cur_startI] = cur_start_iter;
				}
			}
		}
		
		for(uint32 seq_matchI = 0; seq_matchI < list_count; seq_matchI++){
			if(seq_matchI > 0){
				mhe_list[listI][seq_matchI] = mhe_list[listI][0];
				MemHashEntry::seq_compare_start = seq_matchI;
				mhe_list[listI][seq_matchI].sort(&MemHashEntry::start_lessthan_ptr);
			}
	
			//add the new start_iters for the current level...
			mem_iter = mhe_list[listI][seq_matchI].begin();
			prev_iter = mem_iter;
			start_iters[seq_matchI].push_back(mem_iter);
			uint32 added_iter_count = 1;
			for(mem_iter++; mem_iter != mhe_list[listI][seq_matchI].end();  mem_iter++){
				uint32 diff = (*mem_iter)->FirstStart() - (*prev_iter)->FirstStart();
				while(diff > 0){
					start_iters[seq_matchI].push_back(mem_iter);
					added_iter_count++;
					diff--;
				}
				prev_iter++;
			}
			for(; added_iter_count < list_count - listI; added_iter_count++)
				start_iters[seq_matchI].push_back(mhe_list[listI][seq_matchI].end());
	
			//the last start iter on each level signifies the end of the level
			start_iters[seq_matchI].push_back(mhe_list[listI][seq_matchI].end());
			cur_starts[seq_matchI] = start_iters[seq_matchI];
		}
	}
	
	list<MemHashEntry*>* alignment_mem_list;
	alignment_mem_list = new list<MemHashEntry*>[list_count];

	//now hash all the mems again
	uint32 mimless_count = 0;
	for(uint32 i=0; i < list_count; i++){
		cout << "Now there are " << mhe_list[i][0].size() << " " << i+2 << "-way matches\n";
		mhe = mhe_list[i][0].begin();
		for(; mhe != mhe_list[i][0].end(); mhe++){
			(*mhe)->SetExtended(true);
			if((*mhe)->GetSupersets().size() == 0 && (*mhe)->GetSubsets().size() == 0)
				mimless_count++;
			if(!AddHashEntry(**mhe))
				cout << "You got problems bustah!\n";
			alignment_mem_list[i].push_back(*mhe);
//			delete *mhe;
		}
	}
	
	//Print some stats on the mims that were created.
	cout << "Created " << m_mim_list.size() << " inexact matches\n";
	cout << "There remain " << mimless_count << " mems without mims.\n";
	
//	MemAligner ma(seq_count);
//	ma.CreateAlignment(alignment_mem_list);

// time for the cleanup crew
	for(uint32 fdfd = 0; fdfd < list_count; fdfd++)
		delete[] mhe_list[fdfd];
	delete[] mhe_list;
}

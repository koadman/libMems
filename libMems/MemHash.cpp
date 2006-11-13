#include "MemHash.h"
#include <list>

MemHash::MemHash(uint32 tablesize){
	table_size = tablesize;
	mer_size = 0;
	seq_count = 0;
	m_mem_count = 0;
	m_collision_count = 0;
	m_repeat_tolerance = DEFAULT_REPEAT_TOLERANCE;
	m_enumeration_tolerance = DEFAULT_ENUMERATION_TOLERANCE;
	//allocate the hash table
	mem_table = new list<MemHashEntry>[table_size];
	for(uint32 i=0; i < table_size; i++)
		mem_table_count.push_back(0);
}

//make sure this calls the destructor on each element
MemHash::~MemHash(){
	if(mem_table != NULL)
		delete[] mem_table;
}

MemHash::MemHash(const MemHash& mh){
	table_size = mh.table_size;
	mer_size = mh.mer_size;
	seq_count = mh.seq_count;
	m_mem_count = mh.m_mem_count;
	m_collision_count = mh.m_collision_count;
	m_repeat_tolerance = mh.m_repeat_tolerance;
	m_enumeration_tolerance = mh.m_enumeration_tolerance;
	//allocate the hash table
	mem_table = new list<MemHashEntry>[table_size];
	for(uint32 i=0; i < table_size; i++){
		mem_table_count.push_back(mh.mem_table_count[i]);
		mem_table[i] = mh.mem_table[i];
	}
}

MemHash* MemHash::Clone() const{
	return new MemHash(*this);
}

/*boolean MemHash::ResolveAmbiguities( uint32 tolerance ){
	list<MemHashEntry>::iterator cur_mhe;
	list<MemHashEntry>::iterator next_mhe;
	list<MemHashEntry> start_list;
	int offset;
	boolean aligned;
	int64 diff;
	for(uint32 i=0; i < table_size; i++){
		cur_mhe = mem_table[i].begin();
		start_list = mem_table[i];
		MemHashEntry::seq_compare_start = 0;
		start_list.sort(&start_lessthan);
		for(; cur_mhe != start_list.end(); mhe++){
			next_mhe = cur_mhe;
			uint32 startI = cur_mhe.FirstStart();
			while(next_mhe != start_list.end()  && 
			      next_mhe[startI] < cur_mhe[startI]+cur_mhe.length() || 
			      next_mhe[startI] == MEM_NO_MATCH){
				next_mhe++;
			}
			if(next_mhe == start_list.end())
				continue;
			
			//check if the next mem is within the tolerance range
			//if so then check if it's aligned.
		   cur_mhe.DetailCompare(next_mhe, offset, aligned, diff);
			if(aligned && diff < cur_mhe.length() + tolerance){
				//is the difference caused by an ambiguity or an actual difference?
				for(uint32 mismatchI = 0; mismatchI < diff - cur_mhe.length(); mismatchI++){
					gnSequence bad_base, bad_base2;
					uint32 seqI = 0;
					for(; seqI < seq_count; seqI++){
						if(cur_mhe[seqI] == MEM_NO_MATCH)
							continue;
						bad_base = seq_table[seqI].subseq(cur_mhe[seqI]+cur_mhe.length() + mismatchI);
					}
					for(seqI++; seqI < seq_count; seqI++){
						if(cur_mhe[seqI] == MEM_NO_MATCH)
							continue;
						if(cur_mhe[seqI
						bad_base2 = seq_table[seqI].subseq(cur_mhe[seqI]+cur_mhe.length() + mismatchI);
						if(bad_base.compare(bad_base2) == 0)
							//ambiguity.  extend across it.
							
					}
				}
				//collapse-o-matic!
				
			}
			//
		}
	}
}
*/
boolean MemHash::Create(uint32 mer_mask_size){
	return FindMatches(mer_mask_size, 0, GNSEQI_END);
//	EliminateInclusions();
}

boolean MemHash::GetMemList( vector<MemHashEntry>& mem_list ){
	list<MemHashEntry>::iterator mhe;
	for(uint32 i=0; i < table_size; i++){
		mhe = mem_table[i].begin();
		for(; mhe != mem_table[i].end(); mhe++)
			mem_list.push_back(*mhe);
	}
	return true;
}

boolean MemHash::EnumerateMatches( list<idmer>& match_list ){
	//if we aren't allowing repeats only hash the mers which occur
	//fewer than n times.
	//this may still allow repeats if the matches are between subsets
	//of the genomes..
	if(m_repeat_tolerance == 0 && match_list.size() <= seq_count)
		return HashMatch(match_list);
	//the better soln. is to have a repeat tolerance of 1.
	return EnumerateRepeats( match_list );
}

boolean MemHash::EnumerateRepeats( list<idmer>& match_list ){
	match_list.sort(&idmer_id_lessthan);
	list<idmer>::iterator iter = match_list.begin();
	list<idmer>::iterator iter2 = match_list.begin();
	list<idmer> hash_list;
	iter2++;
	int32 cur_tolerance = m_enumeration_tolerance - 1;
	uint32 cur_count = 1;
	hash_list.push_back(*iter);
	for(; iter2 != match_list.end(); iter++){
		if(cur_count > m_repeat_tolerance)
			return true;
		if(iter->id != iter2->id){
			hash_list.push_back(*iter2);
			cur_tolerance = m_enumeration_tolerance - 1;
			cur_count = 1;
		}else{
			if(cur_tolerance > 0){
				hash_list.push_back(*iter);
				cur_tolerance--;
			}
			cur_count++;
		}
		iter2++;
	}

	if(m_enumeration_tolerance == 1)
		return HashMatch(hash_list);

	return MatchFinder::EnumerateMatches( hash_list );
}

//why have separate hash tables?
// MemHashEntries use GENETICIST coordinates.  They start at 1, not 0.
boolean MemHash::HashMatch(list<idmer>& match_list){
	//check that there is at least one forward component
	match_list.sort(&idmer_id_lessthan);
	// initialize the hash entry
	MemHashEntry mhe(seq_count, mer_size);
	mhe.SetLength(mer_size);
	
	//Fill in the new MemHashEntry and set direction parity if needed.
	list<idmer>::iterator iter = match_list.begin();
	for(; iter != match_list.end(); iter++)
		mhe.SetStart(sarid_table[(uint8)iter->id], iter->position + 1);
	SetDirection(mhe);
	mhe.CalculateOffset();
	if(mhe.SequenceCount() < 2){
		cout << "red fag\n";
	}else if(AddHashEntry(mhe)){
		//memhe was added successfully, extend it
//		cout << "Added " << mhe << "\tOffset " << mhe.Offset() << "\n";
	}
	return true;
}

void MemHash::SetDirection(MemHashEntry& mhe){
	//get the reference direction
	boolean ref_forward;
	uint32 seqI=0;
	for(; seqI < seq_count; seqI++)
		if(mhe[seqI] != MEM_NO_MATCH){
			ref_forward = !(sar_table[seqI]->GetMer(mhe[seqI] - 1) & 0x1);
			break;
		}
	//set directional parity for the rest
	for(seqI++; seqI < seq_count; seqI++)
		if(mhe[seqI] != MEM_NO_MATCH)
			if(ref_forward == (sar_table[seqI]->GetMer(mhe[seqI] - 1) & 0x1))
				mhe.SetStart(seqI, -mhe[seqI]);
}

//returns true if the mem was successfully hashed, false otherwise
//in the case of a submem it will always return true.
//the reference genome is now the first one with an entry in the start array.
boolean MemHash::AddHashEntry(MemHashEntry& mhe){
	//first compute which hash table bucket this is going into
	int64 offset = mhe.Offset();

	uint32 bucketI = ((offset % table_size) + table_size) % table_size;
	list<MemHashEntry>::iterator cur_he = mem_table[bucketI].begin();
	list<MemHashEntry>::iterator end_entry = mem_table[bucketI].end();
	//now check for intersection in all hash entries in this bucket
	for(; cur_he != end_entry; cur_he++){
		int64 coffset = cur_he->Offset();
		if(coffset < offset){
			break;
		}else if(coffset == offset){
			//offsets are the same, check for containment...
			if(cur_he->Contains(mhe)){
				// it's a valid intersection.  return false;
				m_collision_count++;
				return false;
			}// if(calignment)
		}
	}
	//if we made it this far there were no collisions
	//insert the mem into the list.
	if(!mhe.SubMem())
		ExtendMatch(mhe);
	mem_table[bucketI].insert(cur_he, mhe);
	mem_table_count[bucketI]++;
	m_mem_count++;
	return true;
}

void MemHash::EliminateInclusions(){
	//can't remove inclusions from two genomes!
	if(seq_count <= 2)
		return;
	
	uint32 debug_level = 10;
	boolean debug_algorithm = false;
	
	list<MemHashEntry>::iterator mem_table_iter;
	list<MemHashEntry*>::iterator mhe, mem_iter, prev_iter;
	//define an array of lists
	//first index is the number of sequences in the match
	//second index is the sequence that the list is sorted on.
	list<MemHashEntry*>** mhe_list;

	uint32 cur_seq_count;
	uint32 list_count = seq_count - 1;
	mhe_list = new list<MemHashEntry*>*[list_count];
	for(uint32 sdf = 0; sdf < list_count; sdf++)
		mhe_list[sdf] = new list<MemHashEntry*>[list_count];
	MemHashEntry* new_mhe;

	//get all the mems out of the hash table and add them to the appropriate list
	for(uint32 i=0; i < table_size; i++){
		mem_table_iter = mem_table[i].begin();
		for(; mem_table_iter != mem_table[i].end(); mem_table_iter++){
			new_mhe = mem_table_iter->Clone();
			cur_seq_count = new_mhe->SequenceCount() - 2;
			mhe_list[cur_seq_count][0].push_back(new_mhe);
		}
		mem_table[i].clear();
	}

	//keep a vector of iterators with the current positions to check in
	//each level of sorted match lists
	vector<list<MemHashEntry*>::iterator>* start_iters;
	vector<list<MemHashEntry*>::iterator>* cur_starts;
	start_iters = new vector<list<MemHashEntry*>::iterator>[list_count];
	cur_starts = new vector<list<MemHashEntry*>::iterator>[list_count];

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
			uint32 cur_minI = ((list_count - checkI - 1) * (list_count - checkI))/2 + (list_count - checkI - 1);
			uint32 cur_maxI = ((list_count - checkI) * (list_count - checkI + 1))/2 + (list_count - checkI);

			boolean erased_start = false;
			mem_iter = mhe_list[listI][0].begin();
			uint32 cur_start_level;
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
				list<MemHashEntry*> intersect_list;
				for(uint32 cur_startI = cur_minI; cur_startI < cur_maxI - 1; cur_startI++){
					//scan until cur_startI is bigger.
					
					if(debug_algorithm)
						if(cur_starts[cur_start_level][cur_startI] != cur_starts[cur_start_level][cur_startI+1])
							cout << "Against " << **cur_starts[cur_start_level][cur_startI] << "\n";
					while(cur_starts[cur_start_level][cur_startI] != cur_starts[cur_start_level][cur_startI+1] &&
					      MemHashEntry::start_compare(**cur_starts[cur_start_level][cur_startI], **mem_iter) < 0){
						cur_starts[cur_start_level][cur_startI]++;
						if(debug_algorithm)
							if(cur_starts[cur_start_level][cur_startI] != cur_starts[cur_start_level][cur_startI+1])
								cout << "Against " << **cur_starts[cur_start_level][cur_startI] << "\n";
					}
					if(cur_starts[cur_start_level][cur_startI] == cur_starts[cur_start_level][cur_startI+1])
						continue;
					//now check for intersections
					if((*mem_iter)->Intersects(**cur_starts[cur_start_level][cur_startI])){
						//add the mem above it to the intersect list
//						intersect_list.push_back(**cur_starts[cur_start_level][cur_startI]);
						
						//delete the inclusion and insert fragments
						MemHashEntry unique_mhe;
						list<MemHashEntry*>::iterator tmp_iter = mem_iter;
						//insert a unique beginning fragment, if it exists
						if((*mem_iter)->GetUniqueStart(**cur_starts[cur_start_level][cur_startI], unique_mhe)){
							mhe_list[listI][0].insert(mem_iter, unique_mhe.Clone());
							if(debug_algorithm)
								cout << "Added unique start " << unique_mhe << "\n";
						}
						//now insert the end fragment
						if((*mem_iter)->GetUniqueEnd(**cur_starts[cur_start_level][cur_startI], unique_mhe)){
							//scan forward to find the insert point for the end mem
							while(tmp_iter != mhe_list[listI][0].end() && 
							      (*tmp_iter)->FirstStart() == unique_mhe.FirstStart() &&
							      MemHashEntry::start_compare(**tmp_iter, unique_mhe) < 0){
								tmp_iter++;
							}
							mhe_list[listI][0].insert(tmp_iter, unique_mhe.Clone());
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
	
	//now hash all the mems again
	for(uint32 i=0; i < list_count; i++){
		mhe = mhe_list[i][0].begin();
		for(; mhe != mhe_list[i][0].end(); mhe++){
			(*mhe)->SetSubMem(true);
			if(!AddHashEntry(**mhe))
				cout << "You got problems bustah!\n";
			delete *mhe;
		}
	}

	for(uint32 fdfd = 0; fdfd < list_count; fdfd++)
		delete[] mhe_list[fdfd];
	delete[] mhe_list;
}
// Question:  Do I try to make mims at this point while I have the contextual information
// of overlapping mems?
// Answer: Probably
// Question:  How do I store these MiMs?
// Answer: Keep an ordered list of mems and add alignment information in the non-matching
// regions later.
// Question: how do i deal with single genome translocations?
// Answer: you don't have to do that here..??
//void MemHash::GetUniquePieces(list<MemHashEntry*>& intersect_list, MemHashEntry* target_mem, list<MemHashEntry*> unique_list){
	
//}

#include "MemSubsets.h"
#include "MemHashEntry.h"
#include <list>
#include <set>
#include <map>

MemSubsets::MemSubsets( const MemSubsets& ms ){
	debug_algorithm = ms.debug_algorithm;
	mhe = ms.mhe;
	mem_iter = ms.mem_iter;
	prev_iter = ms.prev_iter;
	mhe_list = ms.mhe_list;
	list_count = ms.list_count;
	new_mie = ms.new_mie;
	start_iters = ms.start_iters;
	cur_starts = ms.cur_starts;
	listI = ms.listI;
	cur_start_iter = ms.cur_start_iter;
	next_start_iter = ms.next_start_iter;
	matching_iter = ms.matching_iter;
	erased_start = ms.erased_start;
}

void MemSubsets::Subsets( list<MemHashEntry*>& mem_list ){
	//can't remove inclusions from one mem or two genomes!
	if( mem_list.size() <= 1 || (*mem_list.begin())->SeqCount() <= 2 )
		return;
	
	uint32 debug_level = 100;
	debug_algorithm = false;
	
	list<MemHashEntry*>::iterator mem_list_iter;

	list_count = (*mem_list.begin())->SeqCount() - 1;
	mhe_list = new list<MemHashEntry*>*[list_count];
	for(uint32 sdf = 0; sdf < list_count; sdf++)
		mhe_list[sdf] = new list<MemHashEntry*>[list_count];

	//get all the mems out of the mem_list and add them to the list which
	//corresponds to their match count in the first sort level.
	mem_list_iter = mem_list.begin();
	for(; mem_list_iter != mem_list.end(); mem_list_iter++){
		mhe_list[ (*mem_list_iter)->Multiplicity()-2 ][ 0 ].push_back( *mem_list_iter );
	}
	mem_list.clear();
	
	//first index is the sequence that mems are sorted on
	//second index is a pointer to each successive match level's start level in that list.
	//got it?  good.
	start_iters = new vector<list<MemHashEntry*>::iterator>[list_count];
	cur_starts = new vector<list<MemHashEntry*>::iterator>[list_count];
	
	//sort the highest match level on each sequence.
	mhe_list[ list_count - 1 ][ 0 ].sort(&MemHashEntry::start_lessthan_ptr);
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
	for(listI = list_count - 2; listI >= 0; listI--){
		//sort the current level
		MemHashEntry::seq_compare_start = 0;
		mhe_list[listI][0].sort(&MemHashEntry::start_lessthan_ptr);

		//find inclusions starting with the highest sequence match level
		for(int32 checkI = list_count - 1; checkI > listI; checkI--){
			//the following two vars are the ranges of start iterators to check against in this match level
			uint32 cur_minI = quadratic_li(list_count - checkI - 1) + (list_count - checkI - 1);
			uint32 cur_maxI = quadratic_li(list_count - checkI) + (list_count - checkI);

			erased_start = false;
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
				//note: we don't always need to scan through this entire loop.  It should be fixed such that
				//each if FirstStart() of the cur_start is > FirstStart() of mem_iter we skip the loop.
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
							// they evenly overlap.  
							// handle it appropriately
							HandleOverlap( *matching_iter, *mem_iter );
						}
						matching_iter++;
					}
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

	//now add all the mems to the mem_list again
	for(uint32 i=0; i < list_count; i++){
		mhe = mhe_list[i][0].begin();
		for(; mhe != mhe_list[i][0].end(); mhe++){
			(*mhe)->SetExtended(true);
			mem_list.push_back(*mhe);
		}
	}

// time for the cleanup crew
	for(uint32 fdfd = 0; fdfd < list_count; fdfd++)
		delete[] mhe_list[fdfd];
	delete[] mhe_list;

}


void MemSubsets::EliminateLinkedInclusions( list<MemHashEntry*>& mem_list ) {
	if( mem_list.size() < 2 )
		return;	// can't remove inclusions on less than 2 matches!

	// first separate them into lists based on multiplicity
	// then scan through each multiplicity level, starting with 3, eliminating inclusions.
	
	unsigned int seq_count;
	list<MemHashEntry*>::iterator mhe_iter = mem_list.begin();
	seq_count = (*mhe_iter)->SeqCount();
	
	//can't remove inclusions from two genomes!
	if(seq_count <= 2)
		return;

	// create seq_count - 2 lists, each of which correspond to matches with multiplicity >= 3
	set<MemHashEntry*>* mhe_list = new set<MemHashEntry*>[seq_count + 1];

	//get all the mems out of the hash table and add them to the list which
	//corresponds to their multiplicity.
	list<MemHashEntry*>::iterator mem_table_iter = mem_list.begin();
	for(; mem_table_iter != mem_list.end(); mem_table_iter++){
		unsigned int cou = (*mem_table_iter)->Multiplicity();
		mhe_list[cou].insert( *mem_table_iter );
	}

	// now remove all inclusions, starting with matches of multiplicity 3 and increasing.
	for( uint32 multI = 3; multI <= seq_count; multI++ ) {
		set<MemHashEntry*>::iterator mhe_iter = mhe_list[ multI ].begin();
		for( ; mhe_iter != mhe_list[ multI ].end(); mhe_iter++ ){
			MemHashEntry* cur_mem = *mhe_iter;
			set<MemHashEntry*> tmp_subsets = cur_mem->GetSubsets();
			set<MemHashEntry*>::iterator sub_iter = tmp_subsets.begin();
			MemHashEntry* unique_part = new MemHashEntry( **mhe_iter );
			for(; sub_iter != tmp_subsets.end(); sub_iter++ ) {
				MemHashEntry* submem = *sub_iter;
				if( submem->GetUniqueStart( *cur_mem, *unique_part ) && unique_part->Length() > 0 ){
					// have the new submem check whether it should still be linked to all of its supersets
					// have it tell its supersets to link to it if the link is still good.
					unique_part->UpdateSupersets();
					mhe_list[ unique_part->Multiplicity() ].insert( unique_part );
					
					//allocate a new unique part for next time
					unique_part = new MemHashEntry( *submem );
				}

				if( submem->GetUniqueEnd( *cur_mem, *unique_part ) && unique_part->Length() > 0 ){
					// have the new submem check whether it should still be linked to all of its supersets
					// have it tell its supersets to link to it if the link is still good.
					unique_part->UpdateSupersets();
					mhe_list[ unique_part->Multiplicity() ].insert( unique_part );
					//allocate a new unique part for next time
					unique_part = new MemHashEntry( *submem );
				}

				// delete submem from everywhere
				mhe_list[ submem->Multiplicity() ].erase( submem );
				submem->UnlinkSelf();
			}

			delete unique_part;
		}
	}
	
	// finally copy all the new mems to mem_list.
	mem_list.clear();
	for( unsigned int seqI = 0; seqI <= seq_count; seqI++ )
		mem_list.insert( mem_list.end(), mhe_list[ seqI ].begin(), mhe_list[ seqI ].end() );

}

void MemSubsets::UnlinkMatch( Match* match, map<MatchID_t, Match*> match_id_map ) {
	// unlink from the supersets
	set<MatchID_t>& supersets = match->Supersets();
	set<MatchID_t>::iterator match_iter = supersets.begin();
	for(; match_iter != supersets.end(); match_iter++ ){
		map<MatchID_t, Match*>::iterator id_iter = match_id_map.find( *match_iter );
		id_iter->second->UnlinkSubset( match->MatchID() );
	}
	supersets.clear();

	// unlink from the subsets
	set<MatchID_t>& subsets = match->Subsets();
	match_iter = subsets.begin();
	for(; match_iter != subsets.end(); match_iter++ ){
		map<MatchID_t, Match*>::iterator id_iter = match_id_map.find( *match_iter );
		id_iter->second->UnlinkSuperset( match->MatchID() );
	}
	subsets.clear();

}

void MemSubsets::UpdateSupersets( Match* match, map<MatchID_t, Match*> match_id_map ) {
	set<MatchID_t>& supersets = match->Supersets();
	set<MatchID_t>::iterator super_iter = supersets.begin();
	
	while(super_iter != supersets.end()) {
		map<MatchID_t, Match*>::iterator id_iter = match_id_map.find( *super_iter );
		Match* supermem = id_iter->second;
		if( match->EvenlyOverlaps( *(id_iter->second) ) ){
			supermem->AddSubset( match->MatchID() );
			super_iter++;
		}else {
			set<MatchID_t>::iterator to_del = super_iter;
			super_iter++;
			supersets.erase( to_del );
		}
	}
}


/**
 *  While this code appears to function correctly, it is so terribly slow that it's completely
 *  useless.
 */
void MemSubsets::EliminateLinkedInclusions( list<Match*>& mem_list ) {
	// create a matchID database
	map<MatchID_t, Match*> match_id_map;
	
	// first separate them into lists based on multiplicity
	// then scan through each multiplicity level, starting with 3, eliminating inclusions.
	
	unsigned int seq_count;
	list<Match*>::iterator mhe_iter = mem_list.begin();
	if( mem_list.size() > 0 )
		seq_count = (*mhe_iter)->SeqCount();
	
	//can't remove inclusions from two genomes!
	if(seq_count <= 2)
		return;

	// create seq_count - 2 lists, each of which correspond to matches with multiplicity >= 3
	set<Match*>* mhe_list = new set<Match*>[seq_count + 1];

	//get all the mems out of the hash table and add them to the list which
	//corresponds to their multiplicity.
	list<Match*>::iterator mem_table_iter = mem_list.begin();
	for(; mem_table_iter != mem_list.end(); mem_table_iter++){
		unsigned int cou = (*mem_table_iter)->Multiplicity();
		mhe_list[cou].insert( *mem_table_iter );
		match_id_map.insert( map<MatchID_t, Match*>::value_type( (*mem_table_iter)->MatchID(), *mem_table_iter ) );
	}

	// now remove all inclusions, starting with matches of multiplicity 3 and increasing.
	for( uint32 multI = 3; multI <= seq_count; multI++ ) {
		set<Match*>::iterator mhe_iter = mhe_list[ multI ].begin();
		for( ; mhe_iter != mhe_list[ multI ].end(); mhe_iter++ ){
			Match* cur_mem = *mhe_iter;
			set<MatchID_t> tmp_subsets = cur_mem->Subsets();
			set<MatchID_t>::iterator sub_iter = tmp_subsets.begin();
			Match* unique_part = new Match( **mhe_iter );
			for(; sub_iter != tmp_subsets.end(); sub_iter++ ) {
				map<MatchID_t, Match*>::iterator sub_id_iter = match_id_map.find( *sub_iter );
				Match* submem = sub_id_iter->second;
				if( submem->GetUniqueStart( *cur_mem, *unique_part ) && unique_part->Length() > 0 ){
					// have the new submem check whether it should still be linked to all of its supersets
					// have it tell its supersets to link to it if the link is still good.
					unique_part->UseNewMatchID();
					match_id_map.insert( map<MatchID_t, Match*>::value_type( unique_part->MatchID(), unique_part ) );
					cout << "Adding: " << *unique_part << endl;
					UpdateSupersets( unique_part, match_id_map );
					if( unique_part->Multiplicity() == multI )
						cout << "Error on the bounty\n";
					mhe_list[ unique_part->Multiplicity() ].insert( unique_part );
					
					//allocate a new unique part for next time
					unique_part = new Match( *submem );
				}

				if( submem->GetUniqueEnd( *cur_mem, *unique_part ) && unique_part->Length() > 0 ){
					// have the new submem check whether it should still be linked to all of its supersets
					// have it tell its supersets to link to it if the link is still good.
					unique_part->UseNewMatchID();
					cout << "Adding: " << *unique_part << endl;
					match_id_map.insert( map<MatchID_t, Match*>::value_type( unique_part->MatchID(), unique_part ) );
					UpdateSupersets( unique_part, match_id_map );
					if( unique_part->Multiplicity() == multI )
						cout << "Error on the bounty\n";
					mhe_list[ unique_part->Multiplicity() ].insert( unique_part );
					//allocate a new unique part for next time
					unique_part = new Match( *submem );
				}

				// delete submem from everywhere
				match_id_map.erase( submem->MatchID() );
				mhe_list[ submem->Multiplicity() ].erase( submem );
				cout << "Deleted " << *submem << endl;
				
				UnlinkMatch( submem, match_id_map );
			}

			delete unique_part;
		}
	}
	
	// finally copy all the new mems to mem_list.
	mem_list.clear();
	for( unsigned int seqI = 0; seqI <= seq_count; seqI++ )
		mem_list.insert( mem_list.end(), mhe_list[ seqI ].begin(), mhe_list[ seqI ].end() );

}



void SubsetLinker::HandleOverlap( MemHashEntry* m1, MemHashEntry* m2 ){
//	link them to each other.
	cout << *m1 << endl;
	cout << *m2 << endl;
	m1->LinkSubset(m2);
}

void SubsetInclusionRemover::HandleOverlap( MemHashEntry* m1, MemHashEntry* m2 ){
	//delete the inclusion and insert fragments
	MemHashEntry unique_mhe, *new_mie;
	list<MemHashEntry*>::iterator tmp_iter = mem_iter;
	MemHashEntry bad_mie, debug_mie = **cur_start_iter;
	bad_mie = **mem_iter;
	//insert a unique beginning fragment, if it exists
	if((*mem_iter)->GetUniqueStart(**cur_start_iter, unique_mhe)){
		new_mie = unique_mhe.Clone();
		mhe_list[listI][0].insert(mem_iter, new_mie);
		if(debug_algorithm)
			cout << "Added unique start " << unique_mhe << "\n";
	}

	debug_mie = **cur_start_iter;
	bad_mie = **mem_iter;
	
	//now insert the end fragment
	if((*mem_iter)->GetUniqueEnd(**cur_start_iter, unique_mhe)){
		new_mie = unique_mhe.Clone();

		//scan forward to find the insert point for the end mem
		while(tmp_iter != mhe_list[listI][0].end() && 
		      (*tmp_iter)->FirstStart() == unique_mhe.FirstStart() &&
		      MemHashEntry::start_compare(**tmp_iter, unique_mhe) < 0){
			tmp_iter++;
		}
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

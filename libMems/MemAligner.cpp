#include "MemAligner.h"
#include "MimHash.h"
#include <fstream>
#include "math.h"

MemAligner::MemAligner(uint32 seqcount){
	seq_count = seqcount;
}

void MemAligner::CreateAlignment(list<MimHashEntry*>* mem_list){
	//convert the mem list to a set of the n-way matches.
	mem_size_set.insert(mem_list[seq_count-2].begin(), mem_list[seq_count-2].end());
	cout << "Creating ordered sets for " << mem_size_set.size() << " mems.\n";

	current_alignment = 0;
	CreateOrderedSet();
	ofstream ordered_set_file("ordered_sets.txt");
	ordered_set_starts.push_back(ordered_supersets.size());
	for(uint32 setI = 0; setI < ordered_set_starts.size()-1; setI++){
		cout << "There are " << ordered_set_starts[setI+1] - ordered_set_starts[setI] << " elements in set " << setI << ".\n"; 
		ordered_set_file << ">> Set " << setI << endl;
		for(uint32 memI = ordered_set_starts[setI]; memI < ordered_set_starts[setI+1]; memI++)
			ordered_set_file << *(ordered_supersets[memI]) << endl;
	}
	ordered_set_file.close();
	
}

float64 MemAligner::DistanceScore(MimHashEntry* m1, MimHashEntry* m2){
	// how da funk is dis going to werk?
	// try the mems on the end first, I guess, working inward.
	// must also ensure that none of the highest match level matches cross each other.
	// 
}

//returns positive if m2 adds new sequence to the end of m1 in each sequence
//returns negative if m2 adds new sequence to the beginning of m1 in each sequence
//returns 0 otherwise
//if the return is positive or negative then it is the euclidean distance between the
//end of one mem and the start of the next.
float64 MemAligner::DistanceScore(MemHashEntry* m1, MemHashEntry* m2){
	uint64 squares_sum = 0;
	uint32 first_start = m1->FirstStart();
	if(first_start != m2->FirstStart())
		return 0;

	//use FirstStart here.
	//if the first one starts later then reverse their order
	if((*m1)[first_start] > (*m2)[first_start])
		return -DistanceScore(m2, m1);

	int64 diff;
	int64 end_diff;
	int64 sub_len;
	float64 score;
	for(uint32 seqI = 0; seqI < m1->SeqCount(); seqI++){
		int64 start_1 = (*m1)[seqI];
		int64 start_2 = (*m2)[seqI];
		//check for consistent increases or decreases across genomes
		//and for similar match direction
		if((start_2 < 0 && start_1 > 0) ||
			(start_2 > 0 && start_1 < 0) ||
			(start_2 == 0 && start_1 != 0) ||
			(start_2 != 0 && start_1 == 0))
			return 0;
		
		diff = start_2 - start_1;

		end_diff = (int64)m2->Length() - (int64)m1->Length();
		sub_len = -(int64)m1->Length();
		
		if(start_2 < 0){
			end_diff = -end_diff;
			sub_len = -(int64)m1->Length();
		}

		end_diff += diff;

		if(diff < 0)
			return 0;
		if(end_diff <= 0)
			return 0;
		squares_sum += (uint64) pow((end_diff + sub_len), 2);
	}
	score = sqrt(squares_sum);
	return score;
}

//This is slow but it works...
void MemAligner::CreateOrderedSet(){
	set<MimHashEntry*> mem_search_set, mem_reverse_set;
	boolean debug_dynprog = false;
	// do the dynamic programming technique thingee.
	while(mem_size_set.size() > 0){
		MimHashEntry* biggest_mem = *mem_size_set.begin();
		MimHashEntry* best_score_mem, *current_mem = biggest_mem;
		multiset<MimHashEntry*, MimLengthCompare>::iterator best_score_mem_iter;
		int64 best_score = 0;  // 0 is invalid, anything else is better
		mem_size_set.erase(mem_size_set.begin());
		cout << mem_size_set.size() << " mems in mem size set\n";
		// score_mem will be the current mem being scored
		set<MimHashEntry*>::iterator score_mem, to_erase;

		// going_forward will be true while we are extending the set forward and false otherwise
		boolean going_forward = true;
		//true during the first round of scoring on the biggest mem only.
		boolean first_scoring = true;

		//start a new ordered superset
		ordered_set_starts.push_back(ordered_supersets.size());
		ordered_supersets.push_back(biggest_mem);
		
		//initialize the search set
		mem_search_set.clear();
		mem_reverse_set.clear();
		cout << mem_search_set.size() << " mems in mem search set before clear\n";
		mem_search_set.insert(mem_size_set.begin(), mem_size_set.end());
		cout << mem_search_set.size() << " mems in mem search set\n";
		
		//find closest mems in the forward direction first, then in backward direction
		while(true){
			score_mem = mem_search_set.begin();
			best_score = 0;
			if(debug_dynprog)
				cout << "Current mem: " << *current_mem << endl;
			for(; score_mem != mem_search_set.end(); score_mem++){
				if(debug_dynprog)
					cout << "score mem: " << **score_mem;
				int64 cur_score = DistanceScore(current_mem, *score_mem);
				if(debug_dynprog)
					cout << "\t" << cur_score << endl;
				if(going_forward && cur_score <= 0 ||
				  !going_forward && cur_score >= 0){
				  	if(score_mem != mem_search_set.begin()){
				  		if(first_scoring)
					  		mem_reverse_set.insert(*score_mem);
					  	to_erase = score_mem--;
						mem_search_set.erase(to_erase);
					}
					continue;
				}
				if((going_forward && cur_score < best_score) ||
				  (!going_forward && cur_score > best_score) ||
				  (best_score == 0)){
				  	best_score_mem = *score_mem;
					best_score = cur_score;
				}
			}

			if(best_score != 0){
//				for(uint32 memI = 0; memI < ordered_supersets.size(); memI++){
//					if(ordered_supersets[memI] == best_score_mem)
//						cout << "Already in set\n";
//				}
				//add the best scoring mem
				if(going_forward)
					ordered_supersets.push_back(best_score_mem);
				else
					ordered_supersets.insert(ordered_supersets.begin() + ordered_set_starts[ordered_set_starts.size()-1], best_score_mem);
				//set it to be the next examined mem
				current_mem = best_score_mem;
				//remove it from the search set
				mem_search_set.erase(current_mem);

				//delete it from the remaining mem set
				best_score_mem_iter = mem_size_set.find(best_score_mem);
				if(best_score_mem != *best_score_mem_iter){
					cout << "Bad search\n";
				}
				if(best_score_mem_iter == mem_size_set.end()){
					cout << "Error wilson\n";
				}else
					mem_size_set.erase(best_score_mem_iter);
			}else if(going_forward){
				//now extend in the other direction
				going_forward = false;
				current_mem = biggest_mem;
				mem_search_set = mem_reverse_set;
			}else
				//return;
				break;
		}
	}
}

/*
void MemAligner::CreateCoreSets(){
	//create MimHashEntries by stringing together core sets.
	//use a recursive technique which is bounded by the next mem in the ordered set.
	for(uint32 setI = 0; setI < ordered_set_starts.size()-1; setI++){
		uint32 memI = ordered_set_starts[setI];
		for(; memI < ordered_set_starts[setI+1]; memI++){
			
		}
	}
}

void MemAligner::LinkForward(MemInclusionEntry* current_mem, MemInclusionEntry* start_bound, MemInclusionEntry* end_bound, MimHashEntry& core_mim){
	//check each of current_mem's overlaps to see if they are consistent with bounding_mem
}
*/

void MemAligner::CreateSet(MimHashEntry* first_mem){
/*	set<MemInclusionEntry*> to_examine;
	set<MemInclusionEntry*>::iterator to_del;
	MemInclusionEntry* next_mem;
	to_examine.insert(first_mem);
	MemInclusionEntry* current_mem;

	//do a depth first search to extract the set
	//first the underlaps
	while(to_examine.size() > 0){
		current_mem = *(to_examine.begin());
		to_examine.erase(to_examine.begin());
		alignment_group.insert(current_mem);
		mem_size_set.erase(current_mem);

		while(current_mem->m_underlaps.size() != 0){
			next_mem = *(current_mem->m_underlaps.begin());
			current_mem->m_underlaps.erase(current_mem->m_underlaps.begin());
			
			//delete this one from next_mem to prevent it from
			//going in a loop
			to_del = next_mem->m_overlaps.find(current_mem);
			next_mem->m_overlaps.erase(to_del);
			
			//add it to the examination list if it hasnt already been seen
			if(alignment_group.find(next_mem) == alignment_group.end())
				to_examine.insert(next_mem);
		}

		//then the overlaps
		while(current_mem->m_overlaps.size() != 0){
			next_mem = *(current_mem->m_overlaps.begin());
			current_mem->m_overlaps.erase(current_mem->m_overlaps.begin());
			
			//delete this one from next_mem to prevent it from
			//going in a loop
			to_del = next_mem->m_underlaps.find(current_mem);
			next_mem->m_underlaps.erase(to_del);
			
			//add it to the examination list if it hasnt already been seen
			if(alignment_group.find(next_mem) == alignment_group.end())
				to_examine.insert(next_mem);
		}
	
		//finally, add this mem to the alignment group and delete it from mem_size_set.
//		alignment_group.push_back(current_mem);

		//also, increment the size counters.
		group_mem_count[current_alignment]++;
		group_coverage_count[current_alignment] += current_mem->Length();
	}
*/
}

//head towards the beginning of the sequence
/*void MemAligner::ExamineUnderlaps(MemInclusionEntry* this_mem, MemInclusionEntry* parent){
	//first thing to do is check if there are any parents which could
	//overlap this mem.  if not, it gets added to the alignment
	if(this_mem->m_underlaps.size() == 0){
		//add this mem to the sequences
		
	}else{
		//find the parent's nearest neighbor to the left and jump to it.
		set<MemInclusionEntry*>::iterator neighbor_iter;
		neighbor_iter = this_mem->m_underlaps.find(parent);
		if(neighbor_iter != this_mem->m_underlaps.begin()){
			
		}
		//how far away is the neighbor?
		if(neighbor_iter->Offset()
	}
}
*/
void MemAligner::ExamineOverlaps(MimHashEntry* child){

}

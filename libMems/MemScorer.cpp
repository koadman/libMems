#include "MemScorer.h"
#include "math.h"

MemScorer::MemScorer(uint32 hit_resolution){
	m_hit_resolution = hit_resolution;
	percent_hits = MATCH_NOT_SCORED;
	subseq_identity = MATCH_NOT_SCORED;
	heuristic_hits = MATCH_NOT_SCORED;
	hit_matrix = NULL;
	identity_matrix = NULL;
}

MemScorer::~MemScorer(){
	if(hit_matrix != NULL){
		for(uint32 i = 0; i < seq_count; i++){
			delete[] hit_matrix[i];
		}
		delete[] hit_matrix;
	}
	if(identity_matrix != NULL){
		for(uint32 i = 0; i < seq_count; i++){
			delete[] identity_matrix[i];
		}
		delete[] identity_matrix;
	}
}

MemScorer* MemScorer::Clone() const{

}

float64 MemScorer::PercentHitsScore(){
	if(percent_hits != MATCH_NOT_SCORED)
		return percent_hits;
	
	uint64 total_length = 0;
	float64 current_percent = 0;
	AllocateMatrix();
	FindMatches(sar_table[0]->GetMerMaskSize());

	for(uint32 i = 0; i < seq_count; i++){
		total_length += sar_table[i]->Length();
		for(uint32 j = i + 1; j < seq_count; j++){
			percent_hits += hit_matrix[i][j];
		}
	}
	percent_hits /= (float64) total_length;
	return percent_hits;
}

float64 MemScorer::SubsequenceIdentityScore(){
	if(subseq_identity != MATCH_NOT_SCORED)
		return subseq_identity;
	
	return 0;
}

//semi-random sampling with replacement
float64 MemScorer::HeuristicHitsScore(){
	if(heuristic_hits != MATCH_NOT_SCORED)
		return heuristic_hits;

	float64 sample_percent = 0.01;
	uint32 sample_size_min = 100;
	uint32 sample_size_range = 1000;
	
	uint32 current_total = 0;
	float64 current_percent = 0;
	AllocateMatrix();
	while(current_percent < sample_percent){
		uint32 current_size = rand() % sample_size_range;
		current_size += sample_size_min;
		uint32 current_pos = rand() % (sar_table[0]->Length() - current_size);
		FindMatches(sar_table[0]->GetMerMaskSize(), current_pos, current_size);
		current_total += current_size;
		current_percent += (float64)current_total / (float64)sar_table[0]->Length();
	}
	for(uint32 i = 0; i < seq_count; i++){
		for(uint32 j = i + 1; j < seq_count; j++){
			heuristic_hits += hit_matrix[i][j];
		}
	}
	heuristic_hits /= (float64) current_total;
	return heuristic_hits;
}

void MemScorer::AllocateMatrix(){

	if(hit_matrix != NULL){
		for(uint32 i = 0; i < seq_count; i++){
			delete[] hit_matrix[i];
		}
		delete[] hit_matrix;
	}
	if(identity_matrix != NULL){
		for(uint32 i = 0; i < seq_count; i++){
			delete[] identity_matrix[i];
		}
		delete[] identity_matrix;
	}
	hit_matrix = new uint32*[seq_count];
	for(uint32 k = 0; k < seq_count; k++){
		hit_matrix[k] = new uint32[seq_count];
		for(uint32 m = 0; m < seq_count; m++)
			hit_matrix[k][m] = 0;
	}
	identity_matrix = new uint32*[seq_count];
	for(uint32 l = 0; l < seq_count; l++){
		identity_matrix[l] = new uint32[m_hit_resolution];
		for(uint32 n = 0; n < m_hit_resolution; n++)
			identity_matrix[l][n] = 0;
	}
}

boolean MemScorer::EnumerateMatches(list<idmer>& match_list){
	//this must call HashMatch on every possible combination of matches in the list.
	
	match_list.sort(&idmer_id_lessthan);
	uint32 id_count = 0;
	vector<uint32> id_start;
	vector<list<idmer>::iterator> id_pos;
	vector<list<idmer>::iterator> id_end;
	list<idmer>::iterator iter = match_list.begin();
	list<idmer>::iterator iter2 = match_list.begin();
	iter2++;
	id_start.push_back(0);
	id_pos.push_back(iter);
	uint32 prev_i = 0;
	for(uint32 i=0; iter2 != match_list.end(); i++){
		if(iter->id != iter2->id){
			uint32 identity_hits = i - prev_i;
			if(identity_hits >= m_hit_resolution)
				identity_hits = m_hit_resolution - 1;
			identity_matrix[sarid_table[iter->id]][identity_hits]++;
			id_start.push_back(i);
			id_pos.push_back(iter2);
		}
		iter++;
		iter2++;
	}
	
	//the following loop iterates through all possible combinations of idmers with
	//different id's and hashes them.
	id_end = id_pos;
	id_end.push_back(match_list.end());
	while(true){
		list<idmer> cur_match;
		for(uint32 k = 0; k < id_pos.size(); k++){
			cur_match.push_back(*id_pos[k]);
		}
		if(cur_match.size() > seq_count)
			cout << "dshgf";
		HashMatch(cur_match);
		cur_match.clear();

		//increment the iterators (like an odometer)
		uint32 m = id_pos.size() - 1;
		while(true){
			id_pos[m]++;
			if(id_pos[m] == id_end[m+1]){
				if(m == 0)
					return true;
				id_pos[m] = id_end[m];
				m--;
			}else
				break;
		}
	}

	return true;
}

//n^2 across the number of genomes
boolean MemScorer::HashMatch(list<idmer>& match_list){
	
	list<idmer>::iterator firstI = match_list.begin();
	list<idmer>::iterator secondI = match_list.begin();
	secondI++;
	for(; firstI != match_list.end(); firstI++){
		secondI = firstI;
		secondI++;
		for(; secondI != match_list.end(); secondI++)
			hit_matrix[sarid_table[firstI->id]][sarid_table[firstI->id]]++;
	}
}

#include "MemScorer.h"
#include <cmath>

MemScorer::MemScorer(uint32 hit_resolution){
	m_hit_resolution = hit_resolution;
	percent_hits = MATCH_NOT_SCORED;
	subseq_identity = MATCH_NOT_SCORED;
	heuristic_hits = MATCH_NOT_SCORED;
	hit_matrix = NULL;
	identity_matrix = NULL;
	detail_list = NULL;
}

MemScorer::MemScorer( const MemScorer& ms ) : MemHash(ms) {
	m_hit_resolution = ms.m_hit_resolution;
	percent_hits = ms.percent_hits;
	subseq_identity = ms.subseq_identity;
	heuristic_hits = ms.heuristic_hits;
	hit_matrix = ms.hit_matrix;
	identity_matrix = ms.identity_matrix;
	detail_list = ms.detail_list;
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
	if(detail_list != NULL)
		delete[] detail_list;
}

MemScorer* MemScorer::Clone() const{
	return new MemScorer( *this );
}

float64 MemScorer::PercentHitsScore(){
	if(percent_hits != MATCH_NOT_SCORED)
		return percent_hits;
	
	uint64 total_length = 0;
	AllocateMatrix();
	MatchFinder::FindMatches();

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

	AllocateMatrix();
	float64 sample_percent = 0.01;
//	uint32 sample_size_min = 1000;
	gnSeqI total_len = 0;
	for(uint32 sarI = 0; sarI < sar_table.size(); sarI++)
		total_len += sar_table[sarI]->Length();
	
	gnSeqI sample_size = (gnSeqI) (sample_percent * total_len);
	vector<gnSeqI> start_points, search_len;
	GetBreakpoint(0, sample_size, search_len);
	for(uint32 seqI = 0; seqI < sar_table.size(); seqI++){
		start_points.push_back(0);
	}
	MatchFinder::FindMatches(start_points, search_len);
	
	//for each mem, count the number of bases covered in each genome.
	for(uint32 bucketI = 0; bucketI < table_size; bucketI++){
		set<MemHashEntry*, MheCompare>::const_iterator mem_iter = mem_table[bucketI].begin();
		for(; mem_iter != mem_table[bucketI].end(); mem_iter++){
			for(uint32 seqI = 0; seqI < seq_count; seqI++){
				if((**mem_iter)[seqI] != MEM_NO_MATCH){
					for(uint32 seq2 = 0; seq2 < seq_count; seq2++){
						if((**mem_iter)[seq2] != MEM_NO_MATCH){
							identity_matrix[seqI][seq2] += (*mem_iter)->Length();
						}
					}
				}
			}
		}
	}

	//now do it the detailed way if there is space.
	if(seq_count <= MAX_DETAIL_SEQS && detail_list != NULL){
		for(uint32 bucketI = 0; bucketI < table_size; bucketI++){
			set<MemHashEntry*, MheCompare>::const_iterator mem_iter = mem_table[bucketI].begin();
			for(; mem_iter != mem_table[bucketI].end(); mem_iter++){
				//get the sequence matching entry into match bits
				uint32 match_bits = 0;
				for(uint32 seqI = 0; seqI < seq_count; seqI++){
					match_bits <<= 1;
					if((**mem_iter)[seqI] != MEM_NO_MATCH)
						match_bits |= 0x1;
				}
				detail_list[match_bits] += (*mem_iter)->Length();
			}
		}
	}
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
//		FindMatches(current_pos, current_size);
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
	if(detail_list != NULL)
		delete[] detail_list;
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
	if(seq_count <= MAX_DETAIL_SEQS){
		uint64 list_size = (uint64) pow(2,(long double)seq_count);
		detail_list = new uint32[list_size];
		for(uint32 i=0; i < list_size; i++)
			detail_list[i] = 0;
	}
}


//n^2 across the number of genomes
boolean MemScorer::HashMatch(list<idmer>& match_list){
	
	list<idmer>::iterator firstI = match_list.begin();
	list<idmer>::iterator secondI;
	for(; firstI != match_list.end(); firstI++){
		secondI = firstI;
		secondI++;
		for(; secondI != match_list.end(); secondI++)
			hit_matrix[firstI->id][secondI->id]++;
	}
	return MemHash::HashMatch(match_list);
}

void MemScorer::PrintIDMatrix(ostream& os){
	if(identity_matrix == NULL)
		return;
	for(uint32 i=0; i < seq_count; i++){
		for(uint32 j=0; j < seq_count; j++){
			os << "\t" << identity_matrix[i][j];
		}
		os << "\n";
	}
}

void MemScorer::PrintHitMatrix(ostream& os){
	if(hit_matrix == NULL)
		return;
	for(uint32 i=0; i < seq_count; i++){
		for(uint32 j=0; j < seq_count; j++){
			os << "\t" << hit_matrix[i][j];
		}
		os << "\n";
	}
}

void MemScorer::PrintDetailList(ostream& os){
	uint64 detailMax = (uint64) pow(2, (long double)seq_count);
	for(uint64 detailI = 0; detailI < detailMax; detailI++){
		string seq_string;
		uint64 curDetail = detailI;
		for(uint32 seqI = 0; seqI < seq_count; seqI++){
			if(curDetail & 0x1)
				seq_string = "1 " + seq_string;
			else
				seq_string = "0 " + seq_string;
			curDetail >>= 1;
		}
		os << seq_string << '\t' << detail_list[detailI] << '\n';
	}
}

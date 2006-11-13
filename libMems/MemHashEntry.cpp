#include "MemHashEntry.h"

uint32 MemHashEntry::seq_compare_start = 0;

boolean MemHashEntry::offset_lessthan(const MemHashEntry& a, const MemHashEntry& b){
	return a.m_offset < b.m_offset;
}

boolean MemHashEntry::start_lessthan_ptr(const MemHashEntry* a, const MemHashEntry* b){
	int32 start_diff = a->FirstStart() - b->FirstStart();
	if(start_diff == 0){
		uint32 m_count = a->m_start.size();
		m_count = m_count <= b->m_start.size() ? m_count : b->m_start.size();
		for(uint32 seqI = seq_compare_start; seqI < m_count; seqI++){
			int64 a_start = a->m_start[seqI], b_start = b->m_start[seqI];
			if(a_start < 0)
				a_start = -a_start + a->m_length - a->m_mersize;
			if(b_start < 0)
				b_start = -b_start + b->m_length - b->m_mersize;
			int64 diff = a_start - b_start;
			if(a_start == MEM_NO_MATCH || b_start == MEM_NO_MATCH)
				continue;
			else if(diff == 0)
				continue;
			else
				return diff < 0;
		}
	}
	return start_diff < 0;
}

boolean MemHashEntry::end_lessthan(const MemHashEntry& a, const MemHashEntry& b){
	int32 start_diff = a.FirstStart() - b.FirstStart();
	int32 length_diff = a.m_length - b.m_length;
	if(start_diff == 0){
		uint32 m_count = a.m_start.size();
		m_count = m_count <= b.m_start.size() ? m_count : b.m_start.size();
		for(uint32 seqI = seq_compare_start; seqI < m_count; seqI++){
			int64 a_start = a.m_start[seqI], b_start = b.m_start[seqI];
			if(a_start < 0)
				a_start = -a_start + a.m_length - a.m_mersize;
			if(b_start < 0)
				b_start = -b_start + b.m_length - b.m_mersize;
			int64 diff = length_diff + a_start - b_start;
			if(a_start == MEM_NO_MATCH || b_start == MEM_NO_MATCH)
				continue;
			else if(diff == 0)
				continue;
			else
				return diff < 0;
		}
	}
	return start_diff < 0;
}

//ignores mem_no_matches
int64 MemHashEntry::start_compare(const MemHashEntry& a, const MemHashEntry& b){
	uint32 m_count = a.m_start.size();
	m_count = m_count <= b.m_start.size() ? m_count : b.m_start.size();
	for(uint32 seqI = 0; seqI < m_count; seqI++){
		int64 a_start = a.m_start[seqI], b_start = b.m_start[seqI];
		if(a_start < 0)
			a_start = -a_start + a.m_length - a.m_mersize;
		if(b_start < 0)
			b_start = -b_start + b.m_length - b.m_mersize;
		int64 diff = a_start - b_start;
		if(a_start == MEM_NO_MATCH || b_start == MEM_NO_MATCH)
			continue;
		else if(diff == 0)
			continue;
		else
			return diff;
	}
	return 0;
}

int MemHashEntry::fuzzy_align(const MemHashEntry& a, const MemHashEntry& b, uint32 tolerance){
 	for(uint32 i = 0; i < a.m_start.size(); i++){
 		int64 cur_end = a.m_start[i];
 		int64 cur_start = b.m_start[i];
 		if(cur_end == MEM_NO_MATCH){
 			if(cur_start == MEM_NO_MATCH)
 				continue;
 			else
 				return -1;
 		}else if(cur_end < 0){
 			if(cur_start < 0){
 				int64 tmp_end = -cur_end;
 				cur_end = -cur_start + b.m_length;
 				cur_start = tmp_end;
 			}else
 				return -1;
 		}else
	 		cur_end += a.m_length;
 		if(cur_end > cur_start)
 			return -1;
 		if(cur_end + tolerance < cur_start)
 		 	return 1;
 	}
 	return 0;

}

MemHashEntry::MemHashEntry(){
	m_offset = 0;
	m_submem = false;
	m_length = 0;
	m_usedSeqs = 0;
	m_firstStart = UINT32_MAX;
}

MemHashEntry::MemHashEntry(uint32 seq_count, const gnSeqI mersize, boolean submem){
	m_offset = 0;
	m_submem = submem;
	m_mersize = mersize;
	m_length = 0;
	m_usedSeqs = 0;
	m_firstStart = UINT32_MAX;
	for(uint32 i=0; i < seq_count; i++)
		m_start.push_back(MEM_NO_MATCH);
}
MemHashEntry::MemHashEntry(const MemHashEntry& mhe){
	m_offset = mhe.m_offset;
	m_submem = mhe.m_submem;
	m_mersize = mhe.m_mersize;
	m_length = mhe.m_length;
	m_start = mhe.m_start;
	m_usedSeqs = mhe.m_usedSeqs;
	m_firstStart = mhe.m_firstStart;
}
MemHashEntry::~MemHashEntry(){

}
MemHashEntry* MemHashEntry::Clone() const{
	return new MemHashEntry(*this);
}
MemHashEntry& MemHashEntry::operator=(const MemHashEntry& mhe){
	m_offset = mhe.m_offset;
	m_submem = mhe.m_submem;
	m_mersize = mhe.m_mersize;
	m_length = mhe.m_length;
	m_start = mhe.m_start;
	m_usedSeqs = mhe.m_usedSeqs;
	m_firstStart = mhe.m_firstStart;
	return *this;
}

uint32 MemHashEntry::FirstStart() const{
//	if(m_firstStart < UINT32_MAX)
//		return m_firstStart;
	for(uint32 startI=0; startI < m_start.size(); startI++)
		if(m_start[startI] != MEM_NO_MATCH){
			return startI;
			break;
		}
	return UINT32_MAX;
}

void MemHashEntry::HomogenizeParity(vector<int64>& start1, vector<int64>& start2, const uint32 startI) const{
	//make sure the parity of each genome is the same when the start level is different 
	if(start2[startI] > 0 && start1[startI] < 0){
		//invert m_start1...
		for(uint32 j=startI; j < start1.size(); j++)
			start1[j] *= -1;
	}else if(start1[startI] > 0 && start2[startI] < 0){
		//invert m_start2...
		for(uint32 j=startI; j < start2.size(); j++)
			start2[j] *= -1;
	}
}

boolean MemHashEntry::Intersects(const MemHashEntry& mhe) const{
	uint32 i;
	int64 diff_i;
	int64 diff;
	vector<int64> m_start1 = m_start, m_start2 = mhe.m_start;
	uint32 seq_count = m_start2.size();
	//find the first matching sequence...
	for(i=0; i < seq_count; i++){
		if(m_start2[i] == MEM_NO_MATCH || m_start1[i] == MEM_NO_MATCH)
			continue;
		else
			break;
	}
	//make sure the parity of each genome is the same
	if(FirstStart() != mhe.FirstStart())
		HomogenizeParity(m_start1, m_start2, i);
	diff = m_start2[i] - m_start1[i];

	//check for intersection properties
	if(diff >= m_length || -diff >= mhe.m_length)
		return false;

	//everything is ok so far, check for alignment
	int64 diff_rc = (int64)mhe.m_length - (int64)m_length + diff;
	for(i++; i < seq_count; i++){
		//check for consistent alignment between all genomes
		//in the case of revcomp, diff_i must equal m.length - diff
		diff_i = m_start2[i] - m_start1[i];
		//it's ok if either mem has a sequence which the other doesn't
		if(m_start2[i] == MEM_NO_MATCH || m_start1[i] == MEM_NO_MATCH)
			continue;
		else if(m_start2[i] < 0 && diff_rc == diff_i)
			continue;
		else if(diff != diff_i )
			return false;
	}
	//it was aligned and intersects.
	return true;
}

boolean MemHashEntry::GetUniqueStart(const MemHashEntry& mhe, MemHashEntry& unique_mhe) const{
	uint32 i;
	int64 diff_i;
	int64 diff;
	uint32 seq_count = mhe.m_start.size();
	vector<int64> m_start1 = m_start, m_start2 = mhe.m_start;

	//find the first matching sequence...
	for(i=0; i < seq_count; i++){
		if(mhe.m_start[i] == MEM_NO_MATCH || m_start[i] == MEM_NO_MATCH)
			continue;
		else
			break;
	}

	//make sure the parity of each genome is the same
	if(FirstStart() != mhe.FirstStart())
		HomogenizeParity(m_start1, m_start2, i);
	diff = m_start2[i] - m_start1[i];

	//check for overlao properties
	if(diff >= m_length || diff <= 0)
		return false;

	//everything is ok so far, check for alignment
	int64 diff_rc = (int64)mhe.m_length - (int64)m_length + diff;
	for(i++; i < seq_count; i++){
		//check for consistent alignment between all genomes
		//in the case of revcomp, diff_i must equal m.length - diff
		diff_i = m_start2[i] - m_start1[i];
		//it's ok if this mem has a sequence which mhe doesn't
		if(m_start2[i] == MEM_NO_MATCH || m_start1[i] == MEM_NO_MATCH)
			continue;
		else if(m_start2[i] < 0 && diff_rc == diff_i)
			continue;
		else if(diff != diff_i )
			return false;
	}
	//it was aligned and intersects.
	unique_mhe = *this;
	unique_mhe.CropEnd(m_length - diff);
	return true;
}
boolean MemHashEntry::GetUniqueEnd(const MemHashEntry& mhe, MemHashEntry& unique_mhe) const{
	uint32 i;
	int64 diff_i;
	int64 diff;
	uint32 seq_count = mhe.m_start.size();
	vector<int64> m_start1 = m_start, m_start2 = mhe.m_start;

	//find the first matching sequence...
	for(i=0; i < seq_count; i++){
		if(mhe.m_start[i] == MEM_NO_MATCH || m_start[i] == MEM_NO_MATCH)
			continue;
		else
			break;
	}
	//make sure the parity of each genome is the same
	if(FirstStart() != mhe.FirstStart())
		HomogenizeParity(m_start1, m_start2, i);
	diff = m_start2[i] - m_start1[i];

	//check for overlap properties
	if(m_length <= mhe.m_length + diff)
		return false;


	//everything is ok so far, check for alignment
	int64 diff_rc = (int64)mhe.m_length - (int64)m_length + diff;
	for(i++; i < seq_count; i++){
		//check for consistent alignment between all genomes
		//in the case of revcomp, diff_i must equal m.length - diff
		diff_i = m_start2[i] - m_start[i];
		//it's ok if this mem has a sequence which mhe doesn't
		if(m_start2[i] == MEM_NO_MATCH || m_start1[i] == MEM_NO_MATCH)
			continue;
		else if(m_start2[i] < 0 && diff_rc == diff_i)
			continue;
		else if(diff != diff_i )
			return false;
	}
	//it was aligned and intersects.
	unique_mhe = *this;
	unique_mhe.CropStart(mhe.m_length + diff);
	return true;
}

void MemHashEntry::CropStart(gnSeqI crop_amount){
	for(uint32 i=0; i < m_start.size(); i++){
		if(m_start[i] != MEM_NO_MATCH && m_start[i] > 0)
			m_start[i] += crop_amount;
	}
	m_length -= crop_amount;
}

void MemHashEntry::CropEnd(gnSeqI crop_amount){
	for(uint32 i=0; i < m_start.size(); i++){
		if(m_start[i] != MEM_NO_MATCH && m_start[i] < 0)
			m_start[i] -= crop_amount;
	}
	m_length -= crop_amount;
}


boolean MemHashEntry::Contains(const MemHashEntry& mhe) const{
	uint32 i;
	int64 diff_i;
	int64 diff;
	uint32 seq_count = mhe.m_start.size();
	//check for a consistent number of genomes and
	//identical generalized offsets
	if(m_start.size() != seq_count || m_offset != mhe.m_offset)
		return false;

	i = mhe.FirstStart();
	for(; i < seq_count; i++){
		diff = mhe.m_start[i] - m_start[i];
		if(mhe.m_start[i] == MEM_NO_MATCH){
			if(diff != 0)
				return false;
		}else
			break;
	}

	//check for containment properties
	if(diff < 0 || m_length < mhe.m_length + diff)
		return false;

	//everything is ok so far, check for alignment
	int64 diff_rc = -(m_length - diff) + m_mersize;
	for(i++; i < seq_count; i++){
		//check for consistent alignment between all genomes
		//in the case of revcomp, diff_i must equal m.length - diff
		diff_i = mhe.m_start[i] - m_start[i];
		//it's ok if this mem has a sequence which mhe doesn't
		if(mhe.m_start[i] == MEM_NO_MATCH)
			continue;
		else if(mhe.m_start[i] < 0 && diff_rc == diff_i)
			continue;
		else if(diff != diff_i )
			return false;
	}
	//it was contained.
	return true;
}

void MemHashEntry::DetailCompare(const MemHashEntry& mhe, int& offset, boolean& aligned, int64& diff){
	uint32 i;
	int64 diff_i;
	uint32 seq_count = mhe.m_start.size();
	offset = m_offset - mhe.m_offset;
	aligned = false;
	if(m_start.size() != seq_count || offset != 0)
		return;

	//offsets are the same, check if the alignment is the same...
	for(i=0; i < seq_count; i++){
		diff = mhe.m_start[i] - m_start[i];
		if(mhe.m_start[i] == MEM_NO_MATCH){
			if(diff != 0)
				return;
		}else
			break;
	}
	//everything is ok so far, check for alignment
	int64 diff_rc = -(m_length - diff) + m_mersize;
	for(i++; i < seq_count; i++){
		//check for consistent alignment between all genomes
		//in the case of revcomp, diff_i must equal m.length - diff
		diff_i = mhe.m_start[i] - m_start[i];
		if(mhe.m_start[i] == MEM_NO_MATCH && diff_i == 0)
			continue;
		else if(mhe.m_start[i] < 0 && diff_rc == diff_i)
			continue;
		else if(diff != diff_i )
			return;
	}
	aligned = true;
}

void MemHashEntry::CalculateOffset(){
	int64 ref_start;
	uint32 seq_count = m_start.size();
	uint32 seqI = 0;
	m_offset = 0;
	for(; seqI < seq_count; seqI++)
		if(m_start[seqI] != MEM_NO_MATCH){
			ref_start = m_start[seqI];
			m_usedSeqs = 1;
			m_firstStart = seqI;
			break;
		}
	for(seqI++; seqI < seq_count; seqI++)
		if(m_start[seqI] != MEM_NO_MATCH){
			m_offset += m_start[seqI] - ref_start;
			m_usedSeqs++;
		}
}

void MemHashEntry::ExtendStart(gnSeqI extend_amount){
	m_length += extend_amount;
	uint32 seq_count = m_start.size();
	for(uint32 seqI = 0; seqI < seq_count; seqI++)
		if(m_start[seqI] > 0)
			m_start[seqI] -= extend_amount;
}

void MemHashEntry::ExtendEnd(gnSeqI extend_amount){
	m_length += extend_amount;
	uint32 seq_count = m_start.size();
	for(uint32 seqI = 0; seqI < seq_count; seqI++)
		if(m_start[seqI] < 0)
			m_start[seqI] += extend_amount;
}

ostream& operator<<(ostream& os, const MemHashEntry& mhe){ //write to stream.
	os << mhe.m_length;
	for(uint32 i=0; i < mhe.m_start.size(); i++)
		os << '\t' << mhe.m_start[i];
	return os;
}

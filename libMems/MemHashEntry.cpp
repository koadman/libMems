#include "MemHashEntry.h"
#include "gn/gnException.h"
#include "gn/gnDebug.h"

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

boolean MemHashEntry::strict_start_lessthan_ptr(const MemHashEntry* a, const MemHashEntry* b){
	int32 start_diff = a->FirstStart() - b->FirstStart();
	if(start_diff == 0){
		uint32 m_count = a->m_start.size();
		m_count = m_count <= b->m_start.size() ? m_count : b->m_start.size();
		for(uint32 seqI = 0; seqI < m_count; seqI++){
			int64 a_start = a->m_start[seqI], b_start = b->m_start[seqI];
			if(a_start < 0)
				a_start = -a_start + a->m_length - a->m_mersize;
			if(b_start < 0)
				b_start = -b_start + b->m_length - b->m_mersize;
			int64 diff = a_start - b_start;
			if(diff == 0)
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

int64 MemHashEntry::end_to_start_compare(const MemHashEntry& a, const MemHashEntry& b){
	MemHashEntry tmp_a = a;
	tmp_a.CropStart(tmp_a.Length()-1);
	return MemHashEntry::start_compare(tmp_a, b);
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

MemHashEntry::MemHashEntry() : Match(){
	m_extended = false;
	m_mersize = 0;
	m_mim = NULL;
}

MemHashEntry::MemHashEntry(uint32 seq_count, const gnSeqI mersize, MemType m_type)
 : Match( seq_count )
{
	m_extended = m_type == extended;
	m_mersize = mersize;
}

MemHashEntry::MemHashEntry( const Match& mhe ){
	
	m_extended = extended;
	m_mersize = 0;
	m_length = mhe.Length();
	m_start.reserve( mhe.SeqCount() );
	for( uint startI = 0; startI < mhe.SeqCount(); startI++ ){
		m_start.push_back( mhe[ startI ] );
	}
	m_mim = NULL;
	CalculateOffset();
//	m_offsets = mhe.m_offsets;
//	m_offset = mhe.m_offset;
//	m_multiplicity = mhe.Multiplicity();
//	m_firstStart = mhe.m_firstStart;
//	m_matchnumber = mhe.m_matchnumber;
}

MemHashEntry::MemHashEntry(const MemHashEntry& mhe){
	(*this) = mhe;
}

MemHashEntry::~MemHashEntry(){

}
MemHashEntry* MemHashEntry::Clone() const{
	return new MemHashEntry(*this);
}
MemHashEntry& MemHashEntry::operator=(const MemHashEntry& mhe)
{
	m_offsets = mhe.m_offsets;
	m_offset = mhe.m_offset;
	m_extended = mhe.m_extended;
	m_mersize = mhe.m_mersize;
	m_length = mhe.m_length;
	m_start = mhe.m_start;
	m_multiplicity = mhe.m_multiplicity;
	m_firstStart = mhe.m_firstStart;
	m_subsets = mhe.m_subsets;
	m_supersets = mhe.m_supersets;
	m_mim = mhe.m_mim;
	m_matchnumber = mhe.m_matchnumber;
	return *this;
}

boolean MemHashEntry::operator==(const MemHashEntry& mhe) const{
	if( !Match::operator==( mhe ) )
		return false;
	if(m_mersize != mhe.m_mersize)
		return false;
	if(m_extended != mhe.m_extended)
		return false;
	return true;
}

void MemHashEntry::LinkSubset( MemHashEntry* subset ){
	m_subsets.insert( subset );
	subset->AddSuperset( this );
}

void MemHashEntry::AddSuperset( MemHashEntry* overlapper ){
	m_supersets.insert(overlapper);
}

/*
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

//check if this mem evenly overlaps mhe in every sequence for which
//this match is defined.
boolean MemHashEntry::EvenlyOverlaps(const MemHashEntry& mhe) const{
	uint32 i;
	int64 diff_i;
	int64 diff;
	vector<int64> m_start1 = m_start, m_start2 = mhe.m_start;
	uint32 seq_count = m_start2.size();

	//ensure the first sequence in this mem exists in both mems...
	i = FirstStart();
	if(m_start2[i] == MEM_NO_MATCH)
		return false;

	//make sure the parity of each genome is the same
	if(FirstStart() != mhe.FirstStart())
		HomogenizeParity(m_start1, m_start2, i);
	diff = m_start2[i] - m_start1[i];

	//check for even overlap properties
	if(diff >= m_length || -diff >= mhe.m_length)
		return false;

	//everything is ok so far, check for alignment
	int64 diff_rc = (int64)mhe.m_length - (int64)m_length + diff;
	for(i++; i < seq_count; i++){
		//check for consistent alignment between all genomes
		//in the case of revcomp, diff_i must equal m.length - diff
		diff_i = m_start2[i] - m_start1[i];
		//it's ok if mhe has a sequence which this mem doesn't
		//but not vice versa
		if(m_start1[i] == MEM_NO_MATCH)
			continue;
		if(m_start2[i] == MEM_NO_MATCH)
			return false;
		else if(m_start2[i] < 0 && diff_rc == diff_i)
			continue;
		else if(diff != diff_i )
			return false;
	}
	//it was an even overlap
	return true;
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
	//make sure they actually match in at least one genome.
	if(i==seq_count)
		return false;
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
*/
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

	uint fm = i;

	//make sure the parity of each genome is forward
//	if(FirstStart() != mhe.FirstStart())
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
	if( m_start[ fm ] > 0 ){
		unique_mhe.CropEnd( m_length - diff );
	}else{
		unique_mhe.CropEnd( mhe.m_length + diff );
	}
	return true;
}
boolean MemHashEntry::GetUniqueEnd(const MemHashEntry& mhe, MemHashEntry& unique_mhe) const{
	uint32 i;
	int64 diff_i;
	int64 diff;
	uint32 seq_count = mhe.m_start.size();
	vector<int64> m_start1 = m_start, m_start2 = mhe.m_start;
	uint fm;
	
	//find the first matching sequence...
	for(i=0; i < seq_count; i++){
		if(mhe.m_start[i] == MEM_NO_MATCH || m_start[i] == MEM_NO_MATCH)
			continue;
		else
			break;
	}
	fm = i;

	//make sure the parity of each genome is the same
//	if(FirstStart() != mhe.FirstStart())
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
	if( m_start[ fm ] > 0 ){
		unique_mhe.CropStart(mhe.m_length + diff);
	}else{
		unique_mhe.CropStart( m_length - diff );
	}
	return true;
}

void MemHashEntry::UnlinkSubset( MemHashEntry* subset ) {
	m_subsets.erase( subset );
}

void MemHashEntry::UnlinkSuperset( MemHashEntry* superset ) {
	m_supersets.erase( superset );
}

void MemHashEntry::UnlinkSelf() {
	// unlink from the supersets
	set<MemHashEntry*>::iterator match_iter = m_supersets.begin();
	for(; match_iter != m_supersets.end(); match_iter++ )
		(*match_iter)->UnlinkSubset( this );
	m_supersets.clear();

	// unlink from the subsets
	match_iter = m_subsets.begin();
	for(; match_iter != m_subsets.end(); match_iter++ )
		(*match_iter)->UnlinkSuperset( this );
	m_subsets.clear();

}

void MemHashEntry::UpdateSupersets() {
	set<MemHashEntry*>::iterator super_iter = m_supersets.begin();

	while(super_iter != m_supersets.end()) {			
		MemHashEntry* supermem = *super_iter;
		if( EvenlyOverlaps(*supermem) ){
			supermem->m_subsets.insert( this );
			super_iter++;
		}else {
			set<MemHashEntry*>::iterator to_del = super_iter;
			super_iter++;
			m_supersets.erase( to_del );
		}
	}
}

void MemHashEntry::EliminateSubsetInclusions() {
}


int64 MemHashEntry::GapSize(const MemHashEntry& mhe, uint32 seqI){

	int64 start_1 = m_start[seqI];
	int64 start_2 = mhe.m_start[seqI];

	//check for match existence
	//and for similar match direction
	if((start_2 < 0 && start_1 > 0) ||
		(start_2 > 0 && start_1 < 0) ||
		(start_2 == 0 && start_1 != 0) ||
		(start_2 != 0 && start_1 == 0))
		Throw_gnEx(SeqIndexOutOfBounds());
	
	// this code can probably be simplified
	if(start_2 > 0)
		return start_2 - (start_1 + (int64)m_length);
	if(start_2 < 0)
		return -start_1 - (-start_2 + mhe.m_length);
	else
		return 0;
}

/*
// checks if mhe is _perfectly_ contained in this match.
// all offsets in all sequences must be aligned to each other
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
	diff = mhe.m_start[i] - m_start[i];
	if(m_start[i] == MEM_NO_MATCH)
		return false;

	//check for containment properties
	if(diff < 0 || m_length < mhe.m_length + diff)
		return false;

	//everything is ok so far, check for alignment
	int64 diff_rc = (int64)mhe.m_length - (int64)m_length + diff;
	for(i++; i < seq_count; i++){
		//check for consistent alignment between all genomes
		//in the case of revcomp, diff_i must equal diff_rc
		diff_i = mhe.m_start[i] - m_start[i];

		//it's ok if neither matches in a sequence
		if(mhe.m_start[i] == MEM_NO_MATCH && m_start[i] == MEM_NO_MATCH)
			continue;
		else if(mhe.m_start[i] < 0 && diff_rc == diff_i)
			continue;
		else if(diff != diff_i )
			return false;
	}
	//it was contained.
	return true;
}
*/
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
	int64 diff_rc = (int64)mhe.m_length - (int64)m_length + diff;
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


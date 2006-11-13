#include "Match.h"
#include "gn/gnException.h"
#include "gn/gnDebug.h"

MatchID_t Match::GetNewMatchID(){
	static MatchID_t next_id = 0;
	next_id++;
	return next_id;
}


Match::Match(){
	m_offset = 0;
	m_length = 0;
	m_multiplicity = 0;
	m_firstStart = UINT32_MAX;
	match_id = GetNewMatchID();
}

Match::Match( uint32 seq_count ){
	m_offset = 0;
	m_length = 0;
	m_multiplicity = 0;
	m_firstStart = UINT32_MAX;
	m_start.reserve( seq_count );
	for( uint32 seqI = 0; seqI < seq_count; seqI++ )
		m_start.push_back( MEM_NO_MATCH );

	match_id = GetNewMatchID();
}


Match::Match( uint32 seq_count, MatchID_t match_id ){
	m_offset = 0;
	m_length = 0;
	m_multiplicity = 0;
	m_firstStart = UINT32_MAX;
	m_start.reserve( seq_count );
	for( uint32 seqI = 0; seqI < seq_count; seqI++ )
		m_start.push_back( MEM_NO_MATCH );

	this->match_id = match_id;
}

Match::Match(const Match& mhe){
	(*this) = mhe;
}

Match::~Match(){

}

Match* Match::Clone() const{
	return new Match(*this);
}

Match& Match::operator=(const Match& mhe)
{
	m_offsets = mhe.m_offsets;
	m_offset = mhe.m_offset;
	m_length = mhe.m_length;
	m_start = mhe.m_start;
	m_multiplicity = mhe.m_multiplicity;
	m_firstStart = mhe.m_firstStart;
	subsets = mhe.subsets;
	supersets = mhe.supersets;
	m_matchnumber = mhe.m_matchnumber;
	match_id = mhe.match_id;
	return *this;
}

boolean Match::operator==(const Match& mhe) const{
	if(m_length != mhe.m_length)
		return false;
	if(m_start.size() != mhe.m_start.size())
		return false;
	for(uint32 startI = 0; startI < m_start.size(); startI++){
		if(m_start[startI] != mhe.m_start[startI])
			return false;
	}
	return true;
}

int64 Match::End(uint32 endI) const{
	if( m_start[ endI ] > 0 )
		return m_start[endI] + m_length - 1;
	return m_start[ endI ];
}

void Match::Invert(){
	for(uint32 startI = 0; startI < m_start.size(); startI++)
		m_start[startI] = -m_start[startI];
}

void Match::LinkSubset( Match* subset ){
	subsets.insert( subset->MatchID() );
	subset->AddSuperset( match_id );
}

void Match::AddSubset( MatchID_t underlapper ){
	subsets.insert( underlapper );
}

void Match::AddSuperset( MatchID_t overlapper ){
	supersets.insert( overlapper );
}

uint32 Match::FirstStart() const{
	if(m_firstStart != UINT32_MAX)
		return m_firstStart;
	for(uint32 startI=0; startI < m_start.size(); startI++)
		if(m_start[startI] != MEM_NO_MATCH)
			return startI;
	return UINT32_MAX;
}

void Match::HomogenizeParity(vector<int64>& start1, vector<int64>& start2, const uint32 startI) const{
	//make sure the parity of each genome is the forward when the start level is different 
	if( start1[startI] < 0 ){
		//invert m_start1...
		for(uint32 j = startI; j < start1.size(); j++)
			start1[j] *= -1;
	}
	if( start2[startI] < 0 ){
		//invert m_start2...
		for( uint32 j=startI; j < start2.size(); j++)
			start2[j] *= -1;
	}
}

//check if this mem evenly overlaps mhe in every sequence for which
//this match is defined.
boolean Match::EvenlyOverlaps(const Match& mhe) const{
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

boolean Match::Intersects(const Match& mhe) const{
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

boolean Match::GetUniqueStart(const Match& mhe, Match& unique_mhe) const{
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

boolean Match::GetUniqueEnd(const Match& mhe, Match& unique_mhe) const{
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

void Match::UnlinkSubset( MatchID_t subset ) {
	subsets.erase( subset );
}

void Match::UnlinkSuperset( MatchID_t superset ) {
	supersets.erase( superset );
}

int64 Match::GapSize(const Match& mhe, uint32 seqI){

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

void Match::Move( int64 distance ){
	for( uint32 i=0; i < m_start.size(); i++ ){
		if( m_start[i] != MEM_NO_MATCH )
			m_start[i] += distance;
	}
}

void Match::CropStart(gnSeqI crop_amount){
	if( crop_amount > m_length )
		Throw_gnEx( SeqIndexOutOfBounds() );
	for( uint32 i=0; i < m_start.size(); i++ ){
		if( m_start[i]  > 0 )
			m_start[i] += crop_amount;
	}
	m_length -= crop_amount;
}

void Match::CropEnd(gnSeqI crop_amount){
	if( crop_amount > m_length )
		Throw_gnEx( SeqIndexOutOfBounds() );
	for( uint32 i=0; i < m_start.size(); i++ ){
		if( m_start[i]  < 0 )
			m_start[i] -= crop_amount;
	}
	m_length -= crop_amount;
}

// checks if mhe is _perfectly_ contained in this match.
// all offsets in all sequences must be aligned to each other
boolean Match::Contains(const Match& mhe) const{
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

void Match::CalculateOffset(){
	int64 ref_start;
	uint32 seq_count = m_start.size();
	uint32 seqI = 0;
	m_offset = 0;
//	m_matchnumber = valarray<bool>( false, seq_count );
	m_offsets.clear();
	m_offsets.reserve( seq_count );
	for(; seqI < seq_count; seqI++){
		m_offsets.push_back( 0 );
		if(m_start[seqI] != MEM_NO_MATCH){
			ref_start = m_start[seqI];
			m_multiplicity = 1;
			m_firstStart = seqI;
//			m_matchnumber[ seqI ] = true;
			break;
		}
	}
	m_matchnumber = 1;
	for(seqI++; seqI < seq_count; seqI++){
		m_matchnumber <<= 1;
		if(m_start[seqI] > 0){
			m_offsets[ seqI ] = m_start[seqI] - ref_start;
			m_offset += m_offsets[ seqI ];
			m_multiplicity++;
			m_matchnumber |= 1;
		}
		else if( m_start[seqI] < 0 ){
			m_offsets[ seqI ] = m_start[seqI] - (int64)m_length - ref_start;
			m_offset += m_offsets[ seqI ];
			m_multiplicity++;
			m_matchnumber |= 1;
		}
	}
}

void Match::ExtendStart(gnSeqI extend_amount){
	m_length += extend_amount;
	uint32 seq_count = m_start.size();
	for(uint32 seqI = 0; seqI < seq_count; seqI++)
		if(m_start[seqI] > 0)
			m_start[seqI] -= extend_amount;
}

void Match::ExtendEnd(gnSeqI extend_amount){
	m_length += extend_amount;
	uint32 seq_count = m_start.size();
	for(uint32 seqI = 0; seqI < seq_count; seqI++)
		if(m_start[seqI] < 0)
			m_start[seqI] += extend_amount;
}

ostream& operator<<(ostream& os, const Match& mhe){ //write to stream.
	os << mhe.m_length;
	for(uint32 i=0; i < mhe.m_start.size(); i++)
		os << '\t' << mhe.m_start[i];
	return os;
}

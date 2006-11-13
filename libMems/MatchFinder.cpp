#include "MatchFinder.h"

MatchFinder::MatchFinder(){
	mer_size = 0;
	seq_count = 0;
	ambiguity_tolerance = 0;
	m_progress = -1;
	log_stream = NULL;
}

//make sure this calls the destructor on each element
MatchFinder::~MatchFinder(){
}

MatchFinder::MatchFinder(const MatchFinder& mf){
	mer_size = mf.mer_size;
	seq_count = mf.seq_count;
	ambiguity_tolerance = mf.ambiguity_tolerance;

	m_progress = mf.m_progress;
	sar_table = mf.sar_table;
	seq_table = mf.seq_table;
	log_stream = mf.log_stream;
}

void MatchFinder::Clear(){
	mer_size = 0;
	seq_count = 0;
	ambiguity_tolerance = 0;
	m_progress = -1;
	sar_table.clear();
	seq_table.clear();
	log_stream = NULL;
}

void MatchFinder::LogProgress( ostream* os ){
	log_stream = os;
}

boolean MatchFinder::AddSequence( SortedMerList* sar, gnSequence* seq ){
	if(sar == NULL){
		Throw_gnExMsg( NullPointer(), "Null SortedMerList pointer" );
	}
	if(sar == NULL){
		Throw_gnExMsg( NullPointer(), "Null gnSequence pointer" );
	}
	
	//check for consistency between sequence length and sorted mer list lengths
	if(seq != NULL && seq->length() != sar->Length()){
		DebugMsg("MatchFinder::AddSequence: Error mismatched sml and sequence length.");
		return false;
	}
	
	//passed checks, add it to the data structures
	sar_table.push_back(sar);
	seq_count++;
	if(seq != NULL){
		seq_table.push_back(seq);
	}
	
	SMLHeader header = sar->GetHeader();
	alphabet_bits = header.alphabet_bits;
	
	return true;

}

void MatchFinder::GetBreakpoint( uint32 sarI, gnSeqI startI, vector<gnSeqI>& breakpoints ) const{
	breakpoints.clear();
	
	//put the mer to break on in break_mer
	bmer break_mer  = (*GetSar(sarI))[startI];
	uint64 mer_mask = GetSar(sarI)->GetMerMask();
	bmer prev_mer = break_mer;
	//search backwards for the first index of this mer
	while((prev_mer.mer & mer_mask) == (break_mer.mer & mer_mask)){
		if(startI == 0){
			startI--;
			break;
		}
		startI--;
		prev_mer = (*GetSar(sarI))[startI];
	}
	startI++;

	//find the mer's location in the other sorted mer lists
	for(uint32 i=0; i < seq_count; i++){
		if(i == sarI){
			breakpoints.push_back(startI);
		}else{
			gnSeqI cur_start;
			if(GetSar(i)->FindMer(break_mer.mer, cur_start)){
				//we found a match, see how far backwards we can go.
				int64 cur_matchI = cur_start;
				bmer matchmer = (*GetSar(i))[cur_start];
				while(cur_matchI >= 0 && ((matchmer.mer & mer_mask) == (break_mer.mer && mer_mask))){
					cur_matchI--;
					matchmer = (*GetSar(i))[cur_start];
				}
				cur_start = cur_matchI+1;
			}
			breakpoints.push_back(cur_start);
		}
	}
}

boolean MatchFinder::FindMatches(){
	vector<gnSeqI> start_points;
	vector<gnSeqI> search_len;
	for(uint32 i=0; i < sar_table.size(); i++){
		start_points.push_back(0);
		search_len.push_back(GNSEQI_END);
	}
	while( !FindMatches(start_points, search_len) ){;}
	return true;
}

#define MER_REPEAT_LIMIT 10000 // The maximum number of matching mers before they are completely
								// ignored.

/**
 * Template class for an array of vectors which ensures correct deallocation.
 */
template<class T>
class VectorArray{
public:
	/** Allocate a new array of length "bufsize".  The array can be
	 * accessed through the data member
	 */
	VectorArray( uint64 vector_count ){
		v = new vector<T>[vector_count];
	}
	~VectorArray() { 
		if (v != NULL)
			delete[] v;
	}
	/** The actual array of vectors */
	vector<T>* v;
private:
	VectorArray( const VectorArray& a );
	VectorArray& operator=( const VectorArray& a );
	VectorArray(){};
};

//startI must be 0
//At most search_length mers in any one genome will be checked.
boolean MatchFinder::FindMatches(vector<gnSeqI>& start_points, vector<gnSeqI>& search_len){
	//picked a semi-arbitrary number for buffer size.
	uint32 MER_BUFFER_SIZE = 10000;
	vector<uint32> mer_index;   // stores the indexes of the current mers in mer_vector
	vector<uint32> mer_baseindex;   // stores the index in the SortedMerList of each of the first mers in mer_vector
	list<idmer> cur_mers;	// stores the current mers.
	list<idmer> cur_match;	// stores the current matching mers.
	list<uint32> sar_hitlist;	// list of sars to replace
	uint32 read_size;
	
	//make sure there is at least one sequence
	if(sar_table.size() < 1)
		return true;
	
	//check for consistency in mer sizes.
	uint64 mer_mask = sar_table[0]->GetMerMask();
	mer_size = sar_table[0]->MerSize();
	for(uint32 maskI = 0; maskI < sar_table.size(); maskI++){
		if(mer_size != sar_table[maskI]->MerSize()){
			Throw_gnExMsg(InvalidData(), "Unequal mer sizes.");
		}
	}
	
	//check that start_points and end_points are ok.
	if((start_points.size() != sar_table.size()) || (search_len.size() != sar_table.size())){
		Throw_gnExMsg(InvalidData(), "Inconsistent search range specification.");
	}
	
	// initialize subset match data structures
	alpha_map_size = (uint)pow( (double)2, (double)alphabet_bits );
	alpha_map.reserve( alpha_map_size );
	vector< uint32 > tmp_list;
	tmp_list.reserve( seq_count );
	for( uint alphaI = 0; alphaI < alpha_map_size; alphaI++ ){
		alpha_map.push_back( tmp_list );
	}
	
	//allocate buffer space
	// stores arrays of bmers for each sml.
	VectorArray<bmer> mer_vector( sar_table.size() );

	// keep track of the number of mers processed and the total for progress reporting
	uint64 mers_processed = 0;
	uint64 total_mers = 0;

	//initialize the data structures
	idmer newmer;
	for(uint32 n = 0; n < sar_table.size(); n++){
		read_size = MER_BUFFER_SIZE < search_len[n] ? MER_BUFFER_SIZE : search_len[n]; 
		mer_vector.v[n].reserve(read_size);
		sar_table[n]->Read(mer_vector.v[n], read_size, start_points[n]);
		mer_index.push_back(0);
		mer_baseindex.push_back(0);
		newmer.position = mer_vector.v[n][0].position;
		newmer.mer = mer_vector.v[n][0].mer & mer_mask;
		newmer.id = n;
		cur_mers.push_back(newmer);  //cur_mers gets the first mer from each sorted mer list
		total_mers += search_len[n] == GNSEQI_END ? seq_table[n]->length() : search_len[n];
	}
	
	//nobody reads these fucking things.  why am i writing this.because my fucking 
	//roomate needs a goddamn roadmap......   ohhh ecstasy.... haptic pumpkins

	//loop while there is data to hash.
	cur_mers.sort(&idmer_lessthan);
	while(cur_mers.size() > 0){
		list<idmer>::iterator mer_iter = cur_mers.begin();
		sarID_t cur_id = mer_iter->id;
		//first check for matches across genomes.
		if(cur_match.size() > 0){
			if(mer_iter->mer > cur_match.begin()->mer){
				//we are done with this matching.  hash it.
				if(cur_match.size() > 1)
					EnumerateMatches(cur_match);
				cur_match.clear();
			}else if(mer_iter->mer < cur_match.begin()->mer){
				//error checking stuff.
				ErrorMsg("Horrible error occurred!!\n");
			}
		}

		if( cur_match.size() > MER_REPEAT_LIMIT ){
			// scan past the repetitive mers
			// create the lexicographically next mer
			uint64 next_mer = cur_match.begin()->mer;
			next_mer += ~mer_mask + 1;
			gnSeqI next_pos = 0;
			if( !sar_table[0]->FindMer( next_mer, next_pos ))
				next_pos++;
			GetBreakpoint( 0, next_pos, start_points );
			return false;
		}
		//check for matches within the same genome
		gnSeqI merI = mer_index[cur_id];
		boolean buffer_exhausted = merI < mer_vector.v[cur_id].size() ? false : true;
		while(!buffer_exhausted && (mer_iter->mer == (mer_vector.v[cur_id][merI].mer & mer_mask))){
			newmer.position = mer_vector.v[cur_id][merI].position;
			newmer.mer = mer_vector.v[cur_id][merI].mer & mer_mask;
			newmer.id = cur_id;
			cur_match.push_back(newmer);
			merI++;
			mer_index[cur_id]++;
			//check if we've exhausted our buffer
			if(merI == mer_vector.v[cur_id].size())
				buffer_exhausted = true;
		}

		if(buffer_exhausted){
			//if we've exhausted our buffer then refill it
			mer_baseindex[cur_id] += mer_vector.v[cur_id].size();
			
			// update the mers processed
			mers_processed += mer_vector.v[cur_id].size();
			float64 m_oldprogress = m_progress;
			m_progress = ((float64)mers_processed / (float64)total_mers) * PROGRESS_GRANULARITY;
			if( log_stream != NULL ){
				if((int)m_oldprogress != (int)m_progress)
					(*log_stream) << (int)((m_progress / PROGRESS_GRANULARITY) * 100) << "%..";
			}
			uint32 read_size = MER_BUFFER_SIZE;
			if(MER_BUFFER_SIZE + mer_baseindex[cur_id] > search_len[cur_id])
				read_size = search_len[cur_id] - mer_baseindex[cur_id];

			sar_table[cur_id]->Read(mer_vector.v[cur_id], read_size, start_points[cur_id] + mer_baseindex[cur_id]);
			mer_index[cur_id] = 0;
			if(mer_vector.v[cur_id].size() == 0){
				//remove mer_iter so that this sar is forgotten
				cur_mers.erase(mer_iter);
			}
		}else{
			//if we haven't exhausted our buffer then we must have
			//run out of matching mers.
			//remove mer_iter and put in a new idmer with the same id
			cur_mers.erase(mer_iter);
			newmer.position = mer_vector.v[cur_id][merI].position;
			newmer.mer = mer_vector.v[cur_id][merI].mer & mer_mask;
			newmer.id = cur_id;
			mer_iter = cur_mers.begin();
			while(mer_iter->mer < newmer.mer && mer_iter != cur_mers.end())
				mer_iter++;
			cur_mers.insert(mer_iter, newmer);
		}
		
	}
	return true;
}

boolean MatchFinder::EnumerateMatches( list<idmer>& match_list ){
	//this must call HashMatch on every possible combination of matches in the list.
	if(match_list.size() == 2){
		//this is the smallest possible match.  simply hash it.
		return HashMatch(match_list);
	}
	
	match_list.sort(&idmer_id_lessthan);
	vector<uint32> id_start;
	vector<list<idmer>::iterator> id_pos;
	vector<list<idmer>::iterator> id_end;
	list<idmer>::iterator iter = match_list.begin();
	list<idmer>::iterator iter2 = match_list.begin();
	iter2++;
	id_start.push_back(0);
	id_pos.push_back(iter);
	for(uint32 i=0; iter2 != match_list.end(); i++){
		if(iter->id != iter2->id){
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

// takes as input a fully extended mem and returns the subset matches on the lower side
void MatchFinder::FindSubsets(const MemHashEntry& mhe, vector<MemHashEntry>& subset_matches){

	SMLHeader head = GetSar( 0 )->GetHeader();
	uint shift_amt = 64 - head.alphabet_bits;
	uint rshift_amt = head.alphabet_bits * ( mer_size - 1 );

	uint seqI, alphaI;
	
	
	for( seqI = 0; seqI < seq_count; seqI++ ){
		//check that all mers at the new position match
		int64 mer_to_get = mhe[ seqI ];
		if( mer_to_get == MEM_NO_MATCH )
			continue;
		if(mer_to_get < 0){
			mer_to_get *= -1;
			mer_to_get += mhe.Length() - mer_size;
		}

		uint64 cur_mer = GetSar( seqI )->GetMer( mer_to_get - 1 );

		boolean parity;
		if( mhe[ seqI ] < 0 )
			parity = cur_mer & 0x1;
		else
			parity = !(cur_mer & 0x1);

		if( parity ){
			cur_mer >>= shift_amt;
		}else{
			cur_mer <<= rshift_amt;
			cur_mer = ~cur_mer;
			cur_mer >>= shift_amt;
		}

		alpha_map[ cur_mer ].push_back( seqI );

	}
	
	for( alphaI = 0; alphaI < alpha_map_size; alphaI++ ){
		if( alpha_map[ alphaI ].size() < 2 ){
			alpha_map[ alphaI ].clear();
			continue;
		}
		// this is a subset
		MemHashEntry cur_subset = MemHashEntry( mhe.SeqCount(), mhe.MerSize() );
		cur_subset.SetLength( mhe.Length() );

		for( uint subI = 0; subI < alpha_map[ alphaI ].size(); subI++ )
			cur_subset.SetStart( alpha_map[ alphaI ][ subI ], mhe[ alpha_map[ alphaI ][ subI ] ] );
		subset_matches.push_back( cur_subset );
		alpha_map[ alphaI ].clear();
	}
}

/*
// takes as input a fully extended mem and returns the subset matches on the lower side
void MatchFinder::FindSubsets(const MemHashEntry& mhe, vector<MemHashEntry>& subset_matches){
	// stores each unique mer which breaks the match
	uint32 next_subset = mhe.FirstStart();
	boolean found_subset;
	boolean found_match;
	uint64 cur_mer;
	uint64 mer_mask = GetSar(0)->GetMerMask();
	if( mhe[ 0 ] == 1086895 && mhe.Length() == 654 )
		__asm( nop );
	// stores the coordinates of the current subset
	MemHashEntry cur_subset, tmp_mhe(mhe);

	while( next_subset != mhe.SeqCount() ){
		cur_subset = MemHashEntry(tmp_mhe.SeqCount(), tmp_mhe.MerSize());
		cur_subset.SetLength(tmp_mhe.Length());
		cur_subset.SetStart( next_subset, tmp_mhe[next_subset] );

		if( tmp_mhe[ next_subset ] == MEM_NO_MATCH ){
			DebugMsg("Lethal error!\n");
		}

		//check that all mers at the new position match
		int64 mer_to_get = tmp_mhe[ next_subset ];
		if(mer_to_get < 0){
			mer_to_get *= -1;
			mer_to_get += tmp_mhe.Length() - mer_size;
		}
		cur_mer = GetSar( next_subset )->GetMer( mer_to_get - 1 ) & mer_mask;

		found_match = false;
		found_subset = false;
		next_subset++;
		for(uint32 subI = next_subset; subI < tmp_mhe.SeqCount(); subI++){
			if( !found_subset )
				next_subset++;

			mer_to_get = tmp_mhe[ subI ];
			if( mer_to_get == MEM_NO_MATCH )
				continue;
			else if(mer_to_get < 0){
				//Convert the cur_seqs[i] entry since negative implies reverse complement
				mer_to_get *= -1;
				mer_to_get += tmp_mhe.Length() - mer_size;
			}

			if(cur_mer != (GetSar( subI )->GetMer(mer_to_get - 1) & mer_mask)){
				found_subset = true;
			}else{
				cur_subset.SetStart( subI, tmp_mhe[subI] );
				tmp_mhe.SetStart( subI, MEM_NO_MATCH );
				found_match = true;
			}
		}
		if( found_match )
			subset_matches.push_back( cur_subset );
		if( found_subset )
			next_subset--;
	}
}
*/
//mems which span the origin will be hashed a second time
void MatchFinder::ExtendMatch(MemHashEntry& mhe, vector<MemHashEntry>& subset_matches, gnSeqI max_backward, gnSeqI max_forward){
	uint64 cur_mer;
	uint64 mer_mask = GetSar(0)->GetMerMask();

	//which sequences are used in this match?
	uint32* cur_seqs = new uint32[mhe.SeqCount()];
	uint32 used_seqs = 0;
	for(uint32 seqI = 0; seqI < mhe.SeqCount(); seqI++){
		if(mhe[seqI] != MEM_NO_MATCH){
			cur_seqs[used_seqs] = seqI;
			used_seqs++;
		}
	}
	//First extend backwards then extend forwards.  The following loop does them both.
	gnSeqI jump_size = mer_size;
	for(uint32 directionI = 0; directionI < 4; directionI++){
		//how far can we go?	
		//first calculate the maximum amount of traversal
		//then do fewer comparisons.
		int64 maxlen = GNSEQI_END;
		if(directionI == 0)
			maxlen = max_backward;
		else if(directionI == 1)
			maxlen = max_forward;
		else
			maxlen = mer_size;

		for(uint32 maxI = 0; maxI < used_seqs; maxI++)
			if(GetSar(cur_seqs[maxI])->IsCircular()){
				if(GetSar(cur_seqs[maxI])->Length() < maxlen)
					maxlen = GetSar(cur_seqs[maxI])->Length();
			}else if(mhe[cur_seqs[maxI]] < 0){
				int64 rc_len = GetSar(cur_seqs[maxI])->Length() - mhe.Length() + mhe[cur_seqs[maxI]] + 1;
				if( rc_len < maxlen)
					maxlen = rc_len;
			}else if(mhe[cur_seqs[maxI]] - 1 < maxlen)
				maxlen = mhe[cur_seqs[maxI]] - 1;
		uint32 j=0;
		uint32 i = used_seqs;	// set to used_seqs in case maxlen is already less than jump size.
		// FIX ME:
		// when maxlen runs out then we must drop the offending sequences and process the subset.
		while(maxlen - jump_size >= 0){
			mhe.SetLength(mhe.Length() + jump_size);
			maxlen -= jump_size;
			for(j=0; j < used_seqs; j++)
				if(mhe[cur_seqs[j]] > 0){
					mhe.SetStart(cur_seqs[j], mhe[cur_seqs[j]] - jump_size);
					if(mhe[cur_seqs[j]] <= 0)
						mhe.SetStart(cur_seqs[j], mhe[cur_seqs[j]] + GetSar(cur_seqs[j])->Length());
				}
			//check that all mers at the new position match
			int64 mer_to_get = mhe[cur_seqs[0]];
			if(mer_to_get < 0){
				mer_to_get *= -1;
				mer_to_get += mhe.Length() - mer_size;
			}
			cur_mer = GetSar(cur_seqs[0])->GetMer(mer_to_get - 1);
			boolean parity;
			if( mhe[cur_seqs[0]] < 0 )
				parity = cur_mer & 0x1;
			else
				parity = !(cur_mer & 0x1);
			cur_mer &= mer_mask;

			for(i=1; i < used_seqs; i++){
				mer_to_get = mhe[cur_seqs[i]];
				if(mer_to_get < 0){
					//Convert the cur_seqs[i] entry since negative implies reverse complement
					mer_to_get *= -1;
					mer_to_get += mhe.Length() - mer_size;
				}
				uint64 comp_mer = GetSar(cur_seqs[i])->GetMer(mer_to_get - 1);
				boolean comp_parity;				
				if( mhe[cur_seqs[i]] < 0 )
					comp_parity = comp_mer & 0x1;
				else
					comp_parity = !(comp_mer & 0x1);
				comp_mer &= mer_mask;
				
				if(cur_mer != comp_mer || parity != comp_parity ){
					if( directionI > 1 && used_seqs > 2 )
						FindSubsets( mhe, subset_matches );
					maxlen = 0;
					break;
				}
			}
		}
		//this stuff cleans up if there was a mismatch
		if(i < used_seqs){
			mhe.SetLength(mhe.Length() - jump_size);
			for(;j > 0; j--)
				if(mhe[cur_seqs[j - 1]] >= 0)
					mhe.SetStart(cur_seqs[j - 1], mhe[cur_seqs[j - 1]] + jump_size);
		}
		//Invert the sequence directions so that we extend in the other direction
		//next time through the loop.  The second time we do this we are setting
		//sequence directions back to normal.
		mhe.Invert();

		//if we've already been through twice then decrease the jump size
		if(directionI >= 1)
			jump_size = 1;
	}

	// set the subsets so their reference sequence is always positive
	for(uint32 subsetI = 0; subsetI < subset_matches.size(); subsetI++){
		if( subset_matches[subsetI][subset_matches[subsetI].FirstStart()] < 0 )
			subset_matches[subsetI].Invert();
		subset_matches[subsetI].CalculateOffset();
	}

	delete[] cur_seqs;
}

boolean MatchFinder::MatchAmbiguities(MemHashEntry& mhe, uint32 match_size){
	if(ambiguity_tolerance == 0)
		return false;
			//check that all mers at the new position match
	//which sequences are used in this match?
	uint32* cur_seqs = new uint32[mhe.SeqCount()];
	uint32 used_seqs = 0;
	for(uint32 seqI = 0; seqI < mhe.SeqCount(); seqI++){
		if(mhe[seqI] != MEM_NO_MATCH){
			cur_seqs[used_seqs] = seqI;
			used_seqs++;
		}
	}
	string cur_mer, mer_i;
	gnSequence mer_seq;
	int64 mer_to_get = mhe[cur_seqs[0]];
	if(mer_to_get < 0){
		mer_to_get *= -1;
		mer_to_get += mhe.Length() - mer_size;
	}
	cur_mer = seq_table[cur_seqs[0]]->subseq(mer_to_get, match_size).ToString();
	
	for(uint32 i=1; i < used_seqs; i++){
		mer_to_get = mhe[cur_seqs[i]];
		if(mer_to_get < 0){
			//Convert the cur_seqs[i] entry since negative implies reverse complement
			mer_to_get *= -1;
			mer_to_get += mhe.Length() - mer_size;
		}
		mer_seq = seq_table[cur_seqs[i]]->subseq(mer_to_get, match_size);
		if(mer_seq.compare(cur_mer) != 0){
			delete[] cur_seqs;
			return false;
		}
		mer_i = mer_seq.ToString();
		uint32 ambiguity_count = 0;
		for(uint32 baseI = 0; baseI < match_size; baseI++)
			if(cur_mer[baseI] != mer_i[baseI])
				ambiguity_count++;
		if(ambiguity_count > ambiguity_tolerance){
			delete[] cur_seqs;
			return false;
		}
	}
	delete[] cur_seqs;
	return true;
}


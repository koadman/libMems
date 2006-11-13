#include "HashingMatchFinder.h"

HashingMatchFinder::HashingMatchFinder(){
	mer_size = 0;
	seq_count = 0;
	mer_table_size = DEFAULT_MER_TABLE_SIZE;
	mer_table = new list<idmer>[mer_table_size];
}

HashingMatchFinder::~HashingMatchFinder(){
	//make sure this calls the destructor on each element
	delete[] mer_table;
}

HashingMatchFinder* HashingMatchFinder::Clone() const{

}

boolean HashingMatchFinder::AddSuffixArray( SuffixArray* sar, gnSequence* seq ){
	if(sar == NULL){
		DebugMsg("HashingMatchFinder::AddSuffixArray: Error null parameter.");
		return false;
	}
	SarHeader header = sar->GetHeader();
	for(uint32 i = 0; i < seqid_table.size(); i++){
		if(header.id == seqid_table[i]){
			DebugMsg("HashingMatchFinder::AddSuffixArray: Error seqid already exists.");
			return false;
		}
	}

	//check for consistency between sequence length and suffix array lengths
	if(seq != NULL && seq->length() != sar->Length()){
		DebugMsg("HashingMatchFinder::AddSuffixArray: Error mismatched sar and sequence length.");
		return false;
	}
	
	//passed checks, add it to the data structures
	sarid_table[(uint8)header.id] = sar_table.size();
	sar_table.push_back(sar);
	seq_count++;
	if(seq != NULL){
		seq_table.push_back(seq);
		seqid_table.push_back(header.id);
	}
	
	//put the mers in the hashtable.
	uint32 hash_key = ComputeHashKey(hash_mer);
	mer_table[hash_key].push_back(hash_mer);
	return true;
}

uint32 HashingMatchFinder::ComputeHashKey(idmer& hash_mer){
	//simple modulo will suffice for now.
	return hash_mer.mer	% mer_table_size;
}

boolean HashingMatchFinder::FindMatches(uint32 mer_mask_size, gnSeqI startI, gnSeqI endI){
	//use startI and endI as starting and ending indices of buckets in the hash table which will
	//be used.  This will allow efficient parallelization of the algorithm.
	
	uint32 end_bucket = endI < mer_table_size ? endI : mer_table_size;
	for(uint32 bucketI = startI; bucketI < end_bucket; bucketI++){
		//copy to array? then sort.
		mer_table[bucketI].sort(&idmer_lessthan);
		list<idmer> match_list;
		list<idmer>::iterator cur_mer = mer_table[bucketI].start();
		list<idmer>::iterator start_mer = cur_mer;
		list<idmer>::iterator end_mer = mer_table[bucketI].end();
		uint32 cur_match_count = 0;
		for(cur_mer++; cur_mer != end_mer; cur_mer++){
			if(cur_mer->mer != start_mer->mer){
				if(match_list.size() > 0){
					match_list.push_back(start_mer);
					EnumerateMatches(match_list);
					match_list.clear();
				}
				start_mer = cur_mer;
			}else
				match_list.push_back(*cur_mer);
		}
	}

	return true;
}


boolean HashingMatchFinder::EnumerateMatches( list<idmer>& match_list ){
	//this must call HashMatch on every possible combination of matches in the list.
	if(match_list.size() == 2){
		//this is the smallest possible match.  simply hash it.
		return HashMatch(match_list);
	}
	
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
//mems which span the origin will be hashed a second time

void HashingMatchFinder::ExtendMatch(MemHashEntry& mhe, gnSeqI max_backward, gnSeqI max_forward){
	uint64 cur_mer;
	uint64 mer_mask = sar_table[0]->GetMerMask();

	//which sequences are used in this match?
	uint32* cur_seqs = new uint32[seq_count];
	uint32 used_seqs = 0;
	for(uint32 seqI = 0; seqI < seq_count; seqI++){
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
			if(sar_table[cur_seqs[maxI]]->IsCircular()){
				if(sar_table[cur_seqs[maxI]]->Length() < maxlen)
					maxlen = sar_table[cur_seqs[maxI]]->Length();
			}else if(mhe[cur_seqs[maxI]] < 0){
				int64 rc_len = sar_table[cur_seqs[maxI]]->Length() - mhe.Length() -  mhe[cur_seqs[maxI]] - 1;
				if( rc_len < maxlen)
					maxlen = rc_len;
			}else if(mhe[cur_seqs[maxI]] - 1 < maxlen)
				maxlen = mhe[cur_seqs[maxI]] - 1;
		
		int32 j=0;
		uint32 i;
		while(maxlen - jump_size >= 0){
			mhe.SetLength(mhe.Length() + jump_size);
			maxlen -= jump_size;
			for(j=0; j < used_seqs; j++)
				if(mhe[cur_seqs[j]] > 0){
					mhe.SetStart(cur_seqs[j], mhe[cur_seqs[j]] - jump_size);
					if(mhe[cur_seqs[j]] <= 0)
						mhe.SetStart(cur_seqs[j], mhe[cur_seqs[j]] + sar_table[cur_seqs[j]]->Length());
				}
			//check that all mers at the new position match
			int64 mer_to_get = mhe[cur_seqs[0]];
			if(mer_to_get < 0){
				mer_to_get *= -1;
				mer_to_get += mhe.Length() - mer_size;
			}
			cur_mer = sar_table[cur_seqs[0]]->GetMer(mer_to_get - 1) & mer_mask;
			for(i=1; i < used_seqs; i++){
				mer_to_get = mhe[cur_seqs[i]];
				if(mer_to_get < 0){
					//Convert the cur_seqs[i] entry since negative implies reverse complement
					mer_to_get *= -1;
					mer_to_get += mhe.Length() - mer_size;
				}
				if(cur_mer != (sar_table[cur_seqs[i]]->GetMer(mer_to_get - 1) & mer_mask)){
					maxlen = 0;
					break;
				}
			}
		}
		//this stuff cleans up if there was a mismatch
		if(i < used_seqs){
			mhe.SetLength(mhe.Length() - jump_size);
			j--;
			for(;j >= 0; j--)
				if(mhe[cur_seqs[j]] >= 0)
					mhe.SetStart(cur_seqs[j], mhe[cur_seqs[j]] + jump_size);
		}
		//Invert the sequence directions so that we extend in the other direction
		//next time through the loop.  The second time we do this we are setting
		//sequence directions back to normal.
		for(i = 0; i < used_seqs; i++)
			mhe.SetStart(cur_seqs[i], -mhe[cur_seqs[i]]);
		//if we've already been through twice then decrease the jump size
		if(directionI >= 1)
			jump_size = 1;
	}
	delete[] cur_seqs;
}

boolean HashingMatchFinder::MatchAmbiguities(MemHashEntry& mhe){

}
//mems which span the origin will be hashed a second time

boolean HashingMatchFinder::WarpExtend(MemHashEntry& mhe, gnSeqI max_backward, gnSeqI max_forward){
	uint64 cur_mer;
	uint64 mer_mask = sar_table[0]->GetMerMask();
	//which sequences are used in this match?
	uint32* cur_seqs = new uint32[seq_count];
	uint32 used_seqs = 0;
	for(uint32 seqI = 0; seqI < seq_count; seqI++){
		if(mhe[seqI] != MEM_NO_MATCH){
			cur_seqs[used_seqs] = seqI;
			used_seqs++;
		}
	}
	//First extend backwards then extend forwards.  The following loop does them both.
	gnSeqI jump_size = 1;
	for(uint32 directionI = 0; directionI < 2; directionI++){
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
			if(sar_table[cur_seqs[maxI]]->IsCircular()){
				if(sar_table[cur_seqs[maxI]]->Length() < maxlen)
					maxlen = sar_table[cur_seqs[maxI]]->Length();
			}else if(mhe[cur_seqs[maxI]] < 0){
				int64 rc_len = sar_table[cur_seqs[maxI]]->Length() - mhe.Length() -  mhe[cur_seqs[maxI]] - 1;
				if( rc_len < maxlen)
					maxlen = rc_len;
			}else if(mhe[cur_seqs[maxI]] - 1 < maxlen)
				maxlen = mhe[cur_seqs[maxI]] - 1;
		
		int32 j=0;
		uint32 i;
		uint64 last_mer = 0xFFFFFFFF;
		last_mer <<= 32;
		last_mer |= 0xFFFFFFFF;
		//DANGEROUS BUG!  If there is a mem of T's that is bigger than MER_SIZE
		//it won't be extended.
		boolean first_base = true;
		while(maxlen - jump_size >= 0){
			mhe.SetLength(mhe.Length() + jump_size);
			maxlen -= jump_size;
			for(j=0; j < used_seqs; j++)
				if(mhe[cur_seqs[j]] > 0){
					mhe.SetStart(cur_seqs[j], mhe[cur_seqs[j]] - jump_size);
					if(mhe[cur_seqs[j]] <= 0)
						mhe.SetStart(cur_seqs[j], mhe[cur_seqs[j]] + sar_table[cur_seqs[j]]->Length());
				}
			//check that all mers at the new position match
			int64 mer_to_get = mhe[cur_seqs[0]];
			if(mer_to_get < 0){
				mer_to_get *= -1;
				mer_to_get += mhe.Length() - mer_size;
			}
			cur_mer = sar_table[cur_seqs[0]]->GetMer(mer_to_get - 1) & mer_mask;
			for(i=1; i < used_seqs; i++){
				mer_to_get = mhe[cur_seqs[i]];
				if(mer_to_get < 0){
					//Convert the cur_seqs[i] entry since negative implies reverse complement
					mer_to_get *= -1;
					mer_to_get += mhe.Length() - mer_size;
				}
				if(cur_mer != (sar_table[cur_seqs[i]]->GetMer(mer_to_get - 1) & mer_mask)){
					maxlen = 0;
					break;
				}
			}
			if(last_mer <= cur_mer && maxlen > 0){
				delete[] cur_seqs;
				return false;
			}
				
			last_mer = cur_mer;
		}
		//this stuff cleans up if there was a mismatch
		if(i < used_seqs){
			mhe.SetLength(mhe.Length() - jump_size);
			j--;
			for(;j >= 0; j--)
				if(mhe[cur_seqs[j]] >= 0)
					mhe.SetStart(cur_seqs[j], mhe[cur_seqs[j]] + jump_size);
		}
		//Invert the sequence directions so that we extend in the other direction
		//next time through the loop.  The second time we do this we are setting
		//sequence directions back to normal.
		for(i = 0; i < used_seqs; i++)
			mhe.SetStart(cur_seqs[i], -mhe[cur_seqs[i]]);
	}
	delete[] cur_seqs;
	return true;
}

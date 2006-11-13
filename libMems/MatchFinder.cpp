#include "MatchFinder.h"

MatchFinder::MatchFinder(){
	mer_size = 0;
	seq_count = 0;
	ambiguity_tolerance = 0;
}

//make sure this calls the destructor on each element
MatchFinder::~MatchFinder(){
}

MatchFinder::MatchFinder(MatchFinder& mf){
	mer_size = mf.mer_size;
	seq_count = mf.seq_count;
	ambiguity_tolerance = mf.ambiguity_tolerance;

	sar_table = mf.sar_table;
	seq_table = mf.seq_table;
	seqid_table = mf.seqid_table;
	for(uint32 i=0; i < UINT8_MAX; i++)
		sarid_table[i] = mf.sarid_table[i];
}

boolean MatchFinder::AddSuffixArray( SuffixArray* sar, gnSequence* seq ){
	if(sar == NULL){
		DebugMsg("MatchFinder::AddSuffixArray: Error null parameter.");
		return false;
	}
	SarHeader header = sar->GetHeader();
	for(uint32 i = 0; i < seqid_table.size(); i++){
		if(header.id == seqid_table[i]){
			DebugMsg("MatchFinder::AddSequence: Error seqid already exists.");
			return false;
		}
	}

	//check for consistency between sequence length and suffix array lengths
	if(seq != NULL && seq->length() != sar->Length()){
		DebugMsg("MatchFinder::AddSuffixArray: Error mismatched sar and sequence length.");
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
	return true;
}

boolean MatchFinder::FindMatches(uint32 mer_mask_size, gnSeqI startI, gnSeqI endI){
	//picked a semi-arbitrary number for buffer size.
	uint32 MER_BUFFER_SIZE = 10000;
	vector<bmer>* mer_vector;	// stores arrays of bmers for each sar.
	vector<uint32> mer_index;   // stores the indexes of the current mers in mer_vector
	vector<uint32> mer_baseindex;   // stores the index in the suffixArray of each of the first mers in mer_vector
	list<idmer> cur_mers;	// stores the current mers.
	list<idmer> cur_match;	// stores the current matching mers.
	
	list<uint32> sar_hitlist;	// list of sars to replace
	
	for(uint32 maskI = 0; maskI < sar_table.size(); maskI++){
		if(mer_mask_size > sar_table[maskI]->MerSize()){
			DebugMsg("MatchFinder::FindMatches: Requested mask size is too big.");
			return false;
		}else
			sar_table[maskI]->SetMerMaskSize(mer_mask_size);
	}
	mer_size = mer_mask_size;
	uint64 mer_mask = sar_table[0]->GetMerMask();
	
	//allocate buffer space
	mer_vector = new vector<bmer>[sar_table.size()];
	if(mer_vector == NULL){
		DebugMsg("MatchFinder::FindMatches: Error out of memory.");
		return false;
	}
	
	//initialize the data structures
	for(uint32 n = 0; n < sar_table.size(); n++){
		mer_vector[n].reserve(MER_BUFFER_SIZE);
		sar_table[n]->Read(mer_vector[n], MER_BUFFER_SIZE, 0);
		mer_index.push_back(0);
		mer_baseindex.push_back(0);
		idmer newmer;
		newmer.position = mer_vector[n][0].position;
		newmer.mer = mer_vector[n][0].mer & mer_mask;
		newmer.id = sar_table[n]->GetID();
		cur_mers.push_back(newmer);  //cur_mers gets the first mer from each suffix array
	}

	cur_mers.sort(&idmer_lessthan);
	//loop while there is data to hash.
	while(cur_mers.size() > 0){

		//scan through cur_mers looking for matches
		list<idmer>::iterator iter1 = cur_mers.begin();
		list<idmer>::iterator iter2 = cur_mers.begin();
		iter2++;
		if(cur_match.size() > 0){
			if(cur_match.begin()->mer < iter1->mer){
				//can't extend the match further.  Hash it.
				EnumerateMatches(cur_match);
				cur_match.clear();
			}else if(iter1->mer != iter2->mer || iter2 == cur_mers.end())
				//iter1 matches but iter2 doesnt.  make sure iter1 gets added
				cur_match.push_back(*iter1);
		}

		while(iter1->mer == iter2->mer && iter2 != cur_mers.end()){
//			cout << debugCount++ << " Match at " << iter1->position << '\t' << iter2->position << '\n';
			cur_match.push_back(*iter1);
			sar_hitlist.push_back(sarid_table[(uint8)iter1->id]);
			list<idmer>::iterator to_del = iter1;
			iter1++;
			iter2++;
			cur_mers.erase(to_del);
		}
		if(cur_match.size() > 0)	// there was a match.  get the last mer
			cur_match.push_back(*iter1);
		sar_hitlist.push_back(sarid_table[(uint8)iter1->id]);
		cur_mers.erase(iter1);
		
		//now replenish the cur_matches list 
		while(sar_hitlist.size() > 0){
			uint32 sar_hit = *sar_hitlist.begin();
			sar_hitlist.erase(sar_hitlist.begin());
			mer_index[sar_hit]++;
			
			//check if we need to read more data
			if(mer_index[sar_hit] >= mer_vector[sar_hit].size()){
				//out of data, must read more.
				mer_baseindex[sar_hit] += mer_vector[sar_hit].size();
				sar_table[sar_hit]->Read(mer_vector[sar_hit], MER_BUFFER_SIZE, mer_baseindex[sar_hit]);
				mer_index[sar_hit] = 0;
			}
			if(mer_vector[sar_hit].size() == 0){
				sar_hit++;
				continue;
			}
			//insert the new mer into the cur_mers list
			iter1 = cur_mers.begin();
			bmer merle = mer_vector[sar_hit][mer_index[sar_hit]];
			while(iter1->mer < merle.mer && iter1 != cur_mers.end())
				iter1++;
			idmer newmer;
			newmer.position = merle.position;
			newmer.mer = merle.mer & mer_mask;
			newmer.id = sar_table[sar_hit]->GetID();
			cur_mers.insert(iter1, newmer);
		}
	}
	delete[] mer_vector;
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

void MatchFinder::ExtendMatch(MemHashEntry& mhe, gnSeqI max_backward, gnSeqI max_forward){
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
//					if(!MatchAmbiguities(mhe, jump_size)){
						maxlen = 0;
						break;
//					}
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
boolean MatchFinder::MatchAmbiguities(MemHashEntry& mhe, uint32 match_size){
	if(ambiguity_tolerance == 0)
		return false;
			//check that all mers at the new position match
	//which sequences are used in this match?
	uint32* cur_seqs = new uint32[seq_count];
	uint32 used_seqs = 0;
	for(uint32 seqI = 0; seqI < seq_count; seqI++){
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
boolean MatchFinder::WarpExtend(MemHashEntry& mhe, gnSeqI max_backward, gnSeqI max_forward){
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

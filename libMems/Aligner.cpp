#include "Aligner.h"

int64 abs( int64 a ){
	return a < 0 ? -a : a;
}

Aligner::Aligner( int64 offset_tolerance, double LCB_minimum_density, double LCB_minimum_range ){
	debug = false;
	this->offset_tolerance = offset_tolerance;
	this->LCB_minimum_density = LCB_minimum_density;
	this->LCB_minimum_range = LCB_minimum_range;
}

/**
 * Bob's Cheesy LCB algorithm.  
 * This vehicle must be driven by a domain expert
 * This algorithm is a two step process:
 * First out of place MUMs are removed.
 * This is accomplished by sorting the MUMs on the first sequence
 * and tossing out any MUM whose generalized offset is more than "offset_threshold"
 * away from its neighbors.
 * Each MUM which passes the filter is labeled with a number starting at 0
 * and increasing by one for each MUM.
 * The second step involves sorting the labeled MUMs on their start positions in
 * the other sequences and looking for breakpoints.  Breakpoints are encountered
 * when a label changes by more than 1 between neighboring MUMs.  
 * The label of the last MUM in the ongoing block is recorded in the breakpoint set
 * After searching for breakpoints in each sequence the set of breakpoints is
 * returned.
 */
void Aligner::BobsCheesyLCB( MemList& mlist, int64 offset_threshold, set<uint>& breakpoints ){
	if( mlist.mem_list.size() == 0 )
		return;
	// can only look for breakpoints if there is more than one match!!
	if( mlist.mem_list.size() == 1 ){
		breakpoints.insert( 0 );
		return;
	}

	mlist.mem_list.sort( MemStartComparator( 0 ) );
	list<MemHashEntry*>::iterator mem_iter = mlist.mem_list.begin();
	list<MemHashEntry*>::iterator prev_iter;
	list<MemHashEntry*>::iterator next_iter;
	list<MemHashEntry*>::iterator last_iter = mlist.mem_list.end();
	list<LabeledMem> pair_list;
	uint cur_memI = 0;
	last_iter--;

	LabeledMem first_mem;
	first_mem.mem = *mem_iter;
	first_mem.label = cur_memI++;
	pair_list.push_back( first_mem );
	for(mem_iter++; mem_iter != last_iter; mem_iter++ ){
		prev_iter = mem_iter;
		prev_iter--;
		next_iter = mem_iter;
		next_iter++;
// aaron's criteria commented out...
//		if( abs( (*next_iter)->Offset() - (*prev_iter)->Offset()) < offset_threshold ){
//			if( abs( (*next_iter)->Offset() - (*mem_iter)->Offset()) >= offset_threshold
//			 || abs( (*mem_iter)->Offset() - (*prev_iter)->Offset()) >= offset_threshold ){
// bob's criteria in:
		if( offset_threshold != -1 ){
			if( abs( (*next_iter)->Offset() - (*mem_iter)->Offset()) >= offset_threshold
			 && abs( (*mem_iter)->Offset() - (*prev_iter)->Offset()) >= offset_threshold ){
				// it doesn't conform.  trash it.
				list<MemHashEntry*>::iterator to_del = mem_iter;
				mem_iter--;
				mlist.mem_list.erase( to_del );
				continue;
			}
		}
//		}
		// add this one to the list.
		LabeledMem lm;
		lm.mem = *mem_iter;
		lm.label = cur_memI++;
		pair_list.push_back( lm );
	}
	LabeledMem last_mem;
	last_mem.mem = *mem_iter;
	last_mem.label = cur_memI++;
	pair_list.push_back( last_mem );
	
	breakpoints.clear();
	for( uint seqI = 1; seqI < seq_count; seqI++ ){
		// sort the list on the current genome
		LabeledMemComparator lmc( seqI );
		pair_list.sort( lmc );
		list<LabeledMem>::iterator pair_iter = pair_list.begin();
		uint block_start = pair_iter->label;
		for( pair_iter++; pair_iter != pair_list.end(); pair_iter++ ){
			list<LabeledMem>::iterator pair_prev = pair_iter;
			pair_prev--;
			if( pair_iter->mem->Start( seqI ) < 0 &&
				pair_prev->mem->Start( seqI ) < 0 &&
				pair_prev->label - pair_iter->label == 1 ){
					continue;
			}
			if( pair_iter->mem->Start( seqI ) > 0 &&
				pair_prev->mem->Start( seqI ) > 0 &&
				pair_iter->label - pair_prev->label == 1 ){
					continue;
			}
			// since it didn't meet either of the above
			// criteria it's a breakpoint.  insert the label of the end of the current block
			// note that if it's a reverse complement block, the end label is really the start label
			if( pair_prev->mem->Start( seqI ) < 0 )
				breakpoints.insert( block_start );
			else
				breakpoints.insert( pair_prev->label );
			block_start = pair_iter->label;
		}
		if( pair_iter != pair_list.begin() ){
			pair_iter--;
			breakpoints.insert( pair_iter->label );
		}
	}
}

void Aligner::FilterLCBs( MemList& meml, set<uint>& breakpoints, vector<MemList>& lcb_list, double LCB_minimum_density, double LCB_minimum_range, uint LCB_minimum_matches ){
	// there must be at least one end of a block defined
	if( breakpoints.size() < 1 )
		return;
	// organize the LCBs into different MemList instances
	list<MemHashEntry*>::iterator mem_iter;

	vector<MemHashEntry*> mem_list;
	mem_list.reserve( meml.mem_list.size() );
	mem_list.insert( mem_list.begin(), meml.mem_list.begin(), meml.mem_list.end() );
	meml.mem_list.clear();
	lcb_list.clear();
	set<uint>::iterator break_iter = breakpoints.begin();
	uint prev_break = 0;	// prev_break is the first match in the current block

	for( ; break_iter != breakpoints.end(); break_iter++ ){
		if( *break_iter - prev_break + 1 < LCB_minimum_matches ){
			prev_break = *break_iter + 1;
			continue;
		}

		// add the mems in this interval to a new MemList
		MemList lcb = meml;
		lcb.mem_list.clear();
		lcb.mem_list.insert( lcb.mem_list.begin(), mem_list.begin() + prev_break, mem_list.begin() + *break_iter + 1 );
		prev_break = *break_iter + 1;

		// check minimum density and range
		gnSeqI coverage = 0;
		mem_iter = lcb.mem_list.begin();
		while( mem_iter != lcb.mem_list.end() ){
			coverage += (*mem_iter)->Length();
			mem_iter++;
		}

		boolean big_enough = true;
		boolean dense_enough = true;
//		boolean big_enough = false;
//		boolean dense_enough = false;
		for( uint seqI = 0; seqI < seq_count; seqI++ ){
			mem_iter = lcb.mem_list.end();
			mem_iter--;
			int64 end_pos = (*mem_iter)->End( seqI ) + 1;
			mem_iter = lcb.mem_list.begin();
			int64 start_pos = (*mem_iter)->Start( seqI );
			if( abs( end_pos - start_pos ) < LCB_minimum_range )
				big_enough = false;
//			if( abs( end_pos - start_pos ) >= LCB_minimum_range )
//				big_enough = true;

			double density = (double) coverage / (double)(abs( end_pos - start_pos ));
			if( density < LCB_minimum_density )
				dense_enough = false;
//			if( density >= LCB_minimum_density )
//				dense_enough = true;
		}
		if( !big_enough || !dense_enough )
			continue;

		// add the new MemList to the set.
		lcb_list.push_back( lcb );
	}
}

void Aligner::align( MemList& mlist, vector<MemList>& LCB_list, uint min_match_size ){
	if( mlist.mem_list.size() == 0 ){
		Throw_gnExMsg( AlignerError(), "ZERO size mem list!" );
	}

	ml = mlist;
	seq_count = (*ml.mem_list.begin())->SeqCount();

	// now that matches are created, filter out all subset matches.
	cout << "Starting with " << ml.mem_list.size() << " MUMs\n";

	// filter out subset matches.
	ml.MultiplicityFilter( seq_count );

	// filter out small matches
	list<MemHashEntry*>::iterator mem_iter;
	list<MemHashEntry*>::iterator to_del;
	mem_iter = ml.mem_list.begin();
	while( mem_iter != ml.mem_list.end() ){
		if( (*mem_iter)->Length() < min_match_size ){
			// delete this one.
			to_del = mem_iter;
			mem_iter++;
			ml.mem_list.erase( to_del );
		}else
			mem_iter++;
	}
	
	
	cout << "There are " << ml.mem_list.size() << " " << seq_count << "-way MUMs longer than " << min_match_size << " b.p.\n";
	
	if( ml.mem_list.size() == 0 ){
		Throw_gnExMsg( AlignerError(), "ZERO size mem list!" );
	}
	// now do Bob's cheesy LCB algorithm
	set<uint> breakpoints;
	BobsCheesyLCB( ml, offset_tolerance, breakpoints );

	// organize the LCBs into different MemList instances
	FilterLCBs( ml, breakpoints, LCB_list, LCB_minimum_density, LCB_minimum_range, 0 );
	
	ml.mem_list.clear();
	for( uint lcbI = 0; lcbI < LCB_list.size(); lcbI++ ){
		ml.mem_list.insert( ml.mem_list.end(), LCB_list[ lcbI ].mem_list.begin(), LCB_list[ lcbI ].mem_list.end() );
	}

	LCB_list.clear();
	breakpoints.clear();
	BobsCheesyLCB( ml, offset_tolerance, breakpoints );

	cout << "Filtering leaves " << ml.mem_list.size() << " MUMs\n";
	cout << "And " << breakpoints.size() << " breakpoints\n";

	FilterLCBs( ml, breakpoints, LCB_list, LCB_minimum_density, LCB_minimum_range, 2 );
}

void Aligner::writeLCBlist( const vector<MemList>& LCB_list, ostream& lcb_out ){
	if( LCB_list.size() == 0 )
		return;
	
	const list<MemHashEntry*>& mem_list = LCB_list[0].mem_list;
	
	MemHashEntry* first_mem = *(mem_list.begin());
	unsigned int seq_count = first_mem->SeqCount();
	unsigned int mer_size = first_mem->MerSize();

	lcb_out << "FormatVersion" << '\t' << 2 << "\n";
	lcb_out << "SequenceCount" << '\t' << seq_count << "\n";
	lcb_out << "MerSize" << '\t' << mer_size << "\n";
	for(unsigned int seqI = 0; seqI < seq_count; seqI++){
		lcb_out << "Sequence" << seqI << "File" << '\t';
		if( LCB_list[0].seq_filename.size() > seqI )
			lcb_out << LCB_list[0].seq_filename[seqI];
		else
			lcb_out << "null";
		lcb_out << "\n";
		lcb_out << "Sequence" << seqI << "Length" << '\t';
		if( LCB_list[0].seq_table.size() > seqI )
			lcb_out << LCB_list[0].seq_table[seqI]->length();
		else
			lcb_out << "0";
		lcb_out << "\n";
	}

	uint match_count = 0;
	uint lcbI;
	for( lcbI = 0; lcbI < LCB_list.size(); lcbI++ ){
		match_count += LCB_list[ lcbI ].mem_list.size();
	}
	lcb_out << "MatchCount" << '\t' << match_count << endl;

	for( lcbI = 0; lcbI < LCB_list.size(); lcbI++ ){
		//get all the mems out of the lists and write them out
	    list<MemHashEntry*>::const_iterator mem_iter;
		mem_iter = LCB_list[ lcbI ].mem_list.begin();
		for(; mem_iter != LCB_list[ lcbI ].mem_list.end(); mem_iter++)
			lcb_out << lcbI << '\t' << **mem_iter << "\n";
	}
}

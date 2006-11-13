#include "MatchList.h"
#include "DNAFileSML.h"
#include "MemHash.h"
#include <map>
#include <sstream>
#include "wx/timer.h"

MatchList::MatchList( const MatchList& ml ){
	*this = ml;
}

MatchList& MatchList::operator=( const MatchList& ml ){
	match_list = ml.match_list;
	sml_filename = ml.sml_filename;
	seq_filename = ml.seq_filename;
	sml_table = ml.sml_table;
	seq_table = ml.seq_table;
	return *this;
}

void MatchList::LoadSequences( uint mer_size, ostream* log_stream ){
	if( seq_filename.size() != sml_filename.size() ){
		Throw_gnExMsg( InvalidArgument(), "The number of sequence filenames and sorted mer list filenames must be equal." );
	}
	
	if( seq_filename.size() == 0 )
		return;

	for( uint seqI = 0; seqI < seq_filename.size(); seqI++ ){
		gnSequence* file_sequence = new gnSequence();
		// Load the sequence and tell the user if it loaded successfully
		try{
			file_sequence->LoadSource( seq_filename[ seqI ] );
		}catch( gnException& gne ){
			delete file_sequence;
			if( gne.GetCode() == FileNotOpened() )
				cout << "Error loading file.\n";
			else
				cout << gne;
			return;
		}
		
		seq_table.push_back( file_sequence );
		if( log_stream != NULL ){
			(*log_stream) << "Sequence loaded successfully.\n";
			(*log_stream) << seq_filename[ seqI ] << " " << file_sequence->length() << " base pairs.\n";
		}

		// define a DNAFileSML to store a sorted mer list
		DNAFileSML* file_sml = new DNAFileSML();
		boolean success = true;
		try{
			file_sml->LoadFile( sml_filename[ seqI ] );
		}catch( gnException& gne ){
			success = false;
		}
		boolean recreate = false;
		if(success && (file_sml->MerSize() != mer_size )){
			if( log_stream != NULL )
				(*log_stream) << "Mer size mismatch.  A new sorted mer list will be created.\n";
			recreate = true;
		}
		if(!success || recreate){
			if( log_stream != NULL )
				(*log_stream) << "Creating sorted mer list\n";
			long start_time = wxGetLocalTime();
			file_sml->Create( *file_sequence, mer_size );
			long end_time = wxGetLocalTime();
		 	if( log_stream != NULL )
				(*log_stream) << "Create time was: " << end_time - start_time << " seconds.\n";
		}else if( log_stream != NULL )
			(*log_stream) << "Sorted mer list loaded successfully\n";
		
		sml_table.push_back( file_sml );
	}
}

void MatchList::ReadList(istream& match_file){
	string tag;
	gnSeqI len;
	int64 start;
	unsigned int seq_count;
	
	match_file >> tag;	//format version tag
	if( tag != "FormatVersion" ){
		Throw_gnEx(InvalidFileFormat());
	}
	match_file >> tag;	//format version
	if( tag != "3" ){
		Throw_gnEx(InvalidFileFormat());
	}
	match_file >> tag;	//sequence count tag
	if( tag != "SequenceCount" ){
		Throw_gnEx(InvalidFileFormat());
	}
	match_file >> seq_count;	//sequence count
	if(seq_count < 2){
		Throw_gnEx(InvalidFileFormat());
	}
	
	// read the sequence file names and lengths
	for( unsigned int seqI = 0; seqI < seq_count; seqI++ ){
		match_file >> tag;	// name tag
		match_file >> tag;	// name
		seq_filename.push_back(tag);
//		try{
//			gnSequence *new_seq = new gnSequence();
//			new_seq->LoadSource(tag);
//			seq_table.push_back( new_seq );
//		}catch( gnException& gne );
		match_file >> tag;	// length tag
		match_file >> tag;	// length
	}

	// read the number of matches
	unsigned int match_count;
	match_file >> tag;	// match count tag
	match_file >> match_count;	// match count
	
	// read the matches
	map< MatchID_t, Match* > match_map;
	string cur_line;
	getline( match_file, cur_line );
	while( getline( match_file, cur_line ) ){
		Match mhe( seq_count );
		stringstream line_stream( cur_line );
		
		line_stream >> len;
		mhe.SetLength(len);

		for(uint32 seqI = 0; seqI < seq_count; seqI++){
			line_stream >> start;
			mhe.SetStart(seqI, start);
		}
		
		mhe.CalculateOffset();
		
		MatchID_t match_id;
		line_stream >> match_id;
		
		uint sub_count;
		boolean bad_stream = false;
		line_stream >> sub_count;
		for( uint subI = 0; subI < sub_count; subI++ ){
			//break if the stream ended early
			if(!line_stream.good() ){
				bad_stream = true;
				break;
			}
			MatchID_t sub;
			line_stream >> sub;
			mhe.AddSubset( sub );
		}

		if( bad_stream )
			break;

		uint sup_count;
		line_stream >> sup_count;
		for( uint supI = 0; supI < sup_count; supI++ ){
			//break if the stream ended early
			if(!line_stream.good() ){
				bad_stream = true;
				break;
			}
			MatchID_t sup;
			line_stream >> sup;
			mhe.AddSuperset( sup );
		}
		if( bad_stream )
			break;

		Match* new_match = new Match( mhe );
		match_list.push_back( new_match );
		match_map.insert( map< MatchID_t, Match* >::value_type( match_id, new_match ));
	}
	if( match_count != match_list.size() ){
		Throw_gnEx(InvalidFileFormat());
	}
	
	// now remap the subset and superset links
	list<Match*>::iterator match_iter = match_list.begin();
	map<MatchID_t, Match*>::iterator map_iter;
	for(; match_iter != match_list.end(); match_iter++ ){
		// remap all subsets
		set<MatchID_t>& subsets = (*match_iter)->Subsets();
		set<MatchID_t> new_subsets;
		set<MatchID_t>::iterator sub_iter = subsets.begin();
		for(; sub_iter != subsets.end(); sub_iter++ ){
			map_iter = match_map.find( *sub_iter );
			new_subsets.insert( map_iter->second->MatchID() );
		}
		subsets = new_subsets;

		// remap all supersets
		set<MatchID_t>& supersets = (*match_iter)->Supersets();
		set<MatchID_t> new_supersets;
		set<MatchID_t>::iterator super_iter = supersets.begin();
		for(; super_iter != supersets.end(); super_iter++ ){
			map_iter = match_map.find( *super_iter );
			new_supersets.insert( map_iter->second->MatchID() );
		}
		supersets = new_supersets;
	}
}

void MatchList::WriteList(ostream& match_file) const{
	if( match_list.size() == 0 )
		return;
	Match* first_mem = *(match_list.begin());
	unsigned int seq_count = first_mem->SeqCount();

	match_file << "FormatVersion" << '\t' << 3 << "\n";
	match_file << "SequenceCount" << '\t' << seq_count << "\n";
	for(unsigned int seqI = 0; seqI < seq_count; seqI++){
		match_file << "Sequence" << seqI << "File" << '\t';
		if( seq_filename.size() > seqI )
			match_file << seq_filename[seqI];
		else
			match_file << "null";
		match_file << "\n";
		match_file << "Sequence" << seqI << "Length" << '\t';
		if( seq_table.size() > seqI )
			match_file << seq_table[seqI]->length();
		else
			match_file << "0";
		match_file << "\n";
	}

	match_file << "MatchCount" << '\t' << match_list.size() << endl;

	//get all the mems out of the hash table and write them out
    list<Match*>::const_iterator match_iter;
	match_iter = match_list.begin();
	set<MatchID_t> cur_set;
	set<MatchID_t>::iterator set_iter;
	for(; match_iter != match_list.end(); match_iter++){
		// print the match
		match_file << **match_iter << '\t';

		// print the match id
		match_file << (*match_iter)->MatchID() << '\t';
		
		// print subset id's
		cur_set = (*match_iter)->Subsets();
		match_file << cur_set.size();
		set_iter = cur_set.begin();
		for(; set_iter != cur_set.end(); set_iter++ ){
			match_file << '\t' << *set_iter;
		}

		// print superset id's
		cur_set = (*match_iter)->Supersets();
		match_file << '\t' << cur_set.size();
		set_iter = cur_set.begin();
		for(; set_iter != cur_set.end(); set_iter++ ){
			match_file << '\t' << *set_iter;
		}
		match_file << endl;
	}
}

void MatchList::MultiplicityFilter( unsigned mult ){
	list<Match*>::iterator match_iter;
	list<Match*>::iterator to_del;
	map< MatchID_t, Match* > id_map;
	match_iter = match_list.begin();
	while( match_iter != match_list.end() ){
		id_map.insert( map< MatchID_t, Match* >::value_type( (*match_iter)->MatchID(), *match_iter ) );
		match_iter++;
	}

	match_iter = match_list.begin();
	while( match_iter != match_list.end() ){
		if( (*match_iter)->Multiplicity() != mult ){
			MatchID_t cur_id = (*match_iter)->MatchID();
			// delete sub/superset links
			set<MatchID_t>& subsets = (*match_iter)->Subsets();
			set<MatchID_t>::iterator s_iter = subsets.begin();
			for( ; s_iter != subsets.end(); s_iter++ ){
				map< MatchID_t, Match* >::iterator s_map = id_map.find( *s_iter );
				s_map->second->UnlinkSuperset( cur_id );
			}

			set<MatchID_t>& supersets = (*match_iter)->Supersets();
			s_iter = supersets.begin();
			for( ; s_iter != supersets.end(); s_iter++ ){
				map< MatchID_t, Match* >::iterator s_map = id_map.find( *s_iter );
				s_map->second->UnlinkSubset( cur_id );
			}

			map< MatchID_t, Match* >::iterator map_del = id_map.find( cur_id );
			id_map.erase( map_del );

			// delete this one.
			to_del = match_iter;
			match_iter++;
			match_list.erase( to_del );
		}else
			match_iter++;
	}
}

void MatchList::ExactFilter( valarray<bool>& filter_spec ){
	list<Match*>::iterator match_iter;
	list<Match*>::iterator to_del;
	match_iter = match_list.begin();
	while( match_iter != match_list.end() ){
		uint64 matchnumber = (*match_iter)->MatchNumber();
		valarray<bool> matchnum( false, (*match_iter)->SeqCount());
		for( uint32 seqI = (*match_iter)->SeqCount(); seqI > 0; seqI-- ){
			if( matchnumber & 0x1 )
				matchnum[ seqI - 1 ] = true;
			matchnumber >>= 1;
		}

		valarray<bool> equal = matchnum == filter_spec;
		unsigned msum = equal.min();
		if( msum == false ){
			// delete this one.
			to_del = match_iter;
			match_iter++;
			match_list.erase( to_del );
		}else
			match_iter++;
	}
}

void MatchList::IntersectFilter( valarray<bool>& filter_spec ){
	list<Match*>::iterator match_iter;
	list<Match*>::iterator to_del;
	match_iter = match_list.begin();
	while( match_iter != match_list.end() ){
		uint64 matchnumber = (*match_iter)->MatchNumber();
		valarray<bool> matchnum( false, (*match_iter)->SeqCount());
		for( uint32 seqI = (*match_iter)->SeqCount(); seqI > 0; seqI-- ){
			if( matchnumber & 0x1 )
				matchnum[ seqI - 1 ] = true;
			matchnumber >>= 1;
		}

		matchnum &= filter_spec;
		unsigned msum = matchnum.max();
		if( msum == 0 ){
			// delete this one.
			to_del = match_iter;
			match_iter++;
			match_list.erase( to_del );
		}else
			match_iter++;
	}
}

void MatchList::UnlinkedFirstStartFilter( unsigned start_seq )
{
	list<Match*>::iterator match_iter;
	list<Match*>::iterator to_del;
	match_iter = match_list.begin();
	while( match_iter != match_list.end() ){
		if( (*match_iter)->FirstStart() != start_seq ||
			(*match_iter)->Supersets().size() > 0 ) 
		{
			to_del = match_iter;
			match_iter++;
			match_list.erase( to_del );
			continue;
		}
		match_iter++;
	}
}

void MatchList::ToMemList( list<MemHashEntry*>& mem_list ) const{
	vector<MemHashEntry*> mem_array;
	mem_array.reserve( match_list.size() );

	list<Match*>::const_iterator match_iter = match_list.begin();
	map< MatchID_t, uint > id_map;
	uint matchI = 0;
	for( ; match_iter != match_list.end(); match_iter++ ){
		mem_array.push_back( new MemHashEntry( **match_iter ) );
		id_map.insert( map< MatchID_t, uint >::value_type( (*match_iter)->MatchID(), matchI ) );
		matchI++;
	}

	matchI = 0;
	for( match_iter = match_list.begin() ; match_iter != match_list.end(); match_iter++ ){
		set< MatchID_t >& subsets = (*match_iter)->Subsets();
		set< MatchID_t >::iterator sub_iter = subsets.begin();
		for(; sub_iter != subsets.end(); sub_iter++ ){
			map< MatchID_t, uint >::iterator id_iter = id_map.find( *sub_iter );
			mem_array[ matchI ]->LinkSubset( mem_array[ id_iter->second ] );
		}
		matchI++;
	}
	
	mem_list.insert( mem_list.begin(), mem_array.begin(), mem_array.end() );
}


void MatchList::FromMemList( const list<MemHashEntry*>& mem_list ) {
	
	match_list.clear();
	map< MemHashEntry*, Match* > match_map;
	list<MemHashEntry*>::const_iterator mem_iter = mem_list.begin();
	for( ; mem_iter != mem_list.end(); mem_iter++ ){
		Match* new_match = new Match( **mem_iter );
		match_map.insert( map< MemHashEntry*, Match* >::value_type( *mem_iter, new_match ) );
		match_list.push_back( new_match );
	}

	map< MemHashEntry*, Match* >::iterator match_iter = match_map.begin();
	for( ; match_iter != match_map.end(); match_iter++ ){
		set<MemHashEntry*> subsets = match_iter->first->GetSubsets();
		set<MemHashEntry*>::iterator sub_iter = subsets.begin();
		for( ; sub_iter != subsets.end(); sub_iter++ ){
			map< MemHashEntry*, Match* >::iterator sub_match = match_map.find( *sub_iter );
			match_iter->second->LinkSubset( sub_match->second );
		}
	}
}

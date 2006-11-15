/*******************************************************************************
 * $Id: MatchList.cpp,v 1.22 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MatchList.h"
#include "libMems/DNAFileSML.h"
#include "libMems/DNAMemorySML.h"
#include "libMems/MemHash.h"
#include <map>
#include <sstream>
#include <ctime>

using namespace std;
using namespace genome;
namespace mems {

typedef void* MatchID_t;

MatchList::MatchList( const MatchList& ml ){
	*this = ml;
}

MatchList& MatchList::operator=( const MatchList& ml ){
	vector< Match* >::operator=( ml );
	sml_filename = ml.sml_filename;
	seq_filename = ml.seq_filename;
	sml_table = ml.sml_table;
	seq_table = ml.seq_table;
	return *this;
}

void MatchList::LoadSequences( ostream* log_stream ){
	
	if( seq_filename.size() == 0 )
		return;

	gnSeqI total_len = 0;
	for( uint seqI = 0; seqI < seq_filename.size(); seqI++ ){
		gnSequence* file_sequence = new gnSequence();
		// Load the sequence and tell the user if it loaded successfully
		try{
			file_sequence->LoadSource( seq_filename[ seqI ] );
		}catch( gnException& gne ){
			delete file_sequence;
			if( gne.GetCode() == FileNotOpened() )
				cerr << "Error loading " << seq_filename[ seqI ] << endl;
			else
				cerr << gne;
			return;
		}catch( exception& e ){
			delete file_sequence;
			cerr << "Unhandled exception loading " << seq_filename[ seqI ] << endl;
			cerr << "At: " << __FILE__ << ":" << __LINE__ << endl;
			cerr << e.what();
			return;
		}catch( ... ){
			delete file_sequence;
			cerr << "Unknown exception when loading " << seq_filename[ seqI ] << endl;
			return;
		}
		
		total_len += file_sequence->length();
		seq_table.push_back( file_sequence );
		if( log_stream != NULL ){
			(*log_stream) << "Sequence loaded successfully.\n";
			(*log_stream) << seq_filename[ seqI ] << " " << file_sequence->length() << " base pairs.\n";
		}
	}

}

void MatchList::LoadSMLs( uint mer_size, ostream* log_stream ){

	// if the mer_size parameter is 0 then calculate a default mer size for these sequences
	if( mer_size == 0 ){
		mer_size = GetDefaultMerSize( seq_table );
		if( log_stream != NULL ){
			(*log_stream) << "Using weight " << mer_size << " mers for initial seeds\n";
		}
	}

	// load and creates SMLs as necessary
	uint64 default_seed = getSeed( mer_size );
	vector< uint > create_list;
	uint seqI = 0;
	for( seqI = 0; seqI < seq_table.size(); seqI++ ){
		// define a DNAFileSML to store a sorted mer list
		DNAFileSML* file_sml = new DNAFileSML();
		sml_table.push_back( file_sml );

		boolean success = true;
		try{
			file_sml->LoadFile( sml_filename[ seqI ] );
		}catch( gnException& gne ){
			success = false;
			create_list.push_back( seqI );
		}
		boolean recreate = false;
		if(success && (file_sml->Seed() != default_seed )){
			if( log_stream != NULL )
				(*log_stream) << "Default seed mismatch.  A new sorted mer list will be created.\n";
			recreate = true;
			create_list.push_back( seqI );
		}

		if( success && !recreate && log_stream != NULL )
			(*log_stream) << "Sorted mer list loaded successfully\n";
	}

	// free up memory before creating any SMLs
	if( create_list.size() > 0 )
		for( seqI = 0; seqI < sml_table.size(); seqI++ ){
			sml_table[ seqI ]->Clear();
			delete sml_table[ seqI ];
			sml_table[ seqI ] = NULL;
		}
	
	// create any SMLs that need to be created
	for( uint createI = 0; createI < create_list.size(); createI++ ){
		if( log_stream != NULL )
			(*log_stream) << "Creating sorted mer list\n";
		try{

		time_t start_time = time(NULL);
		sml_table[ create_list[ createI ] ] = new DNAFileSML( sml_filename[ create_list[ createI ] ] );
		sml_table[ create_list[ createI ] ]->Create( *seq_table[ create_list[ createI ] ], default_seed );
		time_t end_time = time(NULL);
	 	if( log_stream != NULL )
			(*log_stream) << "Create time was: " << end_time - start_time << " seconds.\n";
		
		}catch(...){
			cerr << "Error creating sorted mer list\n";
			throw;
		}
	}
	
	// reload the other SMLs now that creation has completed
	if( create_list.size() > 0 ){
		for( seqI = 0; seqI < seq_filename.size(); seqI++ ){
			if( sml_table[ seqI ] != NULL )
				continue;
			sml_table[ seqI ] = new DNAFileSML( sml_filename[ seqI ] );
			try{
				((DNAFileSML*)sml_table[ seqI ])->LoadFile( sml_filename[ seqI ] );
			}catch( gnException& gne ){
				cerr << "Error loading sorted mer list\n";
				throw;
			}
		}
	}
}

uint MatchList::GetDefaultMerSize( const vector< gnSequence* >& seq_table ){
	gnSeqI total_len = 0;
	for( uint seqI = 0; seqI < seq_table.size(); seqI++ )
		total_len += seq_table[ seqI ]->length();
	return getDefaultSeedWeight( total_len / seq_table.size() );
}

void MatchList::LoadMFASequences( const string& mfa_filename, uint mer_size, ostream* log_stream, boolean load_smls ) {
	gnSequence file_sequence;
	// Load the sequence and tell the user if it loaded successfully
	try{
		file_sequence.LoadSource( mfa_filename );
	}catch( gnException& gne ){
		if( gne.GetCode() == FileNotOpened() )
			cerr << "Error loading " << mfa_filename << endl;
		else
			cerr << gne;
		return;
	}catch( exception& e ){
		cerr << "Unhandled exception loading " << mfa_filename << endl;
		cerr << "At: " << __FILE__ << ":" << __LINE__ << endl;
		cerr << e.what();
		return;
	}catch( ... ){
		cerr << "Unknown exception when loading " << mfa_filename << endl;
		return;
	}

	seq_filename.clear();
	gnSeqI total_len = 0;
	for( uint contigI = 0; contigI < file_sequence.contigListSize(); contigI++ ){
		gnSequence* contig_seq = new gnSequence( file_sequence.contig( contigI ) );
		seq_filename.push_back( mfa_filename );
//		seq_filename.push_back( file_sequence.contigName( contigI ) );
		if( log_stream != NULL ){
			(*log_stream) << "Sequence loaded successfully.\n";
			(*log_stream) << seq_filename[ contigI ] << " " << contig_seq->length() << " base pairs.\n";
		}
		seq_table.push_back( contig_seq );
	}
	// if the mer_size parameter is 0 then calculate a default mer size for these sequences
	if( mer_size == 0 ){
		mer_size = GetDefaultMerSize( seq_table );
		if( log_stream != NULL ){
			(*log_stream) << "Using " << mer_size << "-mers for initial seeds\n";
		}
	}

	uint64 default_seed = getSeed( mer_size );

	for( uint contigI = 0; contigI < file_sequence.contigListSize(); contigI++ ){
		// define a DNAMemorySML to store a sorted mer list
		if( load_smls ){
			DNAMemorySML* contig_sml = new DNAMemorySML();
			boolean success = true;
			if( log_stream != NULL )
				(*log_stream) << "Creating sorted mer list\n";
			time_t start_time = time(NULL);
			contig_sml->Create( *seq_table[contigI], default_seed );
			time_t end_time = time(NULL);
		 	if( log_stream != NULL )
				(*log_stream) << "Create time was: " << end_time - start_time << " seconds.\n";
			
			sml_table.push_back( contig_sml );
		}
	}
}

void MatchList::Clear() {
	for( uint seqI = 0; seqI < seq_table.size(); seqI++ ){
		if( seq_table[ seqI ] != NULL )
			delete seq_table[ seqI ];
	}
	for( uint seqI = 0; seqI < sml_table.size(); seqI++ ){
		if( sml_table[ seqI ] != NULL )
			delete sml_table[ seqI ];
	}
	vector<Match*>::iterator match_iter = begin();
	for(; match_iter != end(); match_iter++ ){
		(*match_iter)->Free();
		(*match_iter) = NULL;
	}
	seq_table.clear();
	sml_table.clear();
	clear();
	seq_filename.clear();
	sml_filename.clear();
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
		getline( match_file, tag );
		// skip the tab character
		tag = tag.substr( 1 );
		seq_filename.push_back(tag);
//		try{
//			gnSequence *new_seq = new gnSequence();
//			new_seq->LoadSource(tag);
//			seq_table.push_back( new_seq );
//		}catch( gnException& gne );
		match_file >> tag;	// length tag
		gnSeqI seq_len;
		match_file >> seq_len;	// length
		if( seqI < seq_table.size() )
			if( seq_table[ seqI ]->length() != seq_len ){
				cerr << "Warning: Genome sizes in the match list differ.\n";
				cerr << "seq_table[ " << seqI << " ]->length() " << seq_table[ seqI ]->length() << " seq_len: " << seq_len << endl;
			}
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
			mhe.AddSubset( (Match*)sub );
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
			mhe.AddSuperset( (Match*)sup );
		}
		if( bad_stream )
			break;
		
		Match* new_match = mhe.Copy();
		push_back( new_match );
		match_map.insert( map< MatchID_t, Match* >::value_type( match_id, new_match ));
	}
	if( match_count != size() ){
		Throw_gnEx(InvalidFileFormat());
	}
	
	// now remap the subset and superset links
	RemapSubsetMatchAddresses( match_map, *this );

}

void MatchList::WriteList(ostream& match_file) const{
	if( size() == 0 )
		return;
	Match* first_mem = *(begin());
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

	match_file << "MatchCount" << '\t' << size() << endl;

	//get all the mems out of the hash table and write them out
	vector<Match*>::const_iterator match_iter;
	match_iter = begin();
	set<Match*> cur_set;
	set<Match*>::iterator set_iter;
	for(; match_iter != end(); match_iter++){
		// print the match
		match_file << **match_iter << '\t';

		// print the match address
		match_file << (MatchID_t)(*match_iter) << '\t';
		
		// print subset id's
		const set<Match*>& cur_set = (*match_iter)->Subsets();
		match_file << cur_set.size();
		set<Match*>::const_iterator set_iter = cur_set.begin();
		for(; set_iter != cur_set.end(); set_iter++ ){
			match_file << '\t' << (MatchID_t)*set_iter;
		}

		// print superset id's
		const set<Match*>& cur_set2 = (*match_iter)->Supersets();
		match_file << '\t' << cur_set2.size();
		set_iter = cur_set2.begin();
		for(; set_iter != cur_set2.end(); set_iter++ ){
			match_file << '\t' << (MatchID_t)*set_iter;
		}
		match_file << endl;
	}
}

void MatchList::MultiplicityFilter( unsigned mult ){

	vector< Match* > nway_list;
	for( uint memI = 0; memI < size(); memI++ ){
		if( (*this)[ memI ]->Multiplicity() == mult )
			nway_list.push_back( (*this)[ memI ] );
		else{
			(*this)[ memI ]->UnlinkSelf();
			(*this)[ memI ]->Free();
			(*this)[ memI ] = NULL;
		}
	}
	vector< Match* >::operator=( nway_list );
}

void MatchList::LengthFilter( gnSeqI length ){

	size_t cur = 0;
	for( size_t memI = 0; memI < size(); memI++ ){
		if( (*this)[ memI ]->Length() >= length )
			(*this)[cur++] = (*this)[memI];
		else{
			(*this)[ memI ]->UnlinkSelf();
			(*this)[ memI ]->Free();
			(*this)[ memI ] = NULL;
		}
	}
	this->resize(cur);
}

/*
void MatchList::ExactFilter( valarray< bool >& filter_spec ){
	ErrorMsg( "MatchList::ExactFilter() needs to be re-implemented\n" );
	vector<Match*>::iterator match_iter;
	vector<Match*>::iterator to_del;
	match_iter = begin();
	while( match_iter != end() ){
		uint64 matchnumber = (*match_iter)->MatchNumber();
		valarray< bool > matchnum( false, (*match_iter)->SeqCount());
		for( uint32 seqI = (*match_iter)->SeqCount(); seqI > 0; seqI-- ){
			if( matchnumber & 0x1 )
				matchnum[ seqI - 1 ] = true;
			matchnumber >>= 1;
		}

		valarray< bool > equal = matchnum == filter_spec;
		unsigned msum = equal.min();
		if( msum == false ){
			// delete this one.
			to_del = match_iter;
			match_iter++;
			erase( to_del );
		}else
			match_iter++;
	}
}

void MatchList::IntersectFilter( valarray< bool >& filter_spec ){
	ErrorMsg( "MatchList::IntersectFilter() needs to be re-implemented\n" );
	vector<Match*>::iterator match_iter;
	vector<Match*>::iterator to_del;
	match_iter = begin();
	while( match_iter != end() ){
		uint64 matchnumber = (*match_iter)->MatchNumber();
		valarray< bool > matchnum( false, (*match_iter)->SeqCount());
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
			erase( to_del );
		}else
			match_iter++;
	}
}
*/
void MatchList::UnlinkedFirstStartFilter( unsigned start_seq )
{
	ErrorMsg( "MatchList::UnlinkedFirstStartFilter() needs to be re-implemented\n" );
	vector<Match*>::iterator match_iter;
	vector<Match*>::iterator to_del;
	match_iter = begin();
	while( match_iter != end() ){
		if( (*match_iter)->FirstStart() != start_seq ||
			(*match_iter)->Supersets().size() > 0 ) 
		{
			to_del = match_iter;
			match_iter++;
			erase( to_del );
			continue;
		}
		match_iter++;
	}
}

} // namespace mems

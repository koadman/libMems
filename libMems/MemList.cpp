#include "MemList.h"
#include "MemHash.h"

MemList::MemList( const MemList& ml ){
	*this = ml;
}

MemList& MemList::operator=( const MemList& ml ){
	mem_list = ml.mem_list;
	sml_filename = ml.sml_filename;
	seq_filename = ml.seq_filename;
	sml_table = ml.sml_table;
	seq_table = ml.seq_table;
	return *this;
}

void MemList::ReadList(istream& mem_file){
	string tag;
	gnSeqI len;
	int64 start;
	MemHashEntry mhe;
	unsigned int seq_count;
	unsigned int mer_size;
	
	mem_file >> tag;	//format version tag
	if( tag != "FormatVersion" ){
		Throw_gnEx(InvalidFileFormat());
	}
	mem_file >> tag;	//format version
	if( tag != "1" ){
		Throw_gnEx(InvalidFileFormat());
	}
	mem_file >> tag;	//sequence count tag
	if( tag != "SequenceCount" ){
		Throw_gnEx(InvalidFileFormat());
	}
	mem_file >> seq_count;	//sequence count
	if(seq_count < 2){
		Throw_gnEx(InvalidFileFormat());
	}
	mem_file >> tag;	//mer size tag
	if( tag != "MerSize" ){
		Throw_gnEx(InvalidFileFormat());
	}
	mem_file >> mer_size;	//mer size
	if(seq_count < 2){
		Throw_gnEx(InvalidFileFormat());
	}

	mhe = MemHashEntry(seq_count, mer_size, MemHashEntry::extended);
	
	// read the sequence file names and lengths
	for( unsigned int seqI = 0; seqI < seq_count; seqI++ ){
		mem_file >> tag;	// name tag
		mem_file >> tag;	// name
		seq_filename.push_back(tag);
//		try{
//			gnSequence *new_seq = new gnSequence();
//			new_seq->LoadSource(tag);
//			seq_table.push_back( new_seq );
//		}catch( gnException& gne );
		mem_file >> tag;	// length tag
		mem_file >> tag;	// length
	}
	// read the number of matches
	unsigned int match_count;
	mem_file >> tag;	// match count tag
	mem_file >> match_count;	// match count
	
	while(true){
		mem_file >> len;
		if(!mem_file.good() || mem_file.eof() )
			break;
		mhe.SetLength(len);
		for(uint32 seqI = 0; seqI < seq_count; seqI++){
			mem_file >> start;
			mhe.SetStart(seqI, start);
		}
		//break if the stream ended
		if(!mem_file.good() || mem_file.eof() )
			break;
		mhe.CalculateOffset();
		mem_list.push_back( new MemHashEntry(mhe) );
	}
	if( match_count != mem_list.size() ){
		Throw_gnEx(InvalidFileFormat());
	}
}

void MemList::WriteList(ostream& mem_file) const{
	if( mem_list.size() == 0 )
		return;
	MemHashEntry* first_mem = *(mem_list.begin());
	unsigned int seq_count = first_mem->SeqCount();
	unsigned int mer_size = first_mem->MerSize();

	mem_file << "FormatVersion" << '\t' << 1 << "\n";
	mem_file << "SequenceCount" << '\t' << seq_count << "\n";
	mem_file << "MerSize" << '\t' << mer_size << "\n";
	for(unsigned int seqI = 0; seqI < seq_count; seqI++){
		mem_file << "Sequence" << seqI << "File" << '\t';
		if( seq_filename.size() > seqI )
			mem_file << seq_filename[seqI];
		else
			mem_file << "null";
		mem_file << "\n";
		mem_file << "Sequence" << seqI << "Length" << '\t';
		if( seq_table.size() > seqI )
			mem_file << seq_table[seqI]->length();
		else
			mem_file << "0";
		mem_file << "\n";
	}

	mem_file << "MatchCount" << '\t' << mem_list.size() << endl;

	//get all the mems out of the hash table and write them out
    list<MemHashEntry*>::const_iterator mem_iter;
	mem_iter = mem_list.begin();
	for(; mem_iter != mem_list.end(); mem_iter++)
		mem_file << **mem_iter << "\n";
}

void MemList::MultiplicityFilter( unsigned mult ){
	list<MemHashEntry*>::iterator mem_iter;
	list<MemHashEntry*>::iterator to_del;
	mem_iter = mem_list.begin();
	while( mem_iter != mem_list.end() ){
		if( (*mem_iter)->Multiplicity() != mult ){
			// delete this one.
			to_del = mem_iter;
			mem_iter++;
			mem_list.erase( to_del );
		}else
			mem_iter++;
	}
}

void MemList::ExactFilter( valarray<bool>& filter_spec ){
	list<MemHashEntry*>::iterator mem_iter;
	list<MemHashEntry*>::iterator to_del;
	mem_iter = mem_list.begin();
	while( mem_iter != mem_list.end() ){
		uint64 matchnumber = (*mem_iter)->MatchNumber();
		valarray<bool> matchnum( false, (*mem_iter)->SeqCount());
		for( uint32 seqI = (*mem_iter)->SeqCount(); seqI > 0; seqI-- ){
			if( matchnumber & 0x1 )
				matchnum[ seqI - 1 ] = true;
			matchnumber >>= 1;
		}

		valarray<bool> equal = matchnum == filter_spec;
		unsigned msum = equal.min();
		if( msum == false ){
			// delete this one.
			to_del = mem_iter;
			mem_iter++;
			mem_list.erase( to_del );
		}else
			mem_iter++;
	}
}

void MemList::IntersectFilter( valarray<bool>& filter_spec ){
	list<MemHashEntry*>::iterator mem_iter;
	list<MemHashEntry*>::iterator to_del;
	mem_iter = mem_list.begin();
	while( mem_iter != mem_list.end() ){
		uint64 matchnumber = (*mem_iter)->MatchNumber();
		valarray<bool> matchnum( false, (*mem_iter)->SeqCount());
		for( uint32 seqI = (*mem_iter)->SeqCount(); seqI > 0; seqI-- ){
			if( matchnumber & 0x1 )
				matchnum[ seqI - 1 ] = true;
			matchnumber >>= 1;
		}

		matchnum &= filter_spec;
		unsigned msum = matchnum.max();
		if( msum == 0 ){
			// delete this one.
			to_del = mem_iter;
			mem_iter++;
			mem_list.erase( to_del );
		}else
			mem_iter++;
	}
}

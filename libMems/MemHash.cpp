#include "MemHash.h"
#include "gn/gnFilter.h"
#include <list>
#include <map>
#include <sstream>
#include "MemSubsets.h"

MemHash::MemHash() : MatchFinder(){
	table_size = DEFAULT_MEM_TABLE_SIZE;
	mer_size = 0;
	seq_count = 0;
	m_mem_count = 0;
	m_collision_count = 0;
	m_repeat_tolerance = DEFAULT_REPEAT_TOLERANCE;
	m_enumeration_tolerance = DEFAULT_ENUMERATION_TOLERANCE;
	//allocate the hash table
	mem_table = new set<MemHashEntry*, MheCompare>[table_size];
	mem_table_count.reserve( table_size );
	for(uint32 i=0; i < table_size; i++)
		mem_table_count.push_back(0);
}

//make sure this calls the destructor on each element
MemHash::~MemHash(){
	if(mem_table != NULL){
		delete[] mem_table;
		mem_table = NULL;
	}
}

MemHash::MemHash(const MemHash& mh) : MatchFinder(mh){
	table_size = mh.table_size;
	mer_size = mh.mer_size;
	seq_count = mh.seq_count;
	m_mem_count = mh.m_mem_count;
	m_collision_count = mh.m_collision_count;
	m_repeat_tolerance = mh.m_repeat_tolerance;
	m_enumeration_tolerance = mh.m_enumeration_tolerance;
	//allocate the hash table
	mem_table = new set<MemHashEntry*, MheCompare>[table_size];
	for(uint32 i=0; i < table_size; i++){
		mem_table_count.push_back(mh.mem_table_count[i]);
		mem_table[i] = mh.mem_table[i];
	}
}

MemHash* MemHash::Clone() const{
	return new MemHash(*this);
}

void MemHash::Clear() {
	m_mem_count = 0;
	m_collision_count = 0;
	m_repeat_tolerance = DEFAULT_REPEAT_TOLERANCE;
	m_enumeration_tolerance = DEFAULT_ENUMERATION_TOLERANCE;
	//clear the hash table
	for(uint32 listI = 0; listI < table_size; listI++){
		mem_table[listI].clear();
		mem_table_count.push_back(0);
	}
	mem_table_count.clear();
	for(uint32 i=0; i < table_size; i++)
		mem_table_count.push_back(0);
}

void MemHash::SetTableSize(uint32 new_table_size){
	//deallocate the hash table
	if(mem_table != NULL)
		delete[] mem_table;
	//allocate the hash table
	table_size = new_table_size;
	mem_table = new set<MemHashEntry*, MheCompare>[table_size];
	mem_table_count.clear();
	for(uint32 i=0; i < table_size; i++)
		mem_table_count.push_back(0);
}

boolean MemHash::CreateMatches(){
	MatchFinder::FindMatches();
	return true;
}

void MemHash::FindMatches( MatchList& ml ) {
	for( uint32 seqI = 0; seqI < ml.seq_table.size(); seqI++ )
		AddSequence( ml.sml_table[ seqI ], ml.seq_table[ seqI ] );
	MatchFinder::FindMatches();
	
	vector<string> seq_filename = ml.seq_filename;
	vector<string> sml_filename = ml.sml_filename;

	ml = GetMatchList();
	ml.seq_filename = seq_filename;
	ml.sml_filename = sml_filename;
	
}

void MemHash::FindMems( MemList& ml ) {
	for( uint32 seqI = 0; seqI < ml.seq_table.size(); seqI++ )
		AddSequence( ml.sml_table[ seqI ], ml.seq_table[ seqI ] );
	MatchFinder::FindMatches();

	for(uint32 i=0; i < table_size; i++)
		ml.mem_list.insert(ml.mem_list.end(), mem_table[i].begin(), mem_table[i].end() );
}


//this code is all pretty useless for dealing with repeats
void MemHash::BreakOnAmbiguities(){
}

MatchList MemHash::GetMatchList() const{
	MatchList ml;
	ml.seq_table = seq_table;
	ml.sml_table = sar_table;
	
	list<MemHashEntry*> mem_list;
	for(uint32 i=0; i < table_size; i++)
		mem_list.insert(mem_list.end(), mem_table[i].begin(), mem_table[i].end() );
	
	
/*	map< MemHashEntry*, Match* > match_map;
	list<MemHashEntry*>::iterator mem_iter = mem_list.begin();
	for( ; mem_iter != mem_list.end(); mem_iter++ ){
		Match* new_match = new Match( **mem_iter );
		match_map.insert( map< MemHashEntry*, Match* >::value_type( *mem_iter, new_match ) );
		ml.match_list.push_back( new_match );
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
*/
	ml.FromMemList( mem_list );
	return ml;
}

void MemHash::GetMemList( vector<MemHashEntry*>& mem_list ) const{
	mem_list.clear();
	for(uint32 i=0; i < table_size; i++)
		mem_list.insert(mem_list.end(), mem_table[i].begin(), mem_table[i].end() );
}

MemList MemHash::GetMemList() const{
	MemList ml;
	ml.seq_table = seq_table;
	ml.sml_table = sar_table;
	for(uint32 i=0; i < table_size; i++)
		ml.mem_list.insert(ml.mem_list.end(), mem_table[i].begin(), mem_table[i].end() );

	return ml;
}

boolean MemHash::EnumerateMatches( list<idmer>& match_list ){

	match_list.sort(&idmer_id_lessthan);
	list<idmer>::iterator iter = match_list.begin();
	list<idmer>::iterator iter2 = match_list.begin();
	list<idmer> hash_list;
	iter2++;
	int32 cur_tolerance = m_enumeration_tolerance - 1;
	uint32 cur_count = 0;	//the number of mers from a particular sequence
	hash_list.push_back(*iter);
	for(; iter2 != match_list.end(); iter++){
		if(iter->id != iter2->id){
			hash_list.push_back(*iter2);
			cur_tolerance = m_enumeration_tolerance - 1;
			cur_count = 0;
		}else{
			if(cur_tolerance > 0){
				hash_list.push_back(*iter2);
				cur_tolerance--;
			}
			cur_count++;
		}
		iter2++;
		if(cur_count > m_repeat_tolerance)
			return true;
	}
	if(hash_list.size() > 1){
		if(m_enumeration_tolerance == 1)
			return HashMatch(hash_list);
		else
			return MatchFinder::EnumerateMatches( hash_list );
	}
	return true;
}


//why have separate hash tables?
// MemHashEntries use GENETICIST coordinates.  They start at 1, not 0.
boolean MemHash::HashMatch(list<idmer>& match_list){
	//check that there is at least one forward component
	match_list.sort(&idmer_id_lessthan);
	// initialize the hash entry
	MemHashEntry mhe = MemHashEntry(seq_count, mer_size);
	mhe.SetLength(mer_size);
	
	//Fill in the new MemHashEntry and set direction parity if needed.
	list<idmer>::iterator iter = match_list.begin();
	for(; iter != match_list.end(); iter++)
		mhe.SetStart(iter->id, iter->position + 1);
	SetDirection(mhe);
	mhe.CalculateOffset();
	if(mhe.Multiplicity() < 2){
		cout << "red fag " << mhe << "\n";
	}else 
		AddHashEntry(mhe);

	return true;
}

void MemHash::SetDirection(MemHashEntry& mhe){
	//get the reference direction
	boolean ref_forward;
	uint32 seqI=0;
	for(; seqI < mhe.SeqCount(); seqI++)
		if(mhe[seqI] != MEM_NO_MATCH){
			ref_forward = !(GetSar(seqI)->GetMer(mhe[seqI] - 1) & 0x1);
			break;
		}
	//set directional parity for the rest
	for(seqI++; seqI < mhe.SeqCount(); seqI++)
		if(mhe[seqI] != MEM_NO_MATCH)
			if(ref_forward == (GetSar(seqI)->GetMer(mhe[seqI] - 1) & 0x1))
				mhe.SetStart(seqI, -mhe[seqI]);
}

// Tries to add a new mem to the mem hash table
// If the mem already exists in the table, a pointer to it
// is returned.  Otherwise mhe is added and a pointer to
// it is returned.
MemHashEntry* MemHash::AddHashEntry(MemHashEntry& mhe){
	//first compute which hash table bucket this is going into
	int64 offset = mhe.Offset();

	uint32 bucketI = ((offset % table_size) + table_size) % table_size;
    set<MemHashEntry*, MheCompare>::iterator insert_he;
	insert_he = mem_table[bucketI].find(&mhe);
	if( insert_he != mem_table[bucketI].end()){
		m_collision_count++;
		return *insert_he;
	}
	
	//if we made it this far there were no collisions
	//extend the mem into the surrounding region.
	vector<MemHashEntry> subset_matches;
	if( !mhe.Extended() )
		ExtendMatch(mhe, subset_matches);
	MemHashEntry *new_mhe = new MemHashEntry(mhe);

	// can't insert until after the extend!!
	mem_table[bucketI].insert(new_mhe);
	
	// link up the subset matches
	for(uint32 subsetI = 0; subsetI < subset_matches.size(); subsetI++){
		MemHashEntry* submem = AddHashEntry( subset_matches[ subsetI ] );
		new_mhe->LinkSubset( submem );
	}
	
	mem_table_count[bucketI]++;
	m_mem_count++;
	return new_mhe;
}

void MemHash::PrintDistribution(ostream& os) const{
    set<MemHashEntry*, MheCompare>::const_iterator mem_iter;
	gnSeqI base_count;
	for(uint32 i=0; i < mem_table_count.size(); i++){
		mem_iter = mem_table[i].begin();
		base_count = 0;
		for(; mem_iter != mem_table[i].end(); mem_iter++){
			base_count += (*mem_iter)->Length();
		}
		os << i << '\t' << mem_table_count[i] << '\t' << base_count << '\n';
	}
}

void MemHash::LoadFile(istream& mem_file){
	string tag;
	gnSeqI len;
	int64 start;
	MemHashEntry mhe;
	
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
	mhe = MemHashEntry(seq_count, mer_size, MemHashEntry::extended);
	
	// read the sequence file names and lengths
	for( unsigned int seqI = 0; seqI < seq_count; seqI++ ){
		mem_file >> tag;	// name tag
		mem_file >> tag;	// name
		mem_file >> tag;	// length tag
		mem_file >> tag;	// length
	}
	// read the number of matches
	unsigned int match_count;
	mem_file >> tag;	// match count tag
	mem_file >> match_count;	// match count
	
	while(mem_file.good()){
		mem_file >> len;
		if(!mem_file.good())
			break;
		mhe.SetLength(len);
		for(uint32 seqI = 0; seqI < seq_count; seqI++){
			mem_file >> start;
			mhe.SetStart(seqI, start);
		}
		//break if the stream ended
		if(!mem_file.good())
			break;
		mhe.CalculateOffset();
		AddHashEntry(mhe);
	}
	if( match_count != m_mem_count ){
		stringstream mystr;
		mystr << "Expected " << match_count << " mems, but read " << m_mem_count << " mems.";
		Throw_gnExMsg(InvalidFileFormat(), mystr.str().c_str() );
	}
}

void MemHash::WriteFile(ostream& mem_file) const{
	mem_file << "FormatVersion" << '\t' << 1 << "\n";
	mem_file << "SequenceCount" << '\t' << sar_table.size() << "\n";
	for(unsigned int seqI = 0; seqI < seq_count; seqI++){
		mem_file << "Sequence" << seqI << "File";
		gnGenomeSpec* specker = seq_table[seqI]->GetSpec();
		string sourcename = specker->GetName();
		if( sourcename == "" )
			sourcename = "null";
		mem_file << '\t' << sourcename << "\n";
		mem_file << "Sequence" << seqI << "Length";
		mem_file << '\t' << seq_table[seqI]->length() << "\n";
	}
	mem_file << "MatchCount" << '\t' << m_mem_count << endl;
	//get all the mems out of the hash table and write them out
    set<MemHashEntry*, MheCompare>::const_iterator mem_table_iter;
	for(uint32 i=0; i < table_size; i++){
		mem_table_iter = mem_table[i].begin();
		for(; mem_table_iter != mem_table[i].end(); mem_table_iter++)
			mem_file << **mem_table_iter << "\n";
	}
}

#include "RepeatHash.h"
#include <list>

RepeatHash::RepeatHash(){
}

RepeatHash::~RepeatHash(){
}

RepeatHash::RepeatHash(const RepeatHash& mh) : MemHash(mh){

}

RepeatHash* RepeatHash::Clone() const{
	return new RepeatHash(*this);
}

boolean RepeatHash::CreateMatches(){
	if(seq_count == 1)
		return MatchFinder::FindMatches();

	return false;
}

boolean RepeatHash::EnumerateMatches( list<idmer>& match_list ){
	return HashMatch(match_list);
}

//why have separate hash tables?
// RepeatHashEntries use GENETICIST coordinates.  They start at 1, not 0.
boolean RepeatHash::HashMatch(list<idmer>& match_list){
	//check that there is at least one forward component
	match_list.sort(&idmer_id_lessthan);
	// initialize the hash entry
	MemHashEntry mhe = MemHashEntry(seq_count, mer_size);
	mhe.SetLength(mer_size);
	
	//Fill in the new MemHashEntry and set direction parity if needed.
	list<idmer>::iterator iter = match_list.begin();

	uint32 repeatI = 0;
	for(; iter != match_list.end(); iter++)
		mhe.SetStart(repeatI++, iter->position + 1);

	SetDirection( mhe );
	mhe.CalculateOffset();
	if(mhe.Multiplicity() < 2){
		cout << "red fag " << mhe << "\n";
	}else{
		AddHashEntry(mhe);
	}
	return true;
}


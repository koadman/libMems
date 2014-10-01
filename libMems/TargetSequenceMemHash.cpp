/*******************************************************************************
 * $Id: TargetSequenceMemHash.cpp,v 1.3 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/TargetSequenceMemHash.h"
#include <list>

using namespace std;
using namespace genome;
namespace mems {

TargetSequenceMemHash::TargetSequenceMemHash(){
	target = 0;
}


TargetSequenceMemHash::TargetSequenceMemHash(const TargetSequenceMemHash& mh) : MemHash(mh){
	*this = mh;
}

TargetSequenceMemHash& TargetSequenceMemHash::operator=( const TargetSequenceMemHash& mh ){
	target = mh.target;
	return *this;
}

TargetSequenceMemHash* TargetSequenceMemHash::Clone() const{
	return new TargetSequenceMemHash(*this);
}

boolean TargetSequenceMemHash::HashMatch(list<idmer>& match_list){
	//check that there is at least one forward component
	match_list.sort(&idmer_id_lessthan);
	// initialize the hash entry
	MatchHashEntry mhe = MatchHashEntry(seq_count, GetSar(0)->SeedLength());
	mhe.SetLength(GetSar(0)->SeedLength());
	
	//Fill in the new Match and set direction parity if needed.
	list<idmer>::iterator iter = match_list.begin();
	for(; iter != match_list.end(); iter++)
		mhe.SetStart(iter->id, iter->position + 1);
	SetDirection(mhe);
	mhe.CalculateOffset();
	uint64 match_number = 0;

	// Add it if the target sequence is present
	if(mhe.Start(target) != NO_MATCH)
		AddHashEntry(mhe);

	return true;
}

void TargetSequenceMemHash::FindSubsets(const Match& mhe, std::vector<Match>& subset_matches){
	// quick hack: just filter down all subsets to just those which include the target sequence
	// this could be reimplemented for added speed by checking only for subsets with the target
	std::vector<Match> temp_matches;
	MatchFinder::FindSubsets(mhe, temp_matches);
	for(int i = 0; i < temp_matches.size(); i++){
		if(temp_matches[i].Start(target) != NO_MATCH){
			subset_matches.push_back(temp_matches[i]);
		}
	}
}

} // namespace mems

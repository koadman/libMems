#ifndef _RepeatHashCat_h_
#define _RepeatHashCat_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MemHash.h"
#include "libMems/RepeatHash.h"
#include "libMems/RepeatMatchList.h"
#include "libMems/RepeatMatch.h"

namespace mems {

class RepeatHashCat : public RepeatHash
{
public:
	//punt tjt: needed to add this to track where concatenated sequence starts
	std::vector<uint32> concat_contig_start; // number of contigs in each sequence

	RepeatHashCat( void );
	~RepeatHashCat( void );
	RepeatHashCat(const RepeatHashCat& mh);
	std::vector<uint32> concatContigStart(void);
	//vector<uint32> setConcatContigStart(void);
	//virtual RepeatHash* Clone() const;
	virtual boolean CreateMatches();
	virtual void FindMatchesFromPosition( RepeatMatchList& match_list, const std::vector<gnSeqI>& start_points );
	
protected:
	virtual boolean EnumerateMatches( IdmerList& match_list );
	virtual SortedMerList* GetSar(uint32 sarI) const;
	std::set<RepeatMatch*, MheCompare>* mem_table;


};

}	// namespace mems

#endif



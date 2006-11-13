#ifndef _RepeatHash_h_
#define _RepeatHash_h_

#include "MemHash.h"

/**
 * Finds repeats within a single sequence.
 * This class extends the functionality of MemHash to search for repetitive
 * matches within a single sequence.
 */
class RepeatHash : public MemHash{
public:
	RepeatHash();
	~RepeatHash();

	RepeatHash(const RepeatHash& mh);
	virtual RepeatHash* Clone() const;
	virtual boolean CreateMatches();
protected:

	virtual boolean EnumerateMatches( list<idmer>& match_list );
	virtual boolean HashMatch(list<idmer>& match_list);
	virtual SortedMerList* GetSar(uint32 sarI) const;
};


inline
SortedMerList* RepeatHash::GetSar(uint32 sarI) const{
	return sar_table[0];
}

inline
bool idmer_greaterthan(idmer& a_v, idmer& m_v){
	return (a_v.mer < m_v.mer);// ? true : false;
};


#endif //_RepeatHash_h_

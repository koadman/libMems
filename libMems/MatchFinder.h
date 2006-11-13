#ifndef _MatchFinder_h_
#define _MatchFinder_h_

#include "SuffixArray.h"
#include "MemHashEntry.h"
#include <list>

struct idmer{
	gnSeqI	position;	//starting position of this mer in the genome
	sarID_t	id;			//the sequence identifier.
	uint64 	mer; 		//the actual sequence
};

class MatchFinder : public gnClone{
public:
	MatchFinder();
	~MatchFinder();
	MatchFinder(MatchFinder& mf);
	virtual boolean AddSuffixArray( SuffixArray* sar, gnSequence* seq = NULL );
	virtual uint32 SequenceCount(void){return seq_count;};
	virtual void SetAmbiguityTolerance(uint32 ambiguity_tol){ambiguity_tolerance = ambiguity_tol;}
	virtual uint32 AmbiguityTolerance(){return ambiguity_tolerance;}
protected:
	virtual boolean FindMatches(uint32 mer_mask_size, gnSeqI startI = 0, gnSeqI lengthI = GNSEQI_END);
	virtual boolean HashMatch(list<idmer>& match_list) = 0;
	virtual boolean EnumerateMatches(list<idmer>& match_list);
	virtual void ExtendMatch(MemHashEntry& mhe, gnSeqI max_backward = GNSEQI_END, gnSeqI max_forward = GNSEQI_END);
	virtual boolean MatchAmbiguities(MemHashEntry& mhe, uint32 match_size);
	virtual boolean WarpExtend(MemHashEntry& mhe, gnSeqI max_backward = GNSEQI_END, gnSeqI max_forward = GNSEQI_END);
	
	vector<SuffixArray*> sar_table;
	vector<gnSequence*> seq_table;
	vector<int8> seqid_table;
	
	uint32 mer_size;
	uint32 seq_count;
	uint32 sarid_table[UINT8_MAX];
	uint32 ambiguity_tolerance;
};

inline
bool idmer_lessthan(idmer& a_v, idmer& m_v){
	return (a_v.mer < m_v.mer);// ? true : false;
};

//id less than function for STL sort functions
inline
bool idmer_id_lessthan(idmer& a_v, idmer& m_v){
	return (a_v.id < m_v.id);// ? true : false;
};


#endif
#ifndef _MatchFinder_h_
#define _MatchFinder_h_

#include "SortedMerList.h"
#include "MemHashEntry.h"
#include "MemList.h"
#include <list>
#include <iostream>

struct idmer{
	gnSeqI	position;	//starting position of this mer in the genome
	sarID_t	id;			//the sequence identifier.
	uint64 	mer; 		//the actual sequence
};

const unsigned int PROGRESS_GRANULARITY = 100;

/**
 * This pure virtual class implements a general framework for finding
 * exactly matching mers.  It is extended by the MemHash and MemScorer
 * classes.
 * @see MemHash
 * @see MemScorer
 */
class MatchFinder : public gnClone{
public:
	MatchFinder();
	~MatchFinder();
	MatchFinder(const MatchFinder& mf);
	virtual void Clear();
	/**
	 * Adds a sequence to use when searching for exact matches.
	 * @param sar A pointer to the sorted mer list for the new sequence
	 * @param seq A pointer to the gnSequence corresponding to the new sequence.
	 */
	virtual boolean AddSequence( SortedMerList* sar, gnSequence* seq = NULL );
	/**
	 * Given the index of a sequence and an index into the sorted mer list, this function
	 * will search the other sorted mer lists for the same mer.  This function returns the
	 * position of the mer in each sequence in the breakpoints vector.
	 */
	virtual void GetBreakpoint( uint32 sarI, gnSeqI startI, vector<gnSeqI>& breakpoints ) const;
	virtual uint32 Multiplicity(void){return seq_count;};
	/** NOT IMPLEMENTED: Sets the number of ambiguities allowed in a mer match*/
	virtual void SetAmbiguityTolerance(uint32 ambiguity_tol){ambiguity_tolerance = ambiguity_tol;}
	/** @return the number of ambiguities allowed in a mer match */
	virtual uint32 AmbiguityTolerance(){return ambiguity_tolerance;}
	/** @return The progress of the current operation.  Ranges from 0 to 100.  -1 indicates no computation is being performed */
	virtual float GetProgress() const {return m_progress;}

	/**
	 * Logs progress to the designated ostream.  Set to null to skip progress logging.
	 */
	virtual void LogProgress( ostream* os );
protected:
	/** Finds all the matches between the sequences */
	virtual boolean FindMatches();
	/** 
	 * Finds all matches in a designated range of the sequence's sorted mer lists 
	 * @throws InvalidData thrown if the start_points are bad or if the sorted mer lists were sorted on different mer sizes
	 * @return true if completed searching, false if repetitive mers were encountered and FindMatches must be called again.
	 */
	virtual boolean FindMatches(vector<gnSeqI>& start_points, vector<gnSeqI>& search_len);
	/** Called whenever a mer match is found */
	virtual boolean HashMatch(list<idmer>& match_list) = 0;
	virtual boolean EnumerateMatches(list<idmer>& match_list);
	virtual void FindSubsets(const MemHashEntry& mhe, vector<MemHashEntry>& subset_matches);
	virtual void ExtendMatch(MemHashEntry& mhe, vector<MemHashEntry>& subset_matches, gnSeqI max_backward = GNSEQI_END, gnSeqI max_forward = GNSEQI_END);
	virtual boolean MatchAmbiguities(MemHashEntry& mhe, uint32 match_size);

	virtual SortedMerList* GetSar(uint32 sarI) const;
	vector<SortedMerList*> sar_table;
	vector<gnSequence*> seq_table;
	vector<int8> seqid_table;
	
	uint32 mer_size;
	uint32 seq_count;
	uint32 sarid_table[UINT8_MAX];
	uint32 ambiguity_tolerance;
	
	// for subset matches
	vector< vector< uint32 > > alpha_map;
	uint alpha_map_size;
	uint alphabet_bits;
	
	float m_progress;
	ostream* log_stream;
};

/** 
 * InvalidData exceptions are thrown when the input to an algorithm is invalid
 */
CREATE_EXCEPTION( InvalidData );

inline
SortedMerList* MatchFinder::GetSar(uint32 sarI) const{
	return sar_table[sarI];
}

inline
bool idmer_lessthan(idmer& a_v, idmer& m_v){
	return (a_v.mer < m_v.mer);// ? true : false;
};

//id less than function for STL sort functions
inline
bool idmer_id_lessthan(idmer& a_v, idmer& m_v){
	return (a_v.id < m_v.id);// ? true : false;
};


#endif	//_MatchFinder_h_

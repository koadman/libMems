#ifndef _MemScorer_h_
#define _MemScorer_h_

#include "MatchFinder.h"

#define MATCH_NOT_SCORED 2
#define MATCH_DEFAULT_HIT_RESOLUTION 10

class MemScorer : public MatchFinder{
public:
	MemScorer(uint32 hit_resolution = MATCH_DEFAULT_HIT_RESOLUTION);
	~MemScorer();
	virtual MemScorer* Clone() const;
	
	//gives a score between 0 and 1 which is the ratio of matching
	//mers to total mers
	float64 PercentHitsScore();
	
	//Score which should roughly increase when mem size increases
	//Maybe accomplish by extending random matching pieces and collecting
	//their average size.  The likelihood that two independent samples belong
	//to the same mem is mem_size/seq_size.  How do I reward lots of
	//small mems which are close together?
	float64 SubsequenceIdentityScore();
	
	//Heuristic Hits Score
	//A faster version of PercentHitsScore which samples a portion
	//of the suffix arrays for matches
	float64 HeuristicHitsScore();
protected:
	virtual boolean EnumerateMatches(list<idmer>& match_list);
	virtual boolean HashMatch(list<idmer>& match_list);
	void AllocateMatrix();

	float64 percent_hits;
	float64 subseq_identity;
	float64 heuristic_hits;
	
	//3d array seq_count x seq_count x hit_resolution
	uint32 m_hit_resolution;
	uint32 **hit_matrix;
	uint32 **identity_matrix;
};

#endif // _MemScorer_h_

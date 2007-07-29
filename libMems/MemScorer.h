/*******************************************************************************
 * $Id: MemScorer.h,v 1.6 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _MemScorer_h_
#define _MemScorer_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MemHash.h"

namespace mems {

#define MATCH_NOT_SCORED 2
#define MATCH_DEFAULT_HIT_RESOLUTION 10
#define MAX_DETAIL_SEQS 28	///defines the maximum number of sequences that
							///we can do a detailed count of their matchings.
							///the amount of space required is 4 bytes * 2 ^ MAX_DETAIL_SEQS
							///so for MAX_DETAIL_SEQS = 28 we need about a Gigabyte of ram.

class MemScorer : public MemHash
{
public:
	MemScorer(uint32 hit_resolution = MATCH_DEFAULT_HIT_RESOLUTION);
	MemScorer( const MemScorer& ms );
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
	//of the sorted mer lists for matches
	float64 HeuristicHitsScore();

	void PrintIDMatrix(std::ostream& os);
	void PrintHitMatrix(std::ostream& os);
	void PrintDetailList(std::ostream& os);
protected:
	virtual boolean HashMatch(IdmerList& match_list);
	void AllocateMatrix();

	float64 percent_hits;
	float64 subseq_identity;
	float64 heuristic_hits;
	
	//3d array seq_count x seq_count x hit_resolution
	uint32 m_hit_resolution;
	gnSeqI **hit_matrix;
	gnSeqI **identity_matrix;
	gnSeqI *detail_list;
};

}

#endif // _MemScorer_h_

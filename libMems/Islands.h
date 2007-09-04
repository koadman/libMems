/*******************************************************************************
 * $Id: Islands.h,v 1.7 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __Islands_h__
#define __Islands_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include "libMems/SubstitutionMatrix.h"
#include "libMems/IntervalList.h"
#include "libMems/NumericMatrix.h"
#include "libMems/MatchList.h"
#include "libMems/GappedAlignment.h"
#include "libMems/CompactGappedAlignment.h"
#include "libMems/Aligner.h"
#include <boost/multi_array.hpp>
#include "libMems/HomologyHMM/homology.h"
#include "libMems/Scoring.h"

namespace mems {

/**
 * A class to represent an island in an alignment.  Islands are generally
 * large insertions of a region of sequence relative to
 * another sequence.
 */
class Island{
public:
	uint seqI;
	uint seqJ;
	int64 leftI;
	int64 leftJ;
	int64 rightI;
	int64 rightJ;
};

/**
 * Identifies gaps in the alignment between pairs of sequences that are longer than
 * some number of base pairs in length.  Prints islands to an output stream
 */
void simpleFindIslands( IntervalList& iv_list, uint island_size, std::ostream& island_out );
void findIslandsBetweenLCBs( IntervalList& iv_list, uint island_size, std::ostream& island_out );
void simpleFindIslands( IntervalList& iv_list, uint island_size, std::vector< Island >& island_list );

class HssCols{
public:
	uint seqI;
	uint seqJ;
	size_t left_col;
	size_t right_col;
};

typedef std::vector< HssCols > hss_list_t;
typedef boost::multi_array< hss_list_t, 3 > hss_array_t;

typedef HssCols IslandCols;	// use the same structure for island segs

void findHssRandomWalkScoreVector( std::vector< score_t > scores, score_t significance_threshold, hss_list_t& hss_list, uint seqI = 0, uint seqJ = 0, boolean left_homologous = false, boolean right_homologous = false );

template<typename MatchVector>
void findHssRandomWalk( const MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table, const PairwiseScoringScheme& scoring, score_t significance_threshold, hss_array_t& hss_array, boolean left_homologous = false, boolean right_homologous = false );

template<typename MatchVector>
void hssColsToIslandCols( const MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table, std::vector< HssCols >& hss_list, std::vector< IslandCols >& island_col_list );

template<typename MatchVector>
void findHssRandomWalkCga( const MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table, const PairwiseScoringScheme& scoring, score_t significance_threshold, std::vector< CompactGappedAlignment<>* >& hss_list );

template<typename MatchVector>
void findIslandsRandomWalkCga( const MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table, const PairwiseScoringScheme& scoring, score_t significance_threshold, std::vector< CompactGappedAlignment<>* >& island_list );

template<typename MatchVector>
void findIslandsRandomWalk( const MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table, const PairwiseScoringScheme& scoring, score_t significance_threshold, std::vector< Island >& island_list );

/**
 *  Find regions in each sequence that do not belong to any LCB, add them to their own
 * Interval (LCB) in the IntervalList.
 */
void addUnalignedIntervals( IntervalList& iv_list, std::set< uint > seq_set = std::set< uint >(), std::vector<gnSeqI> seq_lengths = std::vector<gnSeqI>() );

/**
 * Identifies stretches of alignment existing in all sequences that doesn't
 * contain a gap larger than a particular size.  Such regions are considered
 * the backbone of the alignment.
 */
void simpleFindBackbone( IntervalList& iv_list, uint backbone_size, uint max_gap_size, std::vector< GappedAlignment >& backbone_regions );

/**
 * writes out a list of backbone regions
 */
void outputBackbone( const std::vector< GappedAlignment >& backbone_regions, std::ostream& backbone_out );

void getGapBounds( std::vector<gnSeqI>& seq_lengths, std::vector< LCB >& adjacencies, uint seqJ, int leftI, int rightI, int64& left_start, int64& right_start );


inline
void findRightEndpoint( size_t seqI, size_t seqJ, score_t significance_threshold, std::vector< score_t >& scores, hss_list_t& hss_list )
{
	// see how long it takes score_sum to go to 0, then scan forward to determine where the hss begins
	score_t score_sum = significance_threshold;
	size_t colI = scores.size();
	for( ; colI > 0; --colI )
	{
		if( scores[colI-1] == INVALID_SCORE )
			continue;

		if( score_sum >= 0 && score_sum + scores[colI-1] < 0 )
		{
			// end of an excursion
			score_sum = 0;
			// backtrack to find the MSC in the other direction?
			// call the entire segment between MSCs the HSS?
			score_t rev_score_sum = 0;
			size_t rev_ladder_point = colI-1;
			size_t rcolI = colI-1;
			for( ; ((int64)rcolI) < scores.size(); ++rcolI )
			{
				if( scores[rcolI] == INVALID_SCORE )
					continue;
				if( rev_score_sum > significance_threshold )
					break;
				if( rev_score_sum >= 0 && rev_score_sum + scores[rcolI] < 0 )
				{
					rev_score_sum = 0;
				}else if( rev_score_sum == 0 && scores[rcolI] > 0 )
				{
					// start a new excursion
					rev_score_sum += scores[rcolI];
					rev_ladder_point = rcolI;
				}else
					rev_score_sum += scores[rcolI];
			}
			// the segment between ladder_point and rev_ladder_point is an HSS
			if( rcolI < scores.size() )
			{
				HssCols ic;
				ic.seqI = seqI;
				ic.seqJ = seqJ;
				ic.left_col = rev_ladder_point;
				ic.right_col = scores.size()-1;
				hss_list.push_back( ic );
			}
			break;
		}else
			score_sum += scores[colI-1];
	}
}


inline
void findHssExcursions( std::vector< score_t > scores, score_t significance_threshold, hss_list_t& hss_list, uint seqI, uint seqJ, boolean left_hss, boolean right_hss )
{
	score_t score_sum = left_hss ? significance_threshold : 0;	// start in an hss if non-homologous
	int64 ladder_point = 0;
	bool fwd_hss = left_hss;

	// scan left to right over the columns to identify HSS
	for( size_t colI = 0; colI <= scores.size(); ++colI )
	{
		if( colI < scores.size() && scores[colI] == INVALID_SCORE )
			continue;

		if( colI == scores.size() || (score_sum >= 0 && score_sum + scores[colI] < 0)  )
		{
			// end of an excursion
			if( colI == scores.size() && right_hss )
				fwd_hss = true;

			score_sum = 0;
			if( fwd_hss )
			{
				// call the entire segment between the current column and the ladder point
				// an excursion
				HssCols ic;
				ic.seqI = seqI;
				ic.seqJ = seqJ;
				ic.left_col = ladder_point;
				if( colI == scores.size() )
					ic.right_col = colI - 1;
				else
					ic.right_col = colI;
				hss_list.push_back( ic );
			}
			fwd_hss = false;
		}else if( score_sum == 0 && scores[colI] > 0 )
		{
			// start a new excursion
			score_sum += scores[colI];
			ladder_point = colI;
		}else
			score_sum += scores[colI];

		if( score_sum > significance_threshold )
			fwd_hss = true;
	}
}


inline
void findMscFromExcursions( std::vector< score_t > scores, score_t significance_threshold, hss_list_t& hss_list, hss_list_t& msc_list, uint seqI, uint seqJ, boolean left_hss, boolean right_hss )
{
	score_t left_end_score = left_hss ? significance_threshold : 0;
	score_t right_end_score = right_hss ? significance_threshold : 0;
	score_t score_sum = left_end_score;	// start in an hss if non-homologous
	int64 ladder_point = 0;
	bool fwd_hss = true;
	if( left_hss )
		scores.front() = significance_threshold;
	if( right_hss )
		scores.back() = significance_threshold;

	// for each excursion in hss_list
	for( size_t exI = 0; exI < hss_list.size(); exI++ )
	{
		// create a vector of score sums
		size_t col_base = hss_list[exI].left_col;
		size_t col_count = hss_list[exI].right_col - hss_list[exI].left_col + 1;
		// find cols with positive scores
		size_t positive_count = 0;
		for( size_t colI = hss_list[exI].left_col; colI < hss_list[exI].right_col + 1; colI++ )
		{
			if( scores[colI] <= 0 || scores[colI] == INVALID_SCORE  )
				continue;
			positive_count++;
		}
		std::vector< size_t > pos_map(positive_count);
		positive_count = 0;
		for( size_t colI = hss_list[exI].left_col; colI < hss_list[exI].right_col + 1; colI++ )
		{
			if( scores[colI] <= 0 || scores[colI] == INVALID_SCORE )
				continue;
			pos_map[positive_count] = colI;
			positive_count++;
		}

		std::vector< score_t > score_sums(positive_count, 0);
		size_t sum_base = 0;
		std::vector<HssCols> cur_msc_list;
		size_t invalid_count = 0;
		for( size_t colI = hss_list[exI].left_col; colI < hss_list[exI].right_col + 1; colI++ )
		{
			// skip this column if it has an invalid score
			if( scores[colI] == INVALID_SCORE )
				continue;
			// otherwise add the current column score to all relevant score sums
			int64 msc_col = -1;
			size_t sumI = sum_base;
			for( ; sumI < score_sums.size(); sumI++ )
			{
				// break the loop if the next positive column is past colI
				if( pos_map[sumI] > colI )
					break;
				// don't worry about this one if it's invalid
				if( score_sums[sumI] < 0 )
					continue;
				score_sums[sumI] += scores[colI];
				// if the local score bottoms out to 0 then this one has failed
				if( score_sums[sumI] < 0 )
					invalid_count++;
				// take the right-most starting column which yields an msc
				if( score_sums[sumI] >= significance_threshold )
					msc_col = sumI;
			}
			// did we find a minimum significant cluster?
			if( msc_col != -1 )
			{
				HssCols ic;
				ic.seqI = seqI;
				ic.seqJ = seqJ;
				ic.left_col = pos_map[msc_col];
				ic.right_col = colI;
				cur_msc_list.push_back( ic );
				sum_base = msc_col + 1;	// any new MSC needs to be to the right of this one's left-end
			}
			// or did all of our local sums become invalid?
//			if( invalid_count == sumI - sum_base )
//			{
//				sum_base = sumI + 1;
//			}
		}
		// merge any overlapping MSCs
		size_t disjoint = cur_msc_list.size() > 0 ? 1 : 0;
		for( size_t mscI = 1; mscI < cur_msc_list.size(); mscI++ )
		{
			if( cur_msc_list[mscI].left_col <= cur_msc_list[mscI-1].right_col )
			{
				cur_msc_list[mscI].left_col = cur_msc_list[mscI-1].left_col;
				cur_msc_list[mscI-1].left_col = (std::numeric_limits<size_t>::max)();
				cur_msc_list[mscI-1].right_col = (std::numeric_limits<size_t>::max)();
			}else
				disjoint++;
		}
		std::vector<HssCols> cur_msc_list2( disjoint );
		size_t mI = 0;
		for( size_t mscI = 0; mscI < cur_msc_list.size(); mscI++ )
		{
			if( cur_msc_list[mscI].left_col != (std::numeric_limits<size_t>::max)() )
				cur_msc_list2[mI++] = cur_msc_list[mscI];
		}

		// TODO: grow MSC boundaries to include all surrounding positively scoring regions
		for( mI = 0; mI < cur_msc_list2.size(); mI++ )
		{
			HssCols& hss = cur_msc_list2[mI];
			// first left, then right
			size_t lcolI = hss.left_col + 1;
			for( ; lcolI > 0 && scores[lcolI - 1] > 0; lcolI-- );
			size_t rcolI = hss.right_col;
			for( ; rcolI < scores.size() && scores[rcolI] > 0; rcolI++ );
			hss.left_col = lcolI;
			hss.right_col = rcolI - 1;
		}
		swap( cur_msc_list, cur_msc_list2 );

		// merge overlapping MSCs again
		// BAD: copied code from above

		disjoint = cur_msc_list.size() > 0 ? 1 : 0;
		for( size_t mscI = 1; mscI < cur_msc_list.size(); mscI++ )
		{
			if( cur_msc_list[mscI].left_col <= cur_msc_list[mscI-1].right_col )
			{
				cur_msc_list[mscI].left_col = cur_msc_list[mscI-1].left_col;
				cur_msc_list[mscI-1].left_col = (std::numeric_limits<size_t>::max)();
				cur_msc_list[mscI-1].right_col = (std::numeric_limits<size_t>::max)();
			}else
				disjoint++;
		}
		mI = msc_list.size();
		msc_list.resize( mI + disjoint );
		for( size_t mscI = 0; mscI < cur_msc_list.size(); mscI++ )
		{
			if( cur_msc_list[mscI].left_col != (std::numeric_limits<size_t>::max)() )
				msc_list[mI++] = cur_msc_list[mscI];
		}
	}
}

static char charmap[128];
inline
char* getCharmap()
{
	static bool initialized = false;
	if(initialized)
		return charmap;
	memset(charmap, 0, 128);
	charmap['a'] = 0;
	charmap['c'] = 1;
	charmap['g'] = 2;
	charmap['t'] = 3;
	charmap['-'] = 4;
	charmap['A'] = 0;
	charmap['C'] = 1;
	charmap['G'] = 2;
	charmap['T'] = 3;
	charmap['-'] = 4;
	initialized = true;
	return charmap;
}
// a mapping from pairwise alignment columns to HomologyHMM emission codes
// row/column indices are as given by the charmap above (ACGT- == 01234).
static char colmap[5][5] = {
//    A   C   G   T   -
	{'1','3','4','5','7'},	// A
	{'3','2','6','4','7'},  // C
	{'4','6','2','3','7'},  // G
	{'5','4','3','1','7'},  // T
	{'7','7','7','7','\0'},  // -
};


inline
void findHssHomologyHMM( std::vector< std::string >& aln_table, hss_list_t& hss_list, uint seqI, uint seqJ, double pGoHomo, double pGoUnrelated,
						boolean left_homologous, boolean right_homologous )
{
	static char* charmap = getCharmap();

	// encode the alignment as column states
	std::string column_states(aln_table[0].size(),'q');
	vector< size_t > col_reference(column_states.size(), (std::numeric_limits<size_t>::max)() );
	size_t refI = 0;
	for( size_t colI = 0; colI < column_states.size(); colI++ )
	{
		char a = charmap[aln_table[seqI][colI]];
		char b = charmap[aln_table[seqJ][colI]];
		column_states[colI] = colmap[a][b];
		if(column_states[colI] != 0 )
			col_reference[refI++] = colI;
	}
	// filter out the gap/gap cols
	std::string::iterator sitr = std::remove(column_states.begin(), column_states.end(), 0);
	column_states.resize(sitr - column_states.begin());

	for( size_t colI = 2; colI < column_states.size(); colI++ )
	{
		if( column_states[colI] == '7' &&
			column_states[colI-1] == '7' &&
			(column_states[colI-2] == '7' || column_states[colI-2] == '8') )
			column_states[colI-1] = '8';
	}
	if( column_states.size() > 1 && column_states[0] == '7' && (column_states[1] == '7' || column_states[1] == '8'))
		column_states[0] = '8';
	if( column_states.size() > 1 && column_states[column_states.size()-1] == '7' && (column_states[column_states.size()-2] == '7'|| column_states[column_states.size()-2] == '8') )
		column_states[column_states.size()-1] = '8';
	// now feed it to the Homology prediction HMM
	string prediction;
	if( right_homologous && !left_homologous )
		std::reverse(column_states.begin(), column_states.end());

	run(column_states, prediction, pGoHomo, pGoUnrelated);

	if( right_homologous && !left_homologous )
		std::reverse(prediction.begin(), prediction.end());
	size_t prev_h = 0;
	size_t i = 1;
	for( ; i < prediction.size(); i++ )
	{
		if( prediction[i] == 'H' && prediction[i-1] == 'N' )
		{
			prev_h = i;
		}
		if( prediction[i] == 'N' && prediction[i-1] == 'H' )
		{
			HssCols hc;
			hc.seqI = seqI;
			hc.seqJ = seqJ;
			hc.left_col = col_reference[prev_h];
			hc.right_col = col_reference[i-1];
			hss_list.push_back(hc);
			prev_h = i;
		}
	}
	// get the last one
	if( prediction[i-1] == 'H' )
	{
		HssCols hc;
		hc.seqI = seqI;
		hc.seqJ = seqJ;
		hc.left_col = col_reference[prev_h];
		hc.right_col = col_reference[i-1];
		hss_list.push_back(hc);
	}
}

inline
void findHssRandomWalkScoreVector( std::vector< score_t > scores, score_t significance_threshold, hss_list_t& hss_list, uint seqI, uint seqJ, boolean left_homologous, boolean right_homologous )
{

	score_t left_end_score = left_homologous ? 0 : significance_threshold;
	score_t right_end_score = right_homologous ? 0 : significance_threshold;
	score_t score_sum = left_end_score;	// start in an hss if non-homologous
	score_t lrh = score_sum;
	int64 ladder_point = -1;
	int64 rev_ladder_point = 0;
	bool fwd_hss = !left_homologous;

	for( size_t colI = 0; colI <= scores.size(); ++colI )
	{
		if( colI < scores.size() && scores[colI] == INVALID_SCORE )
			continue;

		if( colI == scores.size() || (score_sum >= 0 && score_sum + scores[colI] < 0) ||
			(score_sum >= lrh - significance_threshold && score_sum + scores[colI] < lrh - significance_threshold ) )
		{
			if( !fwd_hss && colI == scores.size() && !right_homologous )
			{
				fwd_hss = true;
				if( score_sum <= 0 )
					ladder_point = colI - 1;
			}
			// end of an excursion
			score_sum = 0;
			lrh = 0;
			if( fwd_hss )
			{
				// backtrack to find the MSC in the other direction?
				// call the entire segment between MSCs the HSS?
				score_t rev_score_sum = 0;
				if( colI == scores.size() && !right_homologous )
					rev_score_sum = significance_threshold;
				size_t rev_ladder_point = ladder_point;
				if( colI == scores.size() && !right_homologous )
					rev_ladder_point = colI - 1;
				for( size_t rcolI = colI; ((int64)rcolI) > ladder_point && rcolI > 0; --rcolI )
				{
					if( scores[rcolI-1] == INVALID_SCORE )
						continue;
					if( rev_score_sum > significance_threshold )
						break;
					if( rev_score_sum >= 0 && rev_score_sum + scores[rcolI-1] < 0 )
					{
						rev_score_sum = 0;
					}else if( rev_score_sum == 0 && scores[rcolI-1] > 0 )
					{
						// start a new excursion
						rev_score_sum += scores[rcolI-1];
						rev_ladder_point = rcolI-1;
					}else
						rev_score_sum += scores[rcolI-1];

				}
				// don't make an HSS if there was no reverse HSS, unless we
				// ended the excursion artificially because we hit the end of the block...
				if( ((rev_ladder_point != 0 || ladder_point != -1) ||
					colI == scores.size()) && rev_ladder_point != -1 )
				{
					if( colI == scores.size() && ladder_point == -1 )
					{
						rev_ladder_point = scores.size()-1;
					}
					if( ladder_point == -1 )
						ladder_point = 0;
					// the segment between ladder_point and rev_ladder_point is an HSS
					HssCols ic;
					ic.seqI = seqI;
					ic.seqJ = seqJ;
					ic.left_col = ladder_point;
					ic.right_col = rev_ladder_point;
					hss_list.push_back( ic );
				}
			}
			fwd_hss = false;
		}else if( score_sum == 0 && scores[colI] > 0 )
		{
			// start a new excursion
			score_sum += scores[colI];
			ladder_point = colI;
		}else
			score_sum += scores[colI];

		if( score_sum > significance_threshold )
			fwd_hss = true;
		if( score_sum > lrh )
			lrh = score_sum;
	}
}

template< typename MatchVector >
void findHssRandomWalk_v2( const MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table, const PairwiseScoringScheme& scoring, score_t significance_threshold, hss_array_t& hss_array, boolean left_homologous, boolean right_homologous )
{
	typedef typename MatchVector::value_type MatchType;
	if( iv_list.size() == 0 )
		return;
	uint seq_count = seq_table.size();
	hss_array.resize( boost::extents[seq_count][seq_count][iv_list.size()] );
	for( uint iv_listI = 0; iv_listI < iv_list.size(); iv_listI++ ){
		const MatchType& iv = iv_list[ iv_listI ];
		std::vector< std::string > aln_table;
		GetAlignment( *iv, seq_table, aln_table );
		
		for( uint seqI = 0; seqI < seq_count; seqI++ ){
			uint seqJ;
			for( seqJ = seqI + 1; seqJ < seq_count; seqJ++ ){

				hss_list_t& hss_list = hss_array[seqI][seqJ][iv_listI];
				hss_list.clear();

				std::vector< score_t > scores;
				computeMatchScores( aln_table[seqI], aln_table[seqJ], scoring, scores );
				computeGapScores( aln_table[seqI], aln_table[seqJ], scoring, scores );

				// Invert the scores since we're trying to detect rare bouts of non-homologous sequence
				for( size_t sI = 0; sI < scores.size(); ++sI )
					if( scores[sI] != INVALID_SCORE)
						scores[sI] = -scores[sI];

				hss_list_t excursion_list;
				findHssExcursions( scores, significance_threshold, excursion_list, seqI, seqJ, !left_homologous, !right_homologous );
				findMscFromExcursions( scores, significance_threshold, excursion_list, hss_list, seqI, seqJ, !left_homologous, !right_homologous );
			}
		}
	}
}

template< typename MatchVector >
void findHssRandomWalk( const MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table, const PairwiseScoringScheme& scoring, score_t significance_threshold, hss_array_t& hss_array, boolean left_homologous, boolean right_homologous )
{
	typedef typename MatchVector::value_type MatchType;
	if( iv_list.size() == 0 )
		return;
	uint seq_count = seq_table.size();
	hss_array.resize( boost::extents[seq_count][seq_count][iv_list.size()] );
	for( uint iv_listI = 0; iv_listI < iv_list.size(); iv_listI++ ){
		const MatchType& iv = iv_list[ iv_listI ];
		std::vector< std::string > aln_table;
		GetAlignment( *iv, seq_table, aln_table );
		
		for( uint seqI = 0; seqI < seq_count; seqI++ ){
			uint seqJ;
			for( seqJ = seqI + 1; seqJ < seq_count; seqJ++ ){

				hss_list_t& hss_list = hss_array[seqI][seqJ][iv_listI];
				hss_list.clear();

				std::vector< score_t > scores;
				computeMatchScores( aln_table[seqI], aln_table[seqJ], scoring, scores );
				computeGapScores( aln_table[seqI], aln_table[seqJ], scoring, scores );

				// Invert the scores since we're trying to detect rare bouts of non-homologous sequence
				for( size_t sI = 0; sI < scores.size(); ++sI )
					if( scores[sI] != INVALID_SCORE)
						scores[sI] = -scores[sI];
				findHssRandomWalkScoreVector( scores, significance_threshold, hss_list, seqI, seqJ, left_homologous, right_homologous );
			}
		}
	}
}

template< typename MatchVector >
void findHssHomologyHMM( const MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table, const PairwiseScoringScheme& scoring, hss_array_t& hss_array, double pGoHomo, double pGoUnrelated, boolean left_homologous, boolean right_homologous )
{
	typedef typename MatchVector::value_type MatchType;
	if( iv_list.size() == 0 )
		return;
	uint seq_count = seq_table.size();
	hss_array.resize( boost::extents[seq_count][seq_count][iv_list.size()] );
	for( uint iv_listI = 0; iv_listI < iv_list.size(); iv_listI++ ){
		const MatchType& iv = iv_list[ iv_listI ];
		std::vector< std::string > aln_table;
		GetAlignment( *iv, seq_table, aln_table );
		
		for( uint seqI = 0; seqI < seq_count; seqI++ ){
			uint seqJ;
			for( seqJ = seqI + 1; seqJ < seq_count; seqJ++ ){

				hss_list_t& hss_list = hss_array[seqI][seqJ][iv_listI];
				hss_list.clear();
				findHssHomologyHMM( aln_table, hss_list, seqI, seqJ, pGoHomo, pGoUnrelated, left_homologous, right_homologous );
			}
		}
	}
}


template< typename MatchVector >
void HssColsToIslandCols( const MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table, hss_array_t& hss_array, hss_array_t& island_col_array )
{

	typedef typename MatchVector::value_type MatchType;
	uint seq_count = seq_table.size();
	island_col_array.resize( boost::extents[seq_count][seq_count][iv_list.size()] );
	for( uint iv_listI = 0; iv_listI < iv_list.size(); iv_listI++ ){
		const MatchType& iv = iv_list[ iv_listI ];
		for( uint seqI = 0; seqI < seq_count; seqI++ ){
			uint seqJ;
			for( seqJ = seqI + 1; seqJ < seq_count; seqJ++ ){
				hss_list_t& hss_list = hss_array[seqI][seqJ][iv_listI];
				hss_list_t& island_col_list = island_col_array[seqI][seqJ][iv_listI];
				ComplementHss(iv_list[iv_listI]->AlignmentLength(),hss_list,island_col_list,seqI,seqJ);
			}
		}
	}
}
inline
void ComplementHss( const size_t alignment_length, hss_list_t& hss_list, hss_list_t& island_col_list, uint seqI=0, uint seqJ=0 )
{


	size_t left_col = 0;
	for( size_t hssI = 0; hssI < hss_list.size(); ++hssI )
	{
		if( left_col >= hss_list[hssI].left_col ) 
		{
			left_col = hss_list[hssI].right_col + 1;
			continue;	// handle the case where the HSS starts at col 0
		}
		// ending an island
		IslandCols isle;
		isle.seqI = seqI;
		isle.seqJ = seqJ;
		isle.left_col = left_col;
		isle.right_col = hss_list[hssI].left_col;
		island_col_list.push_back(isle);
		left_col = hss_list[hssI].right_col + 1;
	}

	if( left_col < alignment_length )
	{
		// add the last island
		IslandCols isle;
		isle.seqI = seqI;
		isle.seqJ = seqJ;
		isle.left_col = left_col;
		isle.right_col = alignment_length-1;
		island_col_list.push_back(isle);
	}
}

template< typename MatchVector >
void HssArrayToCga( const MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table, hss_array_t& hss_array, std::vector< CompactGappedAlignment<>* >& cga_list )
{
	typedef typename MatchVector::value_type MatchType;
	uint seq_count = seq_table.size();
	for( uint iv_listI = 0; iv_listI < iv_list.size(); iv_listI++ ){
		const MatchType& iv = iv_list[ iv_listI ];
		
		CompactGappedAlignment<>* iv_cga = dynamic_cast< CompactGappedAlignment<>* >(iv);
		bool allocated = false;
		if( iv_cga == NULL )
		{
			CompactGappedAlignment<> tmp_cga;
			iv_cga = tmp_cga.Copy();
			new (iv_cga) CompactGappedAlignment<>(*iv);
			allocated = true;
		}
		for( uint seqI = 0; seqI < seq_count; seqI++ ){
			for( uint seqJ = seqI + 1; seqJ < seq_count; seqJ++ ){
				hss_list_t& isle_list = hss_array[seqI][seqJ][iv_listI];
				for( size_t curI = 0; curI < isle_list.size(); ++curI )
				{
					// extract a cga
					CompactGappedAlignment<> tmp_cga;
					cga_list.push_back( tmp_cga.Copy() );
					iv_cga->copyRange( *(cga_list.back()), isle_list[curI].left_col, isle_list[curI].right_col - isle_list[curI].left_col + 1 );
					if( cga_list.back()->LeftEnd(0) == NO_MATCH )
					{
						// this one must have been covering an invalid region (gaps aligned to gaps)
						cga_list.back()->Free();
						cga_list.erase( cga_list.end()-1 );
					}
				}
			}
		}
		if( allocated )
			iv_cga->Free();
	}
}

template< typename MatchVector >
void findHssRandomWalkCga( const MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table, const PairwiseScoringScheme& scoring, score_t significance_threshold, std::vector< CompactGappedAlignment<>* >& hss_list )
{
	if( iv_list.size() == 0 )
		return;
	hss_array_t hss_array;
	findHssRandomWalk( iv_list, seq_table, scoring, significance_threshold, hss_array );
	hss_array_t homo_array;
	HssColsToIslandCols( iv_list, seq_table, hss_array, homo_array );
	HssArrayToCga(iv_list, seq_table, homo_array, hss_list);
}
template< typename MatchVector >
void findIslandsRandomWalkCga( const MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table, const PairwiseScoringScheme& scoring, score_t significance_threshold, std::vector< CompactGappedAlignment<>* >& island_list )
{
	if( iv_list.size() == 0 )
		return;
	hss_array_t hss_array;
	findHssRandomWalk( iv_list, seq_table, scoring, significance_threshold, hss_array );
	hss_array_t island_col_array;
//	HssColsToIslandCols( iv_list, hss_array, island_col_array );
	HssArrayToCga(iv_list, seq_table, hss_array, island_list);
}

template< typename MatchVector >
void findIslandsRandomWalk( const MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table, const PairwiseScoringScheme& scoring, score_t significance_threshold, std::vector< Island >& island_list )
{
	if( iv_list.size() == 0 )
		return;
	hss_array_t hss_array;
	findHssRandomWalk( iv_list, seq_table, scoring, significance_threshold, hss_array );

	typedef typename MatchVector::value_type MatchType;
	for( uint iv_listI = 0; iv_listI < iv_list.size(); iv_listI++ ){
		const MatchType& iv = iv_list[ iv_listI ];
		uint seq_count = seq_table.size();
		std::vector< std::string > aln_table;
		GetAlignment( *iv, seq_table, aln_table );
		
		for( uint seqI = 0; seqI < seq_count; seqI++ ){
			uint seqJ;
			for( seqJ = seqI + 1; seqJ < seq_count; seqJ++ ){
				hss_list_t& hss_list = hss_array[seqI][seqJ][iv_listI];

				size_t cur_hss = 0;
				uint columnI = 0;
				gnSeqI curI = 0;
				gnSeqI curJ = 0;
				gnSeqI lastI = 0;
				gnSeqI lastJ = 0;
				for( columnI = 0; columnI < aln_table[0].size(); columnI++ ){
					if( aln_table[ seqI ][ columnI ] != '-' )
						curI++;
					if( aln_table[ seqJ ][ columnI ] != '-' )
						curJ++;

					if( cur_hss < hss_list.size() && 
						columnI == hss_list[cur_hss].left_col )
					{
						// ending an island
						int64 leftI = iv.Start( seqI );
						int64 rightI = leftI < 0 ? leftI - curI : leftI + curI;
						leftI = leftI < 0 ? leftI - lastI : leftI + lastI;
						int64 leftJ = iv.Start( seqJ );
						int64 rightJ = leftJ < 0 ? leftJ - curJ : leftJ + curJ;
						leftJ = leftJ < 0 ? leftJ - lastJ : leftJ + lastJ;
						Island isle;
						isle.seqI = seqI;
						isle.seqJ = seqJ;
						isle.leftI = leftI;
						isle.leftJ = leftJ;
						isle.rightI = rightI;
						isle.rightJ = rightJ;
						island_list.push_back(isle);
					}
					else if( cur_hss < hss_list.size() &&
						columnI == hss_list[cur_hss].right_col )
					{
						// starting an island
						lastI = curI;
						lastJ = curJ;
						cur_hss++;
					}
				}

				// add the last island
				int64 leftI = iv.Start( seqI );
				int64 rightI = leftI < 0 ? leftI - curI : leftI + curI;
				leftI = leftI < 0 ? leftI - lastI : leftI + lastI;
				int64 leftJ = iv.Start( seqJ );
				int64 rightJ = leftJ < 0 ? leftJ - curJ : leftJ + curJ;
				leftJ = leftJ < 0 ? leftJ - lastJ : leftJ + lastJ;
				Island isle;
				isle.seqI = seqI;
				isle.seqJ = seqJ;
				isle.leftI = leftI;
				isle.leftJ = leftJ;
				isle.rightI = rightI;
				isle.rightJ = rightJ;
				island_list.push_back(isle);
			}
		}
	}
}

template< class IntervalListType >
void addUnalignedRegions( IntervalListType& iv_list)
{
	std::vector< AbstractMatch* > new_ivs;
	std::vector< AbstractMatch* > iv_ptrs(iv_list.size());
	for( size_t i = 0; i < iv_list.size(); ++i )
		iv_ptrs[i] = &iv_list[i];
	for( size_t seqI = 0; seqI < iv_list.seq_table.size(); ++seqI )
	{
		SingleStartComparator< AbstractMatch > ssc( seqI );
		std::sort( iv_ptrs.begin(), iv_ptrs.end(), ssc );
		size_t ivI = 0;
		for( ; ivI < iv_ptrs.size(); ++ivI )
			if( iv_ptrs[ivI]->LeftEnd(seqI) != NO_MATCH )
				break;
		std::list< AbstractMatch* > iv_ptr_list;
		iv_ptr_list.insert( iv_ptr_list.end(), iv_ptrs.begin()+ivI, iv_ptrs.end() );
		AddGapMatches( iv_ptr_list, iv_ptr_list.begin(), iv_ptr_list.end(), seqI, 1, iv_list.seq_table[seqI]->length()+1, AbstractMatch::forward, iv_list.seq_table.size() );
		std::list< AbstractMatch* >::iterator iter = iv_ptr_list.begin();
		while( ivI != iv_ptrs.size() && iter != iv_ptr_list.end() )
		{
			if( iv_ptrs[ivI] == *iter )
				ivI++;
			else
				new_ivs.push_back( *iter );
			++iter;
		}
		while( iter != iv_ptr_list.end() )
		{
			new_ivs.push_back( *iter );
			++iter;
		}
	}
	// now add all the new intervals to iv_list
	size_t prev_size = iv_list.size();
	iv_list.resize( iv_list.size() + new_ivs.size() );
	for( size_t newI = 0; newI < new_ivs.size(); ++newI )
	{
		Interval iv( new_ivs.begin() + newI, new_ivs.begin() + newI + 1 );
		iv_list[prev_size + newI] = iv;
	}
}


}

#endif // __Islands_h__

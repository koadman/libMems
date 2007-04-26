/*******************************************************************************
 * $Id: Islands.h,v 1.7 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
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

void allocateDetailList( IntervalList& iv_list, std::vector< std::pair< uint64, uint64 > >& detail_list );
void allocateDetailList( MatchList& iv_list, std::vector< std::pair< uint64, uint64 > >& detail_list );

void getLCBDetailList( Interval& iv, std::vector< std::pair< uint64, uint64 > >& detail_list );

void getLCBDetailList( MatchList& iv, std::vector< std::pair< uint64, uint64 > >& detail_list );

void PrintDetailList( uint seq_count, const std::vector< std::pair< uint64, uint64 > >& detail_list, std::ostream& os, boolean skip_zeros = true);
//void PrintDetailList( uint seq_count, vector< pair< uint64, uint64 > >& detail_list, ostream& os, boolean skip_zeros );


void getGapBounds( std::vector<gnSeqI>& seq_lengths, std::vector< LCB >& adjacencies, uint seqJ, int leftI, int rightI, int64& left_start, int64& right_start );

void TransformDistanceIdentity( NumericMatrix<double>& identity );

void DistanceMatrix( const MatchList& mlist, NumericMatrix<double>& identity );


template< class AbstractMatchVectorType >
void IdentityMatrix( const AbstractMatchVectorType& matches, const std::vector< genome::gnSequence* >& seq_table, NumericMatrix<double>& identity );
template<class AbstractMatchType>
void MatchIdentityMatrix( const AbstractMatchType& amt, const std::vector< genome::gnSequence* >& seq_table, NumericMatrix<double>& identity);

void DistanceMatrix( uint seq_count, const std::vector< std::pair< uint64, uint64 > >& detail_list, NumericMatrix<double>& distance );

void IdentityMatrix( const IntervalList& iv_list, NumericMatrix<double>& identity );
inline
void IdentityMatrix( const IntervalList& iv_list, NumericMatrix<double>& identity )
{
	std::vector< const AbstractMatch* > am_list;
	for( size_t ivI = 0; ivI < iv_list.size(); ivI++ )
		am_list.push_back( &iv_list[ivI] );
	IdentityMatrix( am_list, iv_list.seq_table, identity );
}

template< class AbstractMatchVectorType >
void IdentityMatrix( const AbstractMatchVectorType& matches, const std::vector< genome::gnSequence* >& seq_table, NumericMatrix<double>& identity ){
	if( matches.size() == 0 )
		return;

	uint seq_count = seq_table.size();
	identity = NumericMatrix<double>( seq_count, seq_count );
	identity.init( 0 );
	NumericMatrix<double> possible( seq_count, seq_count );
	possible.init( 0 );
	
	for( uint ivI = 0; ivI < matches.size(); ivI++ ){
		AddToMatchIdentityMatrix( *matches[ ivI ], seq_table, identity );
	}
	for( uint seqI = 0; seqI < seq_count; seqI++ ){
		for( uint seqJ = 0; seqJ < seq_count; seqJ++ ){
			gnSeqI shorter_len = seq_table[seqI]->length() < seq_table[seqJ]->length() ? seq_table[seqI]->length() : seq_table[seqJ]->length();
			possible( seqI, seqJ ) += shorter_len;
		}
	}
	identity /= possible;
}


template< class AbstractMatchVectorType >
void BackboneIdentityMatrix( const AbstractMatchVectorType& matches, const std::vector< genome::gnSequence* >& seq_table, NumericMatrix<double>& identity ){
	if( matches.size() == 0 )
		return;

	size_t seq_count = seq_table.size();
	identity = NumericMatrix<double>( seq_count, seq_count );
	identity.init( 0 );
	
	for( uint ivI = 0; ivI < matches.size(); ivI++ ){
		AddToMatchIdentityMatrix( *matches[ ivI ], seq_table, identity );
	}

	NumericMatrix<double> possible( seq_count, seq_count );
	possible.init( 0 );

	for( size_t mI = 0; mI < matches.size(); ++mI ){
		std::vector< std::string > alignment;
		GetAlignment( *(matches[mI]), seq_table, alignment );
		for( gnSeqI charI = 0; charI < matches[mI]->AlignmentLength(); charI++ ){
			for( size_t seqI = 0; seqI < seq_count; seqI++ ){
				for( size_t seqJ = seqI + 1; seqJ < seq_count; seqJ++ ){
					if( alignment[ seqI ][ charI ] != '-' &&
						alignment[ seqJ ][ charI ] != '-' ){
							possible( seqI, seqJ ) += 1;
					}
				}
			}
		}
	}

	identity /= possible;
}


template<class AbstractMatchType>
void MatchIdentityMatrix( const AbstractMatchType& amt, const std::vector< genome::gnSequence* >& seq_table, NumericMatrix<double>& identity)
{
	if( amt.SeqCount() == 0 )
		return;
	uint seq_count = amt.SeqCount();
	identity = NumericMatrix<double>( seq_count, seq_count );
	identity.init( 0 );
	uint seqI;
	uint seqJ;

	std::vector< std::string > alignment;
	GetAlignment( amt, seq_table, alignment );
	for( gnSeqI charI = 0; charI < amt.AlignmentLength(); charI++ ){
		for( seqI = 0; seqI < seq_count; seqI++ ){
			for( seqJ = seqI + 1; seqJ < seq_count; seqJ++ ){
				if( ( toupper( alignment[ seqI ][ charI ] ) == 
					toupper( alignment[ seqJ ][ charI ] ) ) &&
					alignment[ seqI ][ charI ] != '-' ){
					
						identity( seqI, seqJ ) += 1;
				}
			}
		}
	}

	for( seqI = 0; seqI < seq_count; seqI++ ){
		for( seqJ = seq_count; seqJ > 0; seqJ-- ){
			if( seqI == seqJ - 1 )
				// set the diagonal to identical
				identity( seqI, seqJ - 1 ) = 1;
			else if( seqI < seqJ - 1 ){
				// determine the length of the shorter sequence
				gnSeqI shorter_len = amt.Length( seqI ) < amt.Length( seqJ - 1 ) ? amt.Length( seqI ) : amt.Length( seqJ - 1 );
				// divide through
				identity( seqI, seqJ - 1 ) /= (double)shorter_len;
				// maxes out at 1
				if( identity( seqI, seqJ - 1 ) > 1 )
					identity( seqI, seqJ - 1 ) = 1;
			}else	// copy the other one
				identity( seqI, seqJ - 1 ) = identity( seqJ - 1, seqI );
		}
	}
}



template<class AbstractMatchType>
void AddToMatchIdentityMatrix( const AbstractMatchType& amt, const std::vector< genome::gnSequence* >& seq_table, NumericMatrix<double>& identity)
{
	if( amt.SeqCount() == 0 )
		return;
	uint seq_count = amt.SeqCount();
	uint seqI;
	uint seqJ;

	std::vector< std::string > alignment;
	GetAlignment( amt, seq_table, alignment );
	for( gnSeqI charI = 0; charI < amt.AlignmentLength(); charI++ ){
		for( seqI = 0; seqI < seq_count; seqI++ ){
			for( seqJ = seqI + 1; seqJ < seq_count; seqJ++ ){
				if( ( toupper( alignment[ seqI ][ charI ] ) == 
					toupper( alignment[ seqJ ][ charI ] ) ) &&
					alignment[ seqI ][ charI ] != '-' ){
					
						identity( seqI, seqJ ) += 1;
				}
			}
		}
	}
}

/*
// template specialization for (exact) matches
inline
void AddToMatchIdentityMatrix( const Match& m, const std::vector< genome::gnSequence* >& seq_table, NumericMatrix<double>& identity)
{
	if( m.SeqCount() == 0 )
		return;
	for( uint seqI = 0; seqI < m.SeqCount(); seqI++ )
		if( m.LeftEnd(seqI) != NO_MATCH )
			for( uint seqJ = seqI + 1; seqJ < m.SeqCount(); seqJ++ )
				if( m.LeftEnd(seqJ) != NO_MATCH )
					identity(seqI,seqJ) += m.Length();
}
*/

template< typename MatchVector >
void SingleCopyDistanceMatrix( MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table, NumericMatrix<double>& distance )
{
	uint seq_count = seq_table.size();
	distance = NumericMatrix<double>( seq_count, seq_count );
	distance.init( 0 );
	uint seqI;
	uint seqJ;
	std::vector< std::pair< bitset_t, bitset_t > > tmp_comp( seq_count );
	std::vector< std::vector< std::pair< bitset_t, bitset_t > > > pair_comp( seq_count, tmp_comp );
	for( uint seqI = 0; seqI < seq_count; ++seqI )
	{
		for( uint seqJ = seqI+1; seqJ < seq_count; ++seqJ )
		{
			pair_comp[seqI][seqJ].first.resize( seq_table[seqI]->length(), false );
			pair_comp[seqI][seqJ].second.resize( seq_table[seqJ]->length(), false );
		}
	}
	for( size_t ivI = 0; ivI < iv_list.size(); ++ivI )
	{
		std::vector< bitset_t > aln_table;
		iv_list[ivI]->GetAlignment(aln_table);
		for( uint seqI = 0; seqI < seq_count; ++seqI )
		{
			for( uint seqJ = seqI+1; seqJ < seq_count; ++seqJ )
			{
				gnSeqI seqI_pos = iv_list[ivI]->LeftEnd(seqI);
				gnSeqI seqJ_pos = iv_list[ivI]->LeftEnd(seqJ);
				AbstractMatch::orientation o_i = iv_list[ivI]->Orientation(seqI);
				AbstractMatch::orientation o_j = iv_list[ivI]->Orientation(seqJ);
				if( o_i == AbstractMatch::reverse )
					seqI_pos = iv_list[ivI]->RightEnd(seqI);
				if( o_j == AbstractMatch::reverse )
					seqJ_pos = iv_list[ivI]->RightEnd(seqJ);
				if( seqI_pos == NO_MATCH || seqJ_pos == NO_MATCH )
					continue;
				for( size_t colI = 0; colI < aln_table[seqI].size(); ++colI )
				{
					if( aln_table[seqI].test(colI) && aln_table[seqJ].test(colI) )
					{
						pair_comp[seqI][seqJ].first.set(seqI_pos-1,true);
						pair_comp[seqI][seqJ].second.set(seqJ_pos-1,true);
					}
					if( aln_table[seqI].test(colI) )
						if( o_i == AbstractMatch::forward )
							seqI_pos++;
						else
							seqI_pos--;
					if( aln_table[seqJ].test(colI) )
						if( o_j == AbstractMatch::forward )
							seqJ_pos++;
						else
							seqJ_pos--;
				}
			}
		}
	}
	for( uint seqI = 0; seqI < seq_count; ++seqI )
	{
		for( uint seqJ = seqI+1; seqJ < seq_count; ++seqJ )
		{
			double pI = ((double)pair_comp[seqI][seqJ].first.count())/((double)pair_comp[seqI][seqJ].first.size());
			double pJ = ((double)pair_comp[seqI][seqJ].second.count())/((double)pair_comp[seqI][seqJ].second.size());
			distance(seqI,seqJ) = (pI + pJ) / 2.0;
			distance(seqJ,seqI) = (pI + pJ) / 2.0;
		}
	}
	TransformDistanceIdentity(distance);
}

static const score_t INVALID_SCORE = (std::numeric_limits<score_t>::max)();

//tjtaed: function to compute the SP column score, and cumulative SP score from an alignment
void computeSPScore( const std::vector<std::string>& alignment, const PairwiseScoringScheme& pss, std::vector<score_t>& scores, score_t& score );
//tjt: function to compute the consensus column score, consensus sequence, and cumulative consensus score from an alignment
void computeConsensusScore( const std::vector<std::string>& alignment, const PairwiseScoringScheme& pss, std::vector<score_t>& scores, std::string& consensus, score_t& score );
void computeMatchScores( const std::string& seq1, const std::string& seq2, const PairwiseScoringScheme& scoring, std::vector<score_t>& scores );
void computeGapScores( const std::string& seq1, const std::string& seq2, const PairwiseScoringScheme& scoring, std::vector<score_t>& scores );

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

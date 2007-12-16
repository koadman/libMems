/*******************************************************************************
 * $Id: Backbone.h,v 1.7 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __Backbone_h__
#define __Backbone_h__

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
#include "libMems/Islands.h"
#include <boost/multi_array.hpp>

#include <sstream>
#include <vector>

namespace mems {

// this is a 99.9% score threshold derived from the EVD of
// simulations of homolgous sequence diverged to .75 substitutions per site and .05 indels per site
const mems::score_t DEFAULT_ISLAND_SCORE_THRESHOLD = 2727;

typedef mems::UngappedLocalAlignment< mems::HybridAbstractMatch<> > ULA;
typedef std::vector< std::vector< ULA* > > backbone_list_t;

/**
 * collapse Intervals that are trivially collinear with each other
 */
void collapseCollinear( IntervalList& iv_list );

/**
 * sanity checks for alignment columns that contain only gaps
 */
void checkForAllGapColumns( IntervalList& iv_list );

/**
 * Applies pairwise transitive homology statistics to detect backbone in a single collinear alignment
 * Unaligns any regions found to be non-homologous, returns coordinates of the homologous segments in bb_list
 * @param	m			The input match in which homology detection will be applied
 * @param	seq_table	A sequence table with one gnSequence pointer per match component
 * @param	result		(output) A newly allocated CompactGappedAlignment that contains the resulting alignment of 
 *						homologous sequence.  It is the caller's responsibility to free the memory using AbstractMatch::Free()
 * @param	bb_list		(output) A list of homologous segments among each component of the output match
 * @param	subst_scoring	The pairwise scoring scheme to apply
 * @param	score_threshold	The significance threshold for score drops that will indicate a transition 
 *							from homology to non-homology
 * @param	pGoHomo	        Unrelated to Homologous transition parameter
 * @param	pGoUnrelated	Homologous to Unrelated transition parameter
 * @param	left_homologous	Set to true if the detection code should assume that sequence beyond the left-most alignment
 *							column is homologous sequence
 * @param	right_homologous	Set to true if the detection code should assume that sequence beyond the right-most alignment
 *							column is homologous sequence
 */
void detectAndApplyBackbone( AbstractMatch* m, std::vector< genome::gnSequence* >& seq_table, CompactGappedAlignment<>*& result, backbone_list_t& bb_list, const PairwiseScoringScheme& subst_scoring, score_t score_threshold = DEFAULT_ISLAND_SCORE_THRESHOLD, const float pGoHomo = 0.004, const float pGoUnrelated = 0.004, std::vector<double>* pEmitHomo = NULL, std::vector<double>* pEmitUnrelated = NULL, boolean left_homologous = false, boolean right_homologous = false );

/**
 * Applies pairwise transitive homology statistics to detect backbone in a genome alignment
 * Unaligns any regions found to be non-homologous, returns coordinates of the homologous segments in bb_list
 */
void detectAndApplyBackbone( IntervalList& iv_list, backbone_list_t& bb_list, const PairwiseScoringScheme& subst_scoring, double pGoHomo, double pGoUnrelated, std::vector<double>* pEmitHomo = NULL, std::vector<double>* pEmitUnrelated = NULL,score_t score_threshold = DEFAULT_ISLAND_SCORE_THRESHOLD );

/**
 * Writes a backbone column file.  This file type gets used by the Mauve GUI.
 */
void writeBackboneColumns( std::ostream& bb_out, backbone_list_t& bb_list );

/**
 * Writes a backbone sequence coordinate file.  This file type is easier to analyze with statistical packages.
 */
void writeBackboneSeqCoordinates( backbone_list_t& bb_list, IntervalList& iv_list, std::ostream& bb_out );



typedef std::vector< std::pair< int64, int64 > > bb_seqentry_t;
typedef struct bb_entry_s
{
	bb_seqentry_t bb_seq;
	ULA bb_cols;
	size_t iv;
} bb_entry_t;

inline
void printBbSeq( std::ostream& os, const bb_seqentry_t& bbseq )
{
	for( size_t i = 0; i < bbseq.size(); ++i )
	{
		if( i > 0 )
			os << '\t';
		os << "(" << bbseq[i].first << ", " << bbseq[i].second << ")";
	}
}

inline
void readBackboneSeqFile( std::istream& bbseq_input, std::vector< bb_seqentry_t >& backbone )
{
	std::string cur_line;
	std::getline( bbseq_input, cur_line );	// read off the header line
	while( std::getline( bbseq_input, cur_line ) )
	{
		bb_seqentry_t bb;
		std::stringstream line_str( cur_line );
		int64 lpos = 0;
		while( line_str >> lpos )
		{
			int64 rpos = 0;
			line_str >> rpos;
			bb.push_back( std::make_pair( lpos, rpos ) );
		}
		backbone.push_back(bb);
	}
}

inline
void readBackboneColsFile( std::istream& bbcol_input, std::vector< std::pair< size_t, ULA > >& bb_list )
{
	std::string cur_line;
	while( std::getline( bbcol_input, cur_line ) )
	{
		ULA tmp_ula;
		size_t ivI;
		std::stringstream ss( cur_line );
		ss >> ivI;
		size_t left_col;
		size_t len;
		ss >> left_col;
		ss >> len;
		gnSeqI bbseq;
		while( ss >> bbseq )
		{
			tmp_ula.SetStart( bbseq, left_col );
		}
		tmp_ula.SetLength( len );
		bb_list.push_back( std::make_pair( ivI, tmp_ula ) );
	}
}


}

#endif	// __Backbone_h__


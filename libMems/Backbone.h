/*******************************************************************************
 * $Id: Backbone.h,v 1.7 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
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
 */
void detectAndApplyBackbone( AbstractMatch* m, std::vector< genome::gnSequence* >& seq_table, CompactGappedAlignment<>*& result, backbone_list_t& bb_list, const PairwiseScoringScheme& subst_scoring, score_t score_threshold = DEFAULT_ISLAND_SCORE_THRESHOLD );

/**
 * Applies pairwise transitive homology statistics to detect backbone in a genome alignment
 * Unaligns any regions found to be non-homologous, returns coordinates of the homologous segments in bb_list
 */
void detectAndApplyBackbone( IntervalList& iv_list, backbone_list_t& bb_list, const PairwiseScoringScheme& subst_scoring, score_t score_threshold = DEFAULT_ISLAND_SCORE_THRESHOLD );

/**
 * Writes a backbone column file.  This file type gets used by the Mauve GUI.
 */
void writeBackboneColumns( std::ostream& bb_out, backbone_list_t& bb_list );

/**
 * Writes a backbone sequence coordinate file.  This file type is easier to analyze with statistical packages.
 */
void writeBackboneSeqCoordinates( backbone_list_t& bb_list, IntervalList& iv_list, std::ostream& bb_out );

}

#endif	// __Backbone_h__


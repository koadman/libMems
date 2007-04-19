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
 */
void detectAndApplyBackbone( AbstractMatch* m, std::vector< genome::gnSequence* >& seq_table, backbone_list_t& bb_list, const PairwiseScoringScheme& subst_scoring, score_t score_threshold );

/**
 * Applies pairwise transitive homology statistics to detect backbone in a genome alignment
 * Unaligns any regions found to be non-homologous, returns coordinates of the homologous segments in bb_list
 */
void detectAndApplyBackbone( IntervalList& iv_list, backbone_list_t& bb_list, const PairwiseScoringScheme& subst_scoring, score_t score_threshold );

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


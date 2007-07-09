#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef __GreedyBreakpointElimination_h__
#define __GreedyBreakpointElimination_h__

#include <libMems/AbstractMatch.h>
#include <iostream>
#include <boost/multi_array.hpp>
#include <libMems/PhyloTree.h>
#include <libMems/SubstitutionMatrix.h>
#include <libMems/SeedOccurrenceList.h>
#include <libMems/IntervalList.h>
#include <libMems/LCB.h>

namespace mems {

/**
 * A wrapper that maps a match among extant sequences to a match among ancestral and extant seqs
 */
template <class MatchType>
class LcbTrackingMatch
{ 
public:
	MatchType original_match;
	MatchType node_match;
	size_t match_id;	// used to index into global arrays of lcb_id and score
};
typedef LcbTrackingMatch< mems::AbstractMatch* > TrackingMatch;

/** 
 * This class is used to track relationships between LCBs during the LCB determination process.
 */
template <class MatchType>
class TrackingLCB
{
public:
	TrackingLCB(){}
	TrackingLCB( const TrackingLCB& l ){ *this = l; }
	/** Constructs a TrackingLCB from a pairwise LCB */
	TrackingLCB( const mems::LCB& l ){ *this = l; }
	TrackingLCB& operator=( const mems::LCB& l )
	{
		left_end[0] = l.left_end[0];
		left_end[1] = l.left_end[1];
		right_end[0] = l.right_end[0];
		right_end[1] = l.right_end[1];
		left_adjacency[0] = l.left_adjacency[0];
		left_adjacency[1] = l.left_adjacency[1];
		right_adjacency[0] = l.right_adjacency[0];
		right_adjacency[1] = l.right_adjacency[1];
		lcb_id = l.lcb_id;
		weight = l.weight;
		to_be_deleted = false;
		return *this;
	}
	int64 left_end[2];	/**< The left end position of the LCB in each sequence */
	int64 right_end[2];  /**< The right end position of the LCB in each sequence */
	uint left_adjacency[2];	/**< 'Pointers' (actually IDs) to the LCBs on the left in each sequence */
	uint right_adjacency[2];	/**< 'Pointers' (actually IDs) to the LCBs on the right in each sequence */
	double weight;		/**< The weight (or coverage) of this LCB */
	std::vector< MatchType > matches;
	int lcb_id;			/**< A numerical ID that can be assigned to this LCB */
	bool to_be_deleted;
};

/** indicates an LCB identifier hasn't been assigned or is unknown */
const uint LCB_UNASSIGNED = (std::numeric_limits<uint>::max)();

typedef boost::multi_array< std::vector< TrackingLCB< TrackingMatch* > >, 2 > PairwiseLCBMatrix;


/**
 * computes an anchoring score for the matches contained inside an LCB
 */
template< class MatchVector >
double GetPairwiseAnchorScore( 
		MatchVector& lcb, std::vector< genome::gnSequence* >& seq_table, 
		const mems::PairwiseScoringScheme& subst_scoring, mems::SeedOccurrenceList& sol_1, 
		mems::SeedOccurrenceList& sol_2, bool penalize_gaps = false );


/**
 * Computes all pairwise LCBs from a set of tracking matches
 */
void getPairwiseLCBs( 
	uint nI, 
	uint nJ, 
	uint dI, 
	uint dJ, 
	std::vector< TrackingMatch* >& tracking_matches, 
	std::vector< TrackingLCB<TrackingMatch*> >& t_lcbs,
	boost::multi_array< double, 3 >& tm_score_array,
	boost::multi_array< size_t, 3 >& tm_lcb_id_array );

/** creates an appropriately sized matrix for mapping individual TrackingMatches to their containing LCBs */
void initTrackingMatchLCBTracking( 
  const std::vector< mems::TrackingMatch >& tracking_matches, 
	size_t n1_count, 
	size_t n2_count, 
	boost::multi_array< size_t, 3 >& tm_lcb_id_array );


/** removes an LCB from an LCB list and coalesces surrounding LCBs.  Returns the number of LCBs removed 
 *  After LCBs are removed, the adjacency list should be processed with filterLCBs()
 *  @param	id_remaps	This is populated with a list of LCB ids that were deleted or coalesced and now have a new LCB id
 *                      for each coalesced LCB, an entry of the form <old id, new id> is added, deleted LCBs have
 *						entries of the form <deleted, -1>.  Entries appear in the order operations were performed
 *						and the function undoLcbRemoval() can undo these operations in reverse order
 */
template< class LcbVector >
uint RemoveLCBandCoalesce( size_t lcbI, uint seq_count, 
						  LcbVector& adjacencies, 
						  std::vector< double >& scores, 
						  std::vector< std::pair< uint, uint > >& id_remaps, 
						  std::vector< uint >& impact_list );


void printMatch( mems::AbstractMatch* m, std::ostream& os );

inline
void printMatch( mems::AbstractMatch* m, std::ostream& os )
{
	for( size_t ii = 0; ii < m->SeqCount(); ++ii )
	{
		if( ii > 0 )
			os << '\t';
		os << "(" << m->Start(ii) << "," << m->RightEnd(ii) << ")";
	}
}

void printProgress( uint prev_prog, uint cur_prog, std::ostream& os );


template< typename PairType >
class LabelSort 
{
public:
	LabelSort( uint seqI ) : ssc( seqI ) {};
	bool operator()( const PairType& pt1, const PairType& pt2 )
	{
		return ssc( pt1.first, pt2.first );
	}
private:
	LabelSort();
	mems::SSC<mems::AbstractMatch> ssc;
};

template<class MatchVector>
void IdentifyBreakpoints( MatchVector& mlist, std::vector<gnSeqI>& breakpoints )
{
	if( mlist.size() == 0 )
		return;
	breakpoints = std::vector<gnSeqI>(1, mlist.size()-1);

	mems::SSC<mems::AbstractMatch> ssc(0);
	std::sort( mlist.begin(), mlist.end(), ssc );
	typedef typename MatchVector::value_type value_type;
	typedef std::pair< value_type, size_t > LabelPairType;
	std::vector< LabelPairType > label_list;
	typename MatchVector::iterator cur = mlist.begin();
	typename MatchVector::iterator end = mlist.end();
	size_t i = 0;
	for( ;cur != end; ++cur )
	{
		label_list.push_back( std::make_pair( *cur, i ) );
		++i;
	}

	uint seq_count = mlist[0]->SeqCount();
	// check for breakpoints in each sequence
	for( uint seqI = 1; seqI < seq_count; seqI++ )
	{
		LabelSort< LabelPairType > ls(seqI); 
		std::sort( label_list.begin(), label_list.end(), ls );

		typename std::vector< LabelPairType >::const_iterator prev = label_list.begin();
		typename std::vector< std::pair< typename MatchVector::value_type, size_t > >::const_iterator iter = label_list.begin();
		typename std::vector< std::pair< typename MatchVector::value_type, size_t > >::const_iterator lab_end = label_list.end();

		bool prev_orient = (*prev).first->Orientation(seqI) == (*prev).first->Orientation(0);
		if( !prev_orient )	// if we start in a different orientation than the ref seq there's a bp here
			breakpoints.push_back(prev->second);

		for( ++iter; iter != lab_end; ++iter )
		{
			bool cur_orient = (*iter).first->Orientation(seqI) == (*iter).first->Orientation(0);
			if( prev_orient == cur_orient &&
				( ( prev_orient && (*prev).second + 1 == (*iter).second) ||
				  ( !prev_orient && (*prev).second - 1 == (*iter).second) 
				)
			  )
			{
				prev_orient = cur_orient;
				++prev;
				continue;	// no breakpoint here
			}

			// always add the last match in a new block (scanning from left to right in seq 0)
			if( prev_orient )
				breakpoints.push_back( prev->second );
			if( !cur_orient )
				breakpoints.push_back( iter->second );

			prev_orient = cur_orient;
			++prev;
		}
		if( prev_orient )
			breakpoints.push_back( prev->second );
	}
	std::sort( breakpoints.begin(), breakpoints.end() );
	std::vector<gnSeqI>::iterator uni = std::unique( breakpoints.begin(), breakpoints.end() );
	breakpoints.erase( uni, breakpoints.end() );
}


template< class MatchVector >
void ComputeLCBs_v2( const MatchVector& meml, const std::vector<gnSeqI>& breakpoints, std::vector< MatchVector >& lcb_list )
{
	// there must be at least one end of a block defined
	if( breakpoints.size() < 1 )
		return;
		
	lcb_list.clear();
	
	// organize the LCBs into different MatchVector instances
	std::vector<gnSeqI>::const_iterator break_iter = breakpoints.begin();
	uint prev_break = 0;	// prev_break is the first match in the current block
	MatchVector lcb;
	for( ; break_iter != breakpoints.end(); ++break_iter ){
		// add the new MatchList to the set if it made the cut
		lcb_list.push_back( lcb );
		lcb_list.back().insert( lcb_list.back().end(), meml.begin() + prev_break, meml.begin() + *break_iter + 1 );
		prev_break = *break_iter + 1;
	}
}


template <class MatchVector>
void computeLCBAdjacencies_v3( const std::vector< MatchVector >& lcb_list, std::vector< double >& weights, std::vector< mems::LCB >& adjacencies )
{
	adjacencies.clear(); // start with no LCB adjacencies
	if( lcb_list.size() == 0 )
		return;	// there aren't any LCBs so there aren't any adjacencies!

	uint seq_count = lcb_list.front().front()->SeqCount();
	uint seqI;
	uint lcbI;
	for( lcbI = 0; lcbI < lcb_list.size(); ++lcbI ){
		mems::LCB lcb;
		std::vector<gnSeqI> left_end;
		std::vector<gnSeqI> length;
		std::vector<bool> orientation;
		FindBoundaries( lcb_list[lcbI], left_end, length, orientation );

		lcb.left_adjacency = std::vector<uint>( left_end.size(), -1 );
		lcb.right_adjacency = std::vector<uint>( left_end.size(), -1 );
		lcb.left_end = std::vector<int64>( left_end.size(), 0 );
		lcb.right_end = std::vector<int64>( left_end.size(), 0 );

		for( seqI = 0; seqI < seq_count; seqI++ ){
			// support "ragged edges" on the ends of LCBs
			if( left_end[seqI] == mems::NO_MATCH )
				continue;
			lcb.left_end[seqI] = left_end[seqI];
			lcb.right_end[seqI] = left_end[seqI] + length[seqI];
			if( !orientation[seqI] )
			{
				lcb.left_end[seqI] = -lcb.left_end[seqI];
				lcb.right_end[seqI] = -lcb.right_end[seqI];
			}
		}
		lcb.lcb_id = adjacencies.size();
		lcb.weight = weights[ lcbI ];
		adjacencies.push_back( lcb );
	}

	for( seqI = 0; seqI < seq_count; seqI++ ){
		mems::LCBLeftComparator llc( seqI );
		std::sort( adjacencies.begin(), adjacencies.end(), llc );
		for( lcbI = 1; lcbI + 1 < lcb_list.size(); lcbI++ ){
			adjacencies[ lcbI ].left_adjacency[ seqI ] = adjacencies[ lcbI - 1 ].lcb_id;
			adjacencies[ lcbI ].right_adjacency[ seqI ] = adjacencies[ lcbI + 1 ].lcb_id;
		}
		if( lcbI == lcb_list.size() )
			lcbI--;	// need to decrement when there is only a single LCB

		// set first and last lcb adjacencies to -1
		adjacencies[ 0 ].left_adjacency[ seqI ] = (uint)-1;
		adjacencies[ lcbI ].right_adjacency[ seqI ] = (uint)-1;
		if( lcbI > 0 ){
			adjacencies[ 0 ].right_adjacency[ seqI ] = adjacencies[ 1 ].lcb_id;
			adjacencies[ lcbI ].left_adjacency[ seqI ] = adjacencies[ lcbI - 1 ].lcb_id;
		}
	}
	mems::LCBIDComparator lic;
	std::sort( adjacencies.begin(), adjacencies.end(), lic );

}

/**
 *  Redesign to be more intuitive.  left_adjacency is always left, regardless of LCB orientation
 */
inline
void computeLCBAdjacencies_v3( mems::IntervalList& iv_list, std::vector< double >& weights, std::vector< mems::LCB >& adjacencies ){
	std::vector< std::vector< mems::Interval* > > nivs;
	for( size_t ivI = 0; ivI < iv_list.size(); ivI++ )
		nivs.push_back( std::vector< mems::Interval* >( 1, &iv_list[ivI] ) );
	computeLCBAdjacencies_v3( nivs, weights, adjacencies );
}

/**
 * Takes a set of filtered LCB adjacencies and an unfiltered set of matches as input
 * returns a filtered set of matches that reflects the LCBs found
 */
template< class MatchVector >
void filterMatches_v2( std::vector< mems::LCB >& adjacencies, std::vector< MatchVector >& lcb_list, std::vector< double >& weights, MatchVector& deleted_matches ){
	if( lcb_list.size() < 1 )
		return;
	MatchVector lcb_tmp = lcb_list[ 0 ];
	lcb_tmp.clear();
	std::vector< MatchVector > filtered_lcbs( lcb_list.size(), lcb_tmp );
	uint lcbI;
	for( lcbI = 0; lcbI < adjacencies.size(); lcbI++ ){
		if( adjacencies[ lcbI ].lcb_id == lcbI ){
			filtered_lcbs[ lcbI ].insert( filtered_lcbs[ lcbI ].end(), lcb_list[ lcbI ].begin(), lcb_list[ lcbI ].end() );
			continue;
		}
		if( adjacencies[ lcbI ].lcb_id == -1 ){
			std::cerr << "weird";
			continue; 	// this one was removed
		}
		if( adjacencies[ lcbI ].lcb_id == -2 )
		{
			deleted_matches.insert( deleted_matches.end(), lcb_list[lcbI].begin(), lcb_list[lcbI].end() );
			continue; 	// this one was removed
		}

		// this one points elsewhere
		// search and update the union/find structure for the target
		std::stack< uint > visited_lcbs;
		visited_lcbs.push( lcbI );
		uint cur_lcb = adjacencies[ lcbI ].lcb_id;
		while( adjacencies[ cur_lcb ].lcb_id != cur_lcb ){
			visited_lcbs.push( cur_lcb );
			cur_lcb = adjacencies[ cur_lcb ].lcb_id;
			if( cur_lcb == -1 || cur_lcb == -2 ){
//				std::cerr << "improper hoodidge\n";
				break;	// this one points to an LCB that got deleted
			}
		}
		while( visited_lcbs.size() > 0 ){
			adjacencies[ visited_lcbs.top() ].lcb_id = cur_lcb;
			visited_lcbs.pop();
		}
		// add this LCB's matches to the target LCB.
		if( cur_lcb != -1 && cur_lcb != -2 )
			filtered_lcbs[ cur_lcb ].insert( filtered_lcbs[ cur_lcb ].end(), lcb_list[ lcbI ].begin(), lcb_list[ lcbI ].end() );
		else
			deleted_matches.insert( deleted_matches.end(), lcb_list[lcbI].begin(), lcb_list[lcbI].end() );
	}


	lcb_list.clear();
	std::vector< double > new_weights;
	for( lcbI = 0; lcbI < filtered_lcbs.size(); lcbI++ ){
		if( filtered_lcbs[ lcbI ].size() > 0 ){
			lcb_list.push_back( filtered_lcbs[ lcbI ] );
			new_weights.push_back( weights[lcbI] );
		}
	}

	// sort the matches inside consolidated LCBs
	mems::MatchStartComparator<mems::AbstractMatch> msc( 0 );
	for( lcbI = 0; lcbI < lcb_list.size(); lcbI++ ){
		std::sort( lcb_list[ lcbI ].begin(), lcb_list[ lcbI ].end(), msc );
	}

	// calculate the LCB adjacencies
	weights = new_weights;
	computeLCBAdjacencies_v3( lcb_list, weights, adjacencies );

}





class EvenFasterSumOfPairsBreakpointScorer
{
public:
	EvenFasterSumOfPairsBreakpointScorer( 
		double breakpoint_penalty,
		boost::multi_array<double,2> bp_weight_matrix, 
		boost::multi_array<double,2> conservation_weight_matrix,
		std::vector< TrackingMatch* > tracking_match,
		mems::PairwiseLCBMatrix& pairwise_adjacency_matrix,
		std::vector<node_id_t>& n1_descendants,
		std::vector<node_id_t>& n2_descendants,
		boost::multi_array< double, 3 >& tm_score_array,
		boost::multi_array< size_t, 3 >& tm_lcb_id_array,
		size_t seqI_begin,
		size_t seqI_end,
		size_t seqJ_begin,
		size_t seqJ_end
		);

	/**
	 * Returns the number of possible moves a search algorithm may make from the current 
	 * location in LCB search space.  In this case it's simply the total number of pairwise LCBs
	 */
	size_t getMoveCount();

	/** returns the score of the current state */
	double score();

	/** scores a move */
	double operator()( std::pair< double, size_t >& the_move  );

	/** checks whether a particular move is a valid move */
	bool isValid( std::pair< double, size_t >& the_move );

	bool remove( std::pair< double, size_t >& the_move, std::vector< std::pair< double, size_t > >& new_move_list, size_t& new_move_count );

	/** applies a score difference */
	void applyScoreDifference( boost::multi_array< double, 2 >& lcb_score_diff, boost::multi_array< size_t, 2 >& lcb_removed_count );

	/** undoes a score difference, if it wasn't accepted for example */
	void undoScoreDifference( boost::multi_array< double, 2 >& lcb_score_diff, boost::multi_array< size_t, 2 >& lcb_removed_count );

	/** returns the maximum number of new moves generated by any LCB removal */
	size_t getMaxNewMoveCount();

	/** call to indicate that the given LCB has been removed 
	  * @param really_remove	set to false if the move should merely be checked for validity
	  * returns false if the move was invalid
	  */
	bool remove( std::pair< double, size_t >& the_move, bool really_remove, 
		boost::multi_array< double, 2 >& lcb_score_diff, boost::multi_array< size_t, 2 >& lcb_removed_count, 
		bool score_new_moves, std::vector< std::pair< double, size_t > >& new_move_list, size_t& new_move_count );

	/** returns the final set of TrackingMatch values which remain after applying greedy breakpoint elimination */
	std::vector< mems::TrackingMatch* > getResults();

	/** sanity checks all internal data structures */
	bool validate();

protected:
	double bp_penalty;
	boost::multi_array<double,2> bp_weights;
	boost::multi_array<double,2> conservation_weights;
	std::vector< mems::TrackingMatch* > tracking_matches;
	mems::PairwiseLCBMatrix pairwise_adjacencies;
	std::vector<node_id_t> n1_des;
	std::vector<node_id_t> n2_des;

	boost::multi_array< size_t, 2 > pairwise_lcb_count;
	boost::multi_array< double, 2 > pairwise_lcb_score;

	std::vector< TrackingMatch* > deleted_tracking_matches;

private:
	// avoid continuous size lookup
	const size_t seqI_count;
	const size_t seqJ_count;

	// variables used during score computation
	boost::multi_array< std::vector< std::pair< uint, uint > >, 2 > all_id_remaps;
	boost::multi_array< std::vector< uint >, 2 > full_impact_list;
	boost::multi_array< double, 2 > internal_lcb_score_diff[3];
	boost::multi_array< size_t, 2 > internal_lcb_removed_count[3];
	int using_lsd;
	std::vector< double > lsd_zeros;
	std::vector< size_t > lrc_zeros;
	std::vector< double > bogus_scores;
	std::vector< size_t > my_del_lcbs;
	std::vector< size_t > lcb_ids;

	boost::multi_array< double, 3 >& tm_score_array;
	boost::multi_array< size_t, 3 >& tm_lcb_id_array;

	// limit to a range of sequences
	const size_t seqI_first;
	const size_t seqJ_first;
	const size_t seqI_last;
	const size_t seqJ_last;
};


template< class BreakpointScorerType >
int64 greedyBreakpointElimination_v4( std::vector< mems::LCB >& adjacencies, std::vector< double >& scores, BreakpointScorerType& bp_scorer, std::ostream* status_out, size_t g1_tag = 0, size_t g2_tag = 0 );

template< class SearchScorer >
double greedySearch( SearchScorer& spbs );


/**
 * A breakpoint scorer that applies a fixed penalty for each breakpoint that exists in a set of
 * two or more sequences 
 */
class SimpleBreakpointScorer
{
public:
	SimpleBreakpointScorer( std::vector< LCB >& adjacencies, double breakpoint_penalty );

	size_t getMoveCount();

	double score();

	bool isValid( size_t lcbI, double move_score );

	/** return the relative change in score if lcbI were to be removed */
	double operator()( size_t lcbI );

	/** call to indicate that the given LCB has been removed */
	void remove( uint lcbI, std::vector< std::pair< double, size_t > >& new_moves );

private:
	std::vector< mems::LCB > adjs;
	double bp_penalty;
	std::vector< double > scores;
	double total_weight;
	size_t bp_count;
};


class GreedyRemovalScorer
{
public:
	GreedyRemovalScorer( std::vector< LCB >& adjacencies, double minimum_weight );

	size_t getMoveCount();

	double score();

	bool isValid( size_t lcbI, double move_score );

	/** return the relative change in score if lcbI were to be removed */
	double operator()( size_t lcbI );

	/** call to indicate that the given LCB has been removed */
	void remove( uint lcbI, std::vector< std::pair< double, size_t > >& new_moves );

private:
	std::vector< mems::LCB > adjs;
	double min_weight;
	std::vector< double > scores;
	double total_weight;
};


}	// namespace mems

#endif // __greedyBreakpointElimination_h__


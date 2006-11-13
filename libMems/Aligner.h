#ifndef _Aligner_h_
#define _Aligner_h_

#include "gn/gnSequence.h"
#include "DNAMemorySML.h"
//#include "DNAFileSML.h"
#include "MemList.h"
#include "MemHash.h"

#ifndef uint
typedef unsigned uint;
#endif

/**
 * A mem labeled with a number.
 * Used by LCB construction algorithm 
 */
class LabeledMem{
public:
	MemHashEntry* mem;
	uint label;
};

/**
 * Compares Matches labeled with a number.
 * Used by LCB construction algorithm 
 */
class LabeledMemComparator {
public:
	LabeledMemComparator( uint seq ){
		msc = new MemStartComparator( seq );
	}
	LabeledMemComparator( const LabeledMemComparator& lmc ){
		*this = lmc;
	}
	LabeledMemComparator& operator=( const LabeledMemComparator& lmc ){
		msc = new MemStartComparator( *(lmc.msc) );
		return *this;
	}
	boolean operator()(const LabeledMem& a, const LabeledMem& b) const{
		return (*msc)( a.mem, b.mem );
	}
	~LabeledMemComparator() {
		delete msc;
	}
protected:
	MemStartComparator* msc;
private:
	LabeledMemComparator();
};


/**
 * Used to find locally colinear blocks (LCBs) and do recursive
 * alignments on the blocks
 * To create an alignment one need only use the align method.
 * LCB lists can be written to and read from a file using the methods
 * included in this class.
 * The other methods are available for experimentation.
 */
class Aligner {
public:
	/** Construct a new aligner with the specified values
	 * @param offset_tolerance A measure of how different a match's offset can be from its neighbors
	 *        before it will be thrown out
	 * @param LCB_minimum_density The minimum density that an LCB may have to be considered a valid block
	 *   						  This should be a number between 0 and 1.
	 * @param LCB_minimum_range   The minimum range an LCB must cover to be considered an LCB.
	 * @throws AlignerError may be thrown if an error occurs
	 */
	Aligner( int64 offset_tolerance, double LCB_minimum_density, double LCB_minimum_range );
	
	/**
	 * Performs an alignment.  Takes a MemList as input and outputs a list of LCBs.
	 * @param mlist The list of Mems in which to look for rearrangements.
	 * @param LCB_list The resulting list of LCBs will be placed in this vector.
	 * @param min_match_size The minimum size of matches that should be allowed when determining rearrangements
	 */
	void align( MemList& mlist, vector<MemList>& LCB_list, uint min_match_size = 0 );

	/**
	 * Bob's Cheesy LCB algorithm.  
	 * This vehicle must be driven by a domain expert
	 * This algorithm is a two step process:
	 * First out of place MUMs are removed.
	 * This is accomplished by sorting the MUMs on the first sequence
	 * and tossing out any MUM whose generalized offset is more than "offset_threshold"
	 * away from its neighbors.
	 * Each MUM which passes the filter is labeled with a number starting at 0
	 * and increasing by one for each MUM.
	 * The second step involves sorting the labeled MUMs on their start positions in
	 * the other sequences and looking for breakpoints.  Breakpoints are encountered
	 * when a label changes by more than 1 between neighboring MUMs.  
	 * The label of the last MUM in the ongoing block is recorded in the breakpoint set
	 * After searching for breakpoints in each sequence the set of breakpoints is
	 * returned.
	 */
	void BobsCheesyLCB( MemList& mlist, int64 offset_threshold, set<uint>& breakpoints );
	/**
	 * Extracts a set of LCBs from a MemList using a set of breakpoints as input.
	 * @param meml The MemList containing a set of matches
	 * @param breakpoints A set of integers corresponding to breakpoints in the mem list
	 * @param lcb_list An empty LCB list that will be filled with LCBs when this function is called
	 * @param LCB_minimum_density The minimum percent identity and LCB must have
	 * @param LCB_minimum_range The minimum range in base pairs that an LCB must cover
	 * @param LCB_minimum_matches The minimum number of matches that an LCB must contain
	 */
	void FilterLCBs( MemList& meml, set<uint>& breakpoints, vector<MemList>& lcb_list, double LCB_minimum_density, double LCB_minimum_range, uint LCB_minimum_matches );

	void writeLCBlist( const vector<MemList>& LCB_list, ostream& lcb_out );

protected:
	vector<gnSequence*> sequences;
	vector<SortedMerList*> suffixArrays;
	MemList ml;
	uint32 seq_count;
	boolean debug;
	
	int64 offset_tolerance;
	double LCB_minimum_density;
	double LCB_minimum_range;
};

/**
 * Thrown if some error occurs during alignment
 */
CREATE_EXCEPTION( AlignerError );

#endif // _Aligner_h_

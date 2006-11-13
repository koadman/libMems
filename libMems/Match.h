#ifndef _Match_h_
#define _Match_h_

#include "gnClone.h"
#include <iostream>
#include <vector>
#include <set>
#include <valarray>

static const gnSeqI MEM_NO_MATCH = 0;

typedef unsigned MatchID_t;

/**
 * The Match class stores the location of a matching region of equal size
 * between several sequences.  There are numerous functions in this
 * class which can be used to compare and manipulate this match.
 */

class Match : public gnClone {

public:
	Match();
	/**
	 * Creates a new Match.
	 * @param seq_count The total number of sequences in the alignment
	 * @param mersize The size of the mers used in the sorted mer lists.
	 * @param m_type The type of mem to create, can either be a seed or already extended.
	 * @see MatchType
	 */
	Match(const uint32 seq_count );
	Match(const uint32 seq_count, const MatchID_t match_id );
	Match(const Match& mhe);
	Match& operator=(const Match& mhe);
	virtual Match* Clone() const;
	~Match();
	
	/** comparison operator, compares two matches to see if they are the same */
	boolean operator==(const Match& mhe) const;
	
	MatchID_t MatchID() const { return match_id; }

	/** Allocates a new unique match identifier for this Match */
	void UseNewMatchID(){ match_id = GetNewMatchID(); }

	/** Returns the generalized offset of this match */
	virtual int64 Offset() const{return m_offset;}

	/** Sets the generalized offset of this match to "offset" */
	virtual void SetOffset(int64 offset){m_offset = offset;}

	/** Returns the specified sequence's offset in this match */
	virtual int64 Offset( uint32 seqI ) const{ return m_offsets[ seqI ]; }

	/** Returns the length of this match */
	virtual gnSeqI Length() const{return m_length;}

	/** Sets the length of this match to @param len */
	virtual void SetLength(gnSeqI len){m_length = len;}

	/** Returns the start coordinate of this match in sequence @param startI */
	virtual int64 Start(uint32 startI) const{return m_start[startI];}

	/** Sets the start in sequence @param seqI of this match to @param startI */
	virtual void SetStart(uint32 seqI, int64 startI);

	/** Returns the start coordinate of this match in sequence @param startI */
	virtual int64 operator[](uint32 startI) const{return m_start[startI];}

	/** Returns the start coordinate of this match in sequence @param startI */
	virtual int64& operator[](uint32 startI) {return m_start[startI];}

	/** Returns the last coordinate of this match in sequence @param endI */
	virtual int64 End(uint32 endI) const;

	/** Returns the multiplicity of the match.  e.g. the number of sequences this match occurs in */
	virtual uint32 Multiplicity() const{return m_multiplicity;}

	/** Returns the number of sequences in the alignment which contains this match */
	virtual uint32 SeqCount() const{return m_start.size();}

	/** Returns the index of the first sequence this match occurs in */
	virtual uint32 FirstStart() const;

	/** Returns a numerical representation of which sequences are matched */
	virtual uint64 MatchNumber() const {return m_matchnumber;}

//	/** Returns a numerical representation of which sequences are matched */
//	virtual valarray<bool> MatchNumber() const {return m_matchnumber;}
	
	/**
	 * Links a subset match to this mem and this mem to the subset.
	 * @param subset The subset match which this mem overlaps
	 */
	virtual void LinkSubset( Match* subset );

	/** Returns a copy of set of mems which overlap this mem */
	const set<MatchID_t> Subsets() const { return subsets; }

	/** Returns a copy of set of mems which this mem overlaps */
	const set<MatchID_t> Supersets() const { return supersets; }
	
	/** Returns a copy of set of mems which overlap this mem */
	set<MatchID_t>& Subsets() { return subsets; }

	/** Returns a copy of set of mems which this mem overlaps */
	set<MatchID_t>& Supersets() { return supersets; }

	/**
	 * Inverts the coordinates of this mem.
	 */
	virtual void Invert();
	
	/**
	 * Calculates the generalized offset and other bookkeeping information
	 * for this mem.  This should <b>always</b> be called after changing the start
	 * positions of the mem.
	 */
	virtual void CalculateOffset();
	
	/**
	 * Return the difference between the end of this mem and the start
	 * of mhe in the sequence specified by seqI
	 * @param mhe The mem to compare to.
	 * @param seqI The index of the sequence to compare
	 * @return The gap size.
	 */
	virtual int64 GapSize(const Match& mhe, uint32 seqI);
	
	//warning:  none of the following do bounds checking.
	virtual void Move( int64 distance );
	virtual void CropStart(gnSeqI crop_amount);
	virtual void CropEnd(gnSeqI crop_amount);
	virtual void ExtendStart(gnSeqI extend_amount);
	virtual void ExtendEnd(gnSeqI extend_amount);
	//comparison functions

	/**
	 *	Will return true if this mem contains mhe
	 *  Containment implies that a mem has a length >= the contained
	 *  mem, it has coordinates in every genome the contained mem has,
	 *  the difference in start positions in each genome is the same.
	 * @param mhe The mem to check for containment.
	 * @return True if this mem contains mhe.
	 */
	virtual boolean Contains(const Match& mhe) const;

	/**
	 *  Will return true if this mem intersects mhe
	 *  Intersection implies that in sequences which both matches have 
	 *  coordinates the difference in start positions must always be the
	 *  same and the range of the coordinates overlap.  Each match may have
	 *  coordinates in sequences which the other does not.
	 * @param mhe The Match to check for intersection.
	 * @return True if this Match intersects with mhe.
	 */
	virtual boolean Intersects(const Match& mhe) const;
	virtual boolean GetUniqueStart(const Match& mhe, Match& unique_mhe) const;
	virtual boolean GetUniqueEnd(const Match& mhe, Match& unique_mhe) const;
	
	/**
	 *  Will return true if this mem evenly overlaps mhe
	 *  An even overlap implies that for each genome which this mem has 
	 *  coordinates, the difference in start positions must always be the
	 *  same and this mem's coordinates completely contain mhe's in those
	 *  genomes.
	 *  mhe may have coordinates in genomes which this mem does not.
	 * @param mhe The mem to check for intersection.
	 * @return True if this mem intersects with mhe.
	 */
	virtual boolean EvenlyOverlaps(const Match& mhe) const;
	
	/**
	 * Writes the location of this match to the specified output stream (e.g. cout).
	 */
	friend std::ostream& operator<<(std::ostream& os, const Match& mhe); //write to source.

	virtual void UnlinkSubset( MatchID_t sub_id );
	virtual void UnlinkSuperset( MatchID_t super_id );

	virtual void AddSubset( MatchID_t underlapper );
	virtual void AddSuperset( MatchID_t overlapper );

protected:
	static MatchID_t GetNewMatchID();
	void HomogenizeParity(vector<int64>& start1, vector<int64>& start2, const uint32 startI) const;
	set<MatchID_t> subsets;
	set<MatchID_t> supersets;
	
	MatchID_t match_id;
	int64 m_offset;
	gnSeqI m_length;
	vector<int64> m_start;
	uint32 m_multiplicity;
	uint32 m_firstStart;
//	valarray<bool> m_matchnumber;
	uint64 m_matchnumber;
	vector<int64> m_offsets;
};


inline
void Match::SetStart(uint32 seqI, int64 startI){
	m_start[seqI] = startI;
}

std::ostream& operator<<(std::ostream& os, const Match& mhe); //write to source.
/*
class MheCompare {
public:
	boolean operator()(const Match* a, const Match* b) const{
		int64 a_matchnum = a->MatchNumber();
		int64 b_matchnum = b->MatchNumber();
		if( a_matchnum < b_matchnum ){
			return true;
		}else if( a_matchnum == b_matchnum ){
			//offsets are the same, check for containment...
			if(a->Contains(*b) || b->Contains(*a)){
				return false;
			}else
				return Match::strict_start_lessthan_ptr(a, b);
		}
		return false;
	}
};
*/
class MatchStartComparator {
public:
	MatchStartComparator( unsigned seq = 0 ){
		m_seq = seq;
	}
	MatchStartComparator( MatchStartComparator& msc ){
		m_seq = msc.m_seq;
	}
	// TODO??  make this do a wraparound comparison if all is equal?
	boolean operator()(const Match* a, const Match* b) const{
		int32 start_diff = max( a->FirstStart(), m_seq ) - max( b->FirstStart(), m_seq );
		if(start_diff == 0){
			uint32 m_count = a->SeqCount();
			m_count = m_count <= b->SeqCount() ? m_count : b->SeqCount();
			for(uint32 seqI = m_seq; seqI < m_count; seqI++){
				int64 a_start = a->Start( seqI ), b_start = b->Start( seqI );
				if(a_start < 0)
					a_start = -a_start + a->Length();
				if(b_start < 0)
					b_start = -b_start + b->Length();
				int64 diff = a_start - b_start;
				if(a_start == MEM_NO_MATCH || b_start == MEM_NO_MATCH)
					continue;
				else if(diff == 0)
					continue;
				else
					return diff < 0;
			}
		}
		return start_diff < 0;
	}
private:
	unsigned m_seq;
};

#endif // _Match_h_

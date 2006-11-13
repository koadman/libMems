#include "MatchFinder.h"
#include <list>
#include "MemHashEntry.h"
#include <iostream>

static const uint32 DEFAULT_MEM_TABLE_SIZE = 40000;
static const uint32 DEFAULT_REPEAT_TOLERANCE = 4;
static const uint32 DEFAULT_ENUMERATION_TOLERANCE = 1;

/**	The mem struct represents a matching section between multiple
 *  genomes.  The start array contains the start position of the
 *  match within each genome, or MEM_NO_MATCH if no match exists
 *  in that particular genome.
 */
struct mem{
	gnSeqI length;		///the length of the mem
	int64* start;		///the array of start positions
};

class MemHash : public MatchFinder{
public:
	MemHash(uint32 tablesize = DEFAULT_MEM_TABLE_SIZE);
	~MemHash();
	MemHash(const MemHash& mh);
	virtual MemHash* Clone() const;
	virtual boolean Create(uint32 mer_mask_size);
	virtual boolean GetMemList( vector<MemHashEntry>& mem_list );
	virtual void EliminateInclusions();
	virtual uint32 MemCount(){return m_mem_count;}
	virtual uint32 MemCollisionCount(){return m_collision_count;}
	virtual void MemTableCount(vector<uint32>& table_count){table_count = mem_table_count;}
	virtual void PrintDistribution(ostream& os);

	virtual void SetRepeatTolerance(uint32 repeat_tolerance){m_repeat_tolerance = repeat_tolerance;}
	virtual uint32 GetRepeatTolerance(){return m_repeat_tolerance;}
	virtual void SetEnumerationTolerance(uint32 enumeration_tolerance){m_enumeration_tolerance = enumeration_tolerance;}
	virtual uint32 GetEnumerationTolerance(){return m_enumeration_tolerance;}
protected:
	virtual boolean EnumerateMatches( list<idmer>& match_list );
	virtual boolean EnumerateRepeats( list<idmer>& match_list );
	virtual boolean HashMatch(list<idmer>& match_list);
	virtual void SetDirection(MemHashEntry& mhe);
	virtual boolean AddHashEntry(MemHashEntry& mhe);
	virtual uint32 quadratic_li(uint32 listI){return (listI*(listI+1))/2;}

	uint32 table_size;
	list<MemHashEntry>* mem_table;
	uint32 m_repeat_tolerance;
	uint32 m_enumeration_tolerance;
	uint64 m_mem_count;
	uint64 m_collision_count;
	vector<uint32> mem_table_count;
};

inline
void MemHash::PrintDistribution(ostream& os){
	for(uint32 i=0; i < mem_table_count.size(); i++)
		os << i << '\t' << mem_table_count[i] << '\n';
}

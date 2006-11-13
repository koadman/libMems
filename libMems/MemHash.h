#ifndef _MemHash_h_
#define _MemHash_h_

#include "MatchFinder.h"
#include <set>
#include "MemHashEntry.h"
#include "MemList.h"
#include <iostream>
#include "gn/gnException.h"
#include "MatchList.h"

static const uint32 DEFAULT_MEM_TABLE_SIZE = 40000;
static const uint32 DEFAULT_REPEAT_TOLERANCE = 0;
static const uint32 DEFAULT_ENUMERATION_TOLERANCE = 0;

/**
 * MemHash implements an algorithm for finding exact matches of a certain minimal
 * length in several sequences.
 */
class MemHash : public MatchFinder{
public:
	MemHash();
	~MemHash();
	MemHash(const MemHash& mh);
	virtual MemHash* Clone() const;
	virtual void Clear();
	
	/**
	 * Finds all maximal exact matches in the sequences contained by "match_list"
	 * The resulting list of matches is stored within "match_list"
	 */
	virtual void FindMatches( MatchList& match_list );

	/**
	 * Finds all maximal exact matches in the sequences contained by "mem_list"
	 * The resulting list of mems is stored within "mem_list"
	 */
	virtual void FindMems( MemList& mem_list );

	/**
	 * Generates exact matches for the sequences loaded into this MemHash 
	 */
	virtual boolean CreateMatches();

	/**
	 * Returns the size of the hash table being used. 
	 * @return the size of the hash table being used. 
	 */
	virtual uint32 TableSize() const {return table_size;};
	/**
	 * Sets the size of the hash table to new_table_size.
	 * @param new_table_size The new hash table size
	 */
	virtual void SetTableSize(uint32 new_table_size);
	/**
	 * Creates a new MatchList instance which contains all the matches found by calling Create().
	 */
	virtual MatchList GetMatchList() const;
	/**
	 * Places pointers to the mems that have been found into the vector mem_list
	 * @param mem_list an empty vector.
	 */
	virtual void GetMemList( vector<MemHashEntry*>& mem_list ) const;
	/**
	 * Creates a new MemList instance which contains all the mems found by calling Create().
	 */
	virtual MemList GetMemList() const;
	
	/**
	 * Returns the number of mems found 
	 * @return The number of mems found 
	 */
	virtual uint32 MemCount(){return m_mem_count;}
	/**
	 * Returns the number of mers thrown out because they were contained in an existing mem 
	 * @return The number of mers thrown out because they were contained in an existing mem 
	 */
	virtual uint32 MemCollisionCount(){return m_collision_count;}
	virtual void MemTableCount(vector<uint32>& table_count){table_count = mem_table_count;}
	/**
	 * Prints the number of matches in each hash table bucket to the ostream os.
	 * @param os The stream to print to.
	 */
	virtual void PrintDistribution(ostream& os) const;
	
	/**
	 * Reads in a list of mems from an input stream
	 * @throws A InvalidFileFormat exception if the file format is unknown or the file is corrupt
	 */
	virtual void LoadFile(istream& mem_file);
	/**
	 * Writes the matches stored in this MemHash out to the ostream @param mem_file.
	 */
	virtual void WriteFile(ostream& mem_file) const;

	/**
	 * Sets the permitted repetitivity of match seeds.  
	 * Set @param repeat_tolerance to 0 to generate MUMs, any higher setting will generate MEMs
	 * Many possible combinations of repetitive seed matches may be ignored, depending on the 
	 * setting of the repeat enumeration tolerance.
	 * @see SetEnumerationTolerance
	 * @param repeat_tolerance the permitted repetitivity of match seeds
	 */
	virtual void SetRepeatTolerance(uint32 repeat_tolerance){m_repeat_tolerance = repeat_tolerance;}
	/**
	 * @return the permitted repetitivity of match seeds.  
	 * @see SetRepeatTolerance
	 */
	virtual uint32 GetRepeatTolerance() const{return m_repeat_tolerance;}
	/**
	 * Sets the match seed repeat enumeration tolerance.
	 * When matching mers are found across sequences which also occur several times in any particular
	 * sequence there are several possible match seeds which could be generated.
	 * The enumeration tolerance controls how many of these possibilities are actually used as match
	 * seeds and extended into full matches.  The selection of actual seeds from the realm of possibilities
	 * is essentially arbitrary, though not explicitly randomized.
	 */
	virtual void SetEnumerationTolerance(uint32 enumeration_tolerance){m_enumeration_tolerance = enumeration_tolerance;}
	/**
	 * @return  the match seed repeat enumeration tolerance.
	 * @see SetEnumerationTolerance
	 */
	virtual uint32 GetEnumerationTolerance() const{return m_enumeration_tolerance;}

	/**
	 * Not implemented 
	 */
	virtual void BreakOnAmbiguities();

protected:
	virtual boolean EnumerateMatches( list<idmer>& match_list );
	virtual boolean HashMatch(list<idmer>& match_list);
	virtual void SetDirection(MemHashEntry& mhe);
	virtual MemHashEntry* AddHashEntry(MemHashEntry& mhe);
	virtual uint32 quadratic_li(uint32 listI){return (listI*(listI+1))/2;}
		
	uint32 table_size;
	set<MemHashEntry*, MheCompare>* mem_table;
	uint32 m_repeat_tolerance;
	uint32 m_enumeration_tolerance;
	uint64 m_mem_count;
	uint64 m_collision_count;
	vector<uint32> mem_table_count;
	
};


/**
 * Thrown when a file being read is invalid
 */
CREATE_EXCEPTION(InvalidFileFormat)


#endif //_MemHash_h_

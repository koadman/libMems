#ifndef _SortedMerList_h_
#define _SortedMerList_h_

#include "gnDefs.h"
#include "gnClone.h"
#include "gnDebug.h"
#include "gnSequence.h"
#include "gnException.h"
#include "stdlib.h"
#include <string>
#include <vector>

#define DNA_ALPHA_BITS 2	/**< number of bits to represent each nucleotide of DNA */
#define DNA_MER_SIZE 31		/**< largest possible number of characters in each dna mer ALWAYS ODD */

#define PROTEIN_ALPHA_BITS 5	/**< number of bits to represent each amino acid */
#define PROTEIN_MER_SIZE 12		/**< default number of characters in each protein mer */

#define DESCRIPTION_SIZE 2048	/**< Number of bytes for the freeform text description of an SML */

#define NO_UNIQUE_COUNT UINT32_MAX

typedef int8 sarID_t;

//4 + 8 = 16 (blame C alignment rules.)
struct bmer{
	gnSeqI position;	/**< starting position of this mer in the sequence */
	uint64 mer; 		/**< the actual binary encoded mer */
};

//268 bytes 0x10C
struct SMLHeader{
	uint32 version;						/**< Format version - 4 bytes */
	uint32 alphabet_bits;				/**< Bits per character in the alphabet - 4 bytes */
	uint32 mer_size;					/**< Size of mers used for sorting the list - 4 bytes */
	uint64 length;						/**< length of the sequence before circularity - 8 bytes */
	uint32 unique_mers;					/**< Number of unique mers in the sequence 4 bytes */
	uint32 word_size;					/**< Word size on the machine the sequence was translated */
	boolean little_endian;				/**< Is the byte order little endian?  0==no, !0==yes */
	sarID_t id;							/**< Obsolete ID value - 1 byte, eaten by alignment? */
	boolean circular;					/**< Circularity of sequence - 1 byte */
	uint8 translation_table[UINT8_MAX];	/**< Translation table for ascii characters to binary values -- 256 bytes */
	char description[DESCRIPTION_SIZE]; /**< Freeform text description of sequence data -- 2048 bytes */
};

/**
 * A base class which defines an interface common to all sorted mer lists
 */
class SortedMerList : public gnClone
{
public:
	SortedMerList();
	SortedMerList( const SortedMerList& sa );
	SortedMerList& operator=(const SortedMerList& sa);
	~SortedMerList();

	/**
	 * Creates a new sorted mer list.
	 * This function enumerates each possible mer of the specified size and 
	 * sorts them alphabetically in order to construct a sorted mer list.
	 * @param seq The sequence to create an SML for.
	 * @param mersize The size of the mers to sort on.
	 */
	virtual void Create(const gnSequence& seq, const uint32 mersize = DNA_MER_SIZE);
	/**
	 * Read a range of mers in the sorted mer list.
	 * This function reads a section of data from the sorted mer list starting at 'offset'
	 * and continuing for 'size' mers.  The mers are placed into readVector.  Anything
	 * already in readVector is cleared.  Returns false if there was a problem completing the
	 * read.  If the end of the list is reached, all mers which could be read will be placed
	 * into readVector and false will be returned
	 * @param readVector the vector to read bmers into.
	 * @param size The number of bmers to read.
	 * @param offset The mer index in the sorted mer list to start reading from. 
	 * @return false if a problem was encountered while reading.
	 */
	virtual boolean Read(vector<bmer>& readVector, gnSeqI size, const gnSeqI offset) = 0;
	/**
	 * Merges two SortedMerLists.
	 */
	virtual void Merge(SortedMerList& sa, SortedMerList& sa2) = 0;
	
	/**
	 * Get the mer at the specified index in the sorted mer list.
	 * @param index The index of the mer to return.
	 * @return The specified mer.
	 */
	virtual bmer operator[](gnSeqI index) = 0;
	/**
	 * Get the mer at the specified index in the sorted mer list.
	 * @param position The index of the mer to return.
	 * @return The specified mer.
	 */
	virtual uint64 GetMer(gnSeqI position);
	/**
	 * Searches the SML for a subsequence which matches the query string.
	 * Returns true if one is found, false otherwise.
	 * If no matching mer is found, 'result' contains the index that the query
	 * sequence would be in if it existed in the SML.
	 */
	virtual boolean Find(const string& query_seq, gnSeqI& result);
	/**
	 * Searches the SML for a mer which matches the query mer.
	 * Returns true if one is found, false otherwise.
	 * If no matching mer is found, 'result' contains the index that the query
	 * mer would be in if it existed in the SML.
	 */
	virtual boolean FindMer(const uint64 query_mer, gnSeqI& result);
	/**
	 * Searches the SML for mers which match the query mer.
	 * Puts the indices of all matching mers into the 'result' vector
	 */
	virtual void FindAll(const string& query_seq, vector<gnSeqI> result);
	/**
	 * Returns the number of unique mers in the sequence
	 */
	virtual gnSeqI UniqueMerCount();
	
	/**
	 * Returns a freeform text description of the SML.
	 */
	virtual string Description() const;
	/**
	 * Sets the freeform text description of the SML.
	 */
	virtual void SetDescription(const string& d);
	/**
	 * Returns the size of the mers that this SML was sorted on.
	 */
	virtual uint32 MerSize() const;
	/**
	 * Returns the length of the mer mask.
	 * Some types of sorted mer list support a configurable mer mask size, allowing
	 * the same sorted mer list to behave as though it were sorted on a shorter mer size.
	 * DNA sorted mer lists do not support this feature.
	 */
	virtual uint32 GetMerMaskSize() const;
	/**
	 * Sets the length of the mer mask.
	 * Some types of sorted mer list support a configurable mer mask size, allowing
	 * the same sorted mer list to behave as though it were sorted on a shorter mer size.
	 * DNA sorted mer lists do not support this feature.
	 */
	virtual void SetMerMaskSize(uint32 mer_size);
	/**
	 * Returns the length of the sorted mer list.
	 */
	gnSeqI Length() const;
	/**
	 * Ignore this.
	 */
	virtual sarID_t GetID() const;
	/**
	 * Ignore this.
	 */
	virtual void SetID(const sarID_t d);
	/**
	 * Returns true if this SML is circular.  False otherwise.
	 */
	virtual boolean IsCircular() const;
	/**
	 * Returns a mask which can be bitwise AND'ed to a mer in order to
	 * get only the relevant bits of sequence data without direction bits.
	 */
	virtual uint64 GetMerMask() const;
	/**
	 * Returns a copy of the header information for this SML.
	 */
	virtual SMLHeader GetHeader() const;
	/**
	 * Returns a translation table for DNA sequence which disambiguates each nucleotide.
	 */
	static const uint8* const BasicDNATable();
	/**
	 * Returns a translation table for Protein sequence.
	 */
	static const uint8* const ProteinTable();
	/** 
	 * Places a copy of the binary encoded sequence data into dest.
	 * @param len The length in sequence characters to copy
	 * @param offset The sequence offset to start copying from
	 * @throws IndexOutOfBounds if offset or len are invalid
	 */
	virtual void GetBSequence(uint32* dest, const gnSeqI len, const gnSeqI offset);
protected:
	struct SMLHeader header; /**< stores general information about this sorted mer list */
	uint64 mer_mask;	/**< a mask for the used bits in a mer */
	uint32 mask_size;   /**< the number of characters covered by the mask */
	uint32 *sequence;	/**< Stores the sequence data */
	gnSeqI binary_seq_len;	/**< Stores the length in 32 bit words of the sequence */

	/** Set the sequence data to the seq_len characters in seq_buf */
	virtual void SetSequence(gnSeqC* seq_buf, gnSeqI seq_len);
	/** Fill in the vector of bmers with the initial unsorted bmers for the sequence in seq_buf  */
	virtual void FillSML(gnSeqC* seq_buf, gnSeqI seq_len, boolean circular, vector<bmer>& sml_array);
	virtual void FillSML(const gnSequence& seq, vector<bmer>& sml_array);
	virtual void FillDnaSML(const gnSequence& seq, vector<bmer>& sml_array);
	/** Fill in the vector of positions with the initial unsorted positions for the sequence in seq_buf  */
	virtual void FillSML(gnSeqI seq_len, vector<gnSeqI>& sml_array);
	virtual uint64 GetDnaMer(gnSeqI offset);

	virtual gnSeqI bsearch(const struct bmer& query_mer, const gnSeqI start, const gnSeqI end);
	virtual void translate(uint8* dest, const gnSeqC* src, const gnSeqI len) const;
	virtual void translate32(uint32* dest, const gnSeqC* src, const gnSeqI len) const;
	/**
	 * Shifts an entire array of words left or right by a few bits
	 * @param data A pointer to the array of words
	 * @param bits The number of bits to shift by.  A positive number shifts right and a negative number shifts left.
	 */
	virtual void ShiftWords(uint32* data, uint32 length, int32 bits);
	virtual uint32 CalculateMaxMerSize() const;

	static const uint8* const CreateBasicDNATable();
	static const uint8* const CreateProteinTable();
};

/**
 * Thrown when there is an error creating a sorted mer list.
 */
CREATE_EXCEPTION(SMLCreateError);

/**
 * Thrown when there is an error merging two sorted mer lists.
 */
CREATE_EXCEPTION(SMLMergeError);

class MerCompare {
public:
	MerCompare( SortedMerList* sa ){ sar = sa; }
	boolean operator()(const gnSeqI a, const gnSeqI b) const{
		return sar->GetMer(a) < sar->GetMer(b);
	}
protected:
	SortedMerList* sar;
};

bool bmer_lessthan(const bmer& a_v, const bmer& m_v);
bool bmer_id_lessthan(const bmer& a_v, const bmer& m_v);

int bmer_compare(const void* a_v, const void* m_v);
bool bmer_id_lessthan(const bmer& a_v, const bmer& m_v);

//less than function for STL sort functions
inline
bool bmer_lessthan(const bmer& a_v, const bmer& m_v){
	return (a_v.mer < m_v.mer);// ? true : false;
};

inline
int bmer_compare(const void* a_v, const void* m_v){
	return (int)((int64)(((bmer*)a_v)->mer) - (int64)(((bmer*)m_v)->mer));
}


#endif   //_SortedMerList_h_

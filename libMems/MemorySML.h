#ifndef _MemorySML_h_
#define _MemorySML_h_

#include "gnSequence.h"
#include "SortedMerList.h"


/** The MemorySML is an implementation of sorted mer lists which creates and
 *  stores the sorted mer list entirely in memory.  A MemorySML consumes
 *  roughly 32 + alpha_bits bits of memory per character in the sequences.
 *  For unambiguous DNA sequences 4.25 bytes per base are required.
 */
class MemorySML : public SortedMerList
{
public:
	/** 
	 *  Create an empty MemorySML
	 *  Creates an empty MemorySML with the supplied translation
	 *  table and alphabet bit size.  Defaults to DNA settings
	 *  @param table The array used to translate characters into binary code
	 *  @param alpha_bits The number of bits each character consumes in binary
	 */
	MemorySML(const uint8* table = SortedMerList::BasicDNATable(), const uint32 alpha_bits = DNA_ALPHA_BITS);
	MemorySML(const MemorySML& msa);
	MemorySML& operator=(const MemorySML& msa );
	MemorySML* Clone() const;
	
	~MemorySML();

	virtual void Create(const gnSequence& seq, const uint32 mersize = DNA_MER_SIZE);
	virtual boolean Read(vector<bmer>& readVector, gnSeqI size, gnSeqI offset = 0);
	virtual void Merge(SortedMerList& sa, SortedMerList& sa2);
	
	virtual bmer operator[](gnSeqI index);
	
protected:

//	virtual void FillSML(const gnSeqI seq_len, vector<gnSeqI>& sml_array);
	vector<gnSeqI> positions;

};

#endif   //_MemorySML_h_

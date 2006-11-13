#ifndef _DNAMemorySML_h_
#define _DNAMemorySML_h_

#include "gnSequence.h"
#include "MemorySML.h"


/** The DNAMemorySML is an implementation of sorted mer lists which creates and
 *  stores the sorted mer list entirely in memory.  A DNAMemorySML consumes
 *  roughly 32 + alpha_bits bits of memory per character in the sequences.
 *  For unambiguous DNA sequences 4.25 bytes per base are required.
 */
class DNAMemorySML : public MemorySML
{
public:
	/** 
	 *  Create an empty DNAMemorySML
	 *  Creates an empty DNAMemorySML with the supplied translation
	 *  table and alphabet bit size.  Defaults to DNA settings
	 *  @param table The array used to translate characters into binary code
	 *  @param alpha_bits The number of bits each character consumes in binary
	 */
	DNAMemorySML(const uint8* table = SortedMerList::BasicDNATable(), const uint32 alpha_bits = DNA_ALPHA_BITS);
	DNAMemorySML(const DNAMemorySML& msa);
	DNAMemorySML(const SortedMerList& sa);
	DNAMemorySML& operator=(const DNAMemorySML& msa );
	DNAMemorySML* Clone() const;
	
	
	virtual uint64 GetMer(gnSeqI offset);
	
protected:

	virtual void FillSML(const gnSequence& seq, vector<bmer>& sml_array);

};

#endif   //_DNAMemorySML_h_

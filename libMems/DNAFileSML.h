#ifndef _DNAFileSML_h_
#define _DNAFileSML_h_

#include "FileSML.h"

class DNAFileSML : public FileSML
{
public:
	DNAFileSML();
	
	/** 
	 *  Load or create a DNAFileSML ()
	 *  Attempts to load a DNA sorted mer list from the named file if it exists.
	 *  If the given file does not exist it creates an empty DNAFileSML with 
	 *  the supplied translation table and alphabet bit size.
	 *  @param fname The name of the file to create.
	 *  @param table The array used to translate characters into binary code
	 *  @param alpha_bits The number of bits each character consumes in binary
	 */
	DNAFileSML(const string& fname, const uint8* table = SortedMerList::BasicDNATable(), const uint32 alpha_bits = DNA_ALPHA_BITS);
	DNAFileSML(const SortedMerList& sa);
	DNAFileSML& operator=(const DNAFileSML& msa );
	~DNAFileSML();
	
	DNAFileSML* Clone() const;
	
	virtual uint64 GetMer(gnSeqI position);
	
	virtual uint32 FormatVersion();

protected:
	virtual void FillSML(const gnSequence& seq, vector<bmer>& sml_array);
	virtual uint32 CalculateMaxMerSize() const;
	virtual uint64 GetNeededMemory(gnSeqI len);
};

inline
uint32 DNAFileSML::FormatVersion(){
	static uint32 f_version = 3;
	return f_version;
}

#endif   //_DNAFileSML_h_

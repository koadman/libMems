// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
    #pragma hdrstop
#endif

// for all others, include the necessary headers (this file is usually all you
// need because it includes almost all "standard" wxWindows headers)
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include "gnSequence.h"
#include "DNAMemorySML.h"


DNAMemorySML::DNAMemorySML(const uint8* table, const uint32 alpha_bits) : 
MemorySML( table, alpha_bits )
{}

DNAMemorySML& DNAMemorySML::operator=(const DNAMemorySML& msa ){
	MemorySML::operator=(msa);
	return *this;
}

DNAMemorySML* DNAMemorySML::Clone() const{
	DNAMemorySML *bdsa = new DNAMemorySML();
	(*bdsa) = *this;
	return bdsa;
}

uint64 DNAMemorySML::GetMer(gnSeqI position){
	return GetDnaMer( position );
}

void DNAMemorySML::FillSML(const gnSequence& seq, vector<bmer>& sml_array)
{
	FillDnaSML(seq, sml_array);
}

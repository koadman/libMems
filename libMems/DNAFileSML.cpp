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
#include "gnFilter.h"
#include "DNAFileSML.h"


DNAFileSML::DNAFileSML() : FileSML(){
	FileSML::header.version = FormatVersion();
}

DNAFileSML::~DNAFileSML(){
}

DNAFileSML::DNAFileSML(const string& fname, const uint8* table, const uint32 alpha_bits){
	header.alphabet_bits = alpha_bits;
	memcpy(header.translation_table, table, UINT8_MAX);
	filename = fname;
	header.version = FormatVersion();
}

DNAFileSML& DNAFileSML::operator=(const DNAFileSML& msa ){
	FileSML::operator=(msa);
	return *this;
}

DNAFileSML* DNAFileSML::Clone() const{
	DNAFileSML *bdsa = new DNAFileSML();
	(*bdsa) = *this;
	return bdsa;
}

uint64 DNAFileSML::GetNeededMemory(gnSeqI len){
	uint64 neededmem = (len * FileSML::header.alphabet_bits) / 8;
	//forward and reverse copies of the sequence
	neededmem += len * 2;
	neededmem += sizeof(bmer) * len;
	return neededmem;
}

uint32 DNAFileSML::CalculateMaxMerSize() const{
	return 62 / header.alphabet_bits;
}

uint64 DNAFileSML::GetMer(gnSeqI position){
	return GetDnaMer( position );
}

void DNAFileSML::FillSML(const gnSequence& seq, vector<bmer>& sml_array)
{
	FillDnaSML(seq, sml_array);
}

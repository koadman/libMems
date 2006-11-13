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
#include "BigDnaSar.h"

uint32 BigDnaSar::FORMAT_VERSION = 2;

BigDnaSar::BigDnaSar(){
	//default to BasicDNA settings
	header.version = FORMAT_VERSION;
	format_version = FORMAT_VERSION;
}

BigDnaSar::BigDnaSar(const string& fname, const uint8* table, const uint32 alpha_bits){
	header.alphabet_bits = alpha_bits;
	memcpy(header.translation_table, table, UINT8_MAX);
	filename = fname;
	header.version = FORMAT_VERSION;
	format_version = FORMAT_VERSION;
}

BigDnaSar::BigDnaSar(const string& fname, const SuffixArray& sa){
}

BigDnaSar* BigDnaSar::Clone() const{
	BigDnaSar *bdsa = new BigDnaSar();
	(*bdsa) = *this;
	return bdsa;
}

uint64 BigDnaSar::GetNeededMemory(gnSeqI len){
	uint64 neededmem = 0;
	//forward and reverse copies of the sequence
	neededmem += len * 2;
	neededmem += sizeof(bmer) * len;
	return neededmem;
}

void BigDnaSar::SetMerMaskSize(uint32 mer_size){
	SuffixArray::SetMerMaskSize(mer_size);
	mer_mask |= 3;
}

uint32 BigDnaSar::CalculateMaxMerSize() const{
	return 62 / header.alphabet_bits;
}

void BigDnaSar::FillSar(gnSeqC* seq_buf, gnSeqI seq_len, boolean circular, vector<bmer>& s_array){
	DiskSuffixArray::FillDnaSar(seq_buf, seq_len, circular, s_array);
}


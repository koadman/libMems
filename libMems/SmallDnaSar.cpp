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
#include "SmallDnaSar.h"

uint32 SmallDnaSar::FORMAT_VERSION = 3;

SmallDnaSar::SmallDnaSar(){
	header.version = FORMAT_VERSION;
	format_version = FORMAT_VERSION;
}

SmallDnaSar::~SmallDnaSar(){
}

SmallDnaSar::SmallDnaSar(const string& fname, const uint8* table, const uint32 alpha_bits){
	header.alphabet_bits = alpha_bits;
	memcpy(header.translation_table, table, UINT8_MAX);
	filename = fname;
	header.version = FORMAT_VERSION;
	format_version = FORMAT_VERSION;
}
SmallDnaSar* SmallDnaSar::Clone() const{
	SmallDnaSar *bdsa = new SmallDnaSar();
	(*bdsa) = *this;
	return bdsa;
}

uint64 SmallDnaSar::GetNeededMemory(gnSeqI len){
	uint64 neededmem = (len * header.alphabet_bits) / 8;
	//forward and reverse copies of the sequence
	neededmem += len * 2;
	neededmem += sizeof(bmer) * len;
	return neededmem;
}

boolean SmallDnaSar::Convert(BigDnaSar& sa){
	SarHeader sa_head = sa.GetHeader();
	if(sa_head.mer_size * sa_head.alphabet_bits > 64){
		DebugMsg("SmallDnaSar::Create: Mer size is too large\n");
		return false;
	}

	//Allocate some memory
	uint32 SAR_BUFFER_SIZE = 32000;
	gnSeqI *pos_array = new gnSeqI[SAR_BUFFER_SIZE];
	if(pos_array == NULL){
		DebugMsg("SmallDnaSar::Create: Unable to allocate memory\n");
		return false;
	}

	// Open sarfile for writing
	if(!OpenForWriting())
		return false;

	header = sa_head;
	header.version = FORMAT_VERSION;
	sarfile.write((char*)&header, sizeof(struct SarHeader));
	if(!sarfile.good()){
		DebugMsg("SmallDnaSar::Create: Error writing header to disk.\n");
		sarfile.clear();
		sarfile.open(filename.c_str(), ios::binary | ios::in );
		return false;
	}
	//get binary sequence data and write it out
	gnSeqI binary_seq_len = (header.length * header.alphabet_bits) / 32;
	binary_seq_len++;	//fix for access violations.
	if(sequence != NULL)
		delete[] sequence;

	sequence = new uint32[binary_seq_len];
	if(sequence == NULL || !sa.GetBSequence(sequence, header.length)){
		DebugMsg("SmallDnaSar::Create: Error writing header to disk.\n");
		sarfile.clear();
		sarfile.open(filename.c_str(), ios::binary | ios::in );
		return false;
	}
	
	//write out the suffix array
	vector<bmer> s_array;
	s_array.reserve(SAR_BUFFER_SIZE);
	gnSeqI cur_pos = 0;
	while(cur_pos < header.length){
		sa.Read(s_array, SAR_BUFFER_SIZE, cur_pos);
		gnSeqI read_len = cur_pos;
		cur_pos += s_array.size();
		read_len = cur_pos - read_len;
		for(uint32 i=0; i < read_len; i++)
			pos_array[i] = s_array[i].position;

		sarfile.write((char*)&pos_array, sizeof(gnSeqI) * read_len);
		if(!sarfile.good()){
			DebugMsg("SmallDnaSar::Create: Error writing suffix array to disk.\n");
			sarfile.clear();
			sarfile.open(filename.c_str(), ios::binary | ios::in );
			delete[] pos_array;
			return false;
		}
	}
	delete[] pos_array;
	return true;
}

uint32 SmallDnaSar::CalculateMaxMerSize() const{
	return 62 / header.alphabet_bits;
}

void SmallDnaSar::FillSar(gnSeqC* seq_buf, gnSeqI seq_len, boolean circular, vector<bmer>& s_array){
	DiskSuffixArray::FillDnaSar(seq_buf, seq_len, circular, s_array);
}

uint64 SmallDnaSar::GetMer(gnSeqI offset){
	//check this for access violations.
	uint64 mer_a;
	gnSeqI mer_word, mer_bit;
	uint32 merle;
	//get mer_a
	mer_a = 0;
	mer_word = (offset * header.alphabet_bits) / 32;
	mer_bit = (offset * header.alphabet_bits) % 32;
	mer_a |= sequence[mer_word++];
	mer_a <<= 32;
	mer_a |= sequence[mer_word++];
	if(mer_bit > 0){
		merle = sequence[mer_word];
		merle >>= 32 - mer_bit;
		mer_a <<= mer_bit;
		mer_a |= merle;
	}
	mer_a &= mer_mask;
	
	//find the reverse complement of mer_a and return it if it's
	//smaller
	uint64 mer_b, mer_c = 0;	//mer_c will be the reverse complement
	mer_b = ~mer_a;
	uint32 masq = 0xffffffff;
	masq >>= 32 - header.alphabet_bits;
	for(uint32 i = 0; i < 64; i += header.alphabet_bits){
		mer_c |= mer_b & masq;
		mer_b >>= header.alphabet_bits;
		mer_c <<= header.alphabet_bits;
	}
	mer_c <<= 64 - (header.alphabet_bits * (header.mer_size+1));
	mer_c |= 1;

	return mer_a < mer_c ? mer_a : mer_c;
}



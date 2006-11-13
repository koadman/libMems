#include "MemorySuffixArray.h"

MemorySuffixArray::MemorySuffixArray(gnSequence& gps, uint32 mersize){
	if(gps.length() == 0)
		return;
	if(mersize * header.alphabet_bits > 64){
		DebugMsg("Requested mer size too large");
		return;
	}

// Initialization and setup
		
	gnSeqC* seq_buf;
	length = gps.length();
	alpha_bits = 2;
	mer_size = mersize;
	mer_bytes = mersize * alpha_bits / 8;
	if((mer_size * alpha_bits) % 8 != 0)
		mer_bytes++;
	description = "";
	trans = BasicDNATable;
	
// Memory allocation	
	
/*	trans = new gnSeqC[UINT8_MAX];
	if(trans == NULL){
		DebugMsg("Error allocating translation table.");
		return;
	}
*/
	orientation = new uint8[length/8 + 1];
	if(orientation == NULL){
		DebugMsg("Error allocating orientation table.");
		return;
	}

	seq_buf = new char[seq_len];
	if(seq_buf == NULL){
		DebugMsg("Error allocating sequence memory.");
		return;
	}

	sequence = new uint8[(length * alpha_bits) / 8 + 1];
	if(sequence == NULL){
		DebugMsg("Error allocating sequence memory.");
		return;
	}


	DebugMsg("Making memory suffix array");

// first read and convert the sequence data
	gps.ToArray(seq_buf, length);
	translate(sequence, seq_buf, length);
	delete seq_buf;
	
// fill in the suffix array	
	for(gnSeqI i = 0; i < length; i++){
		s_array[i].position = i;
		s_array[i].id = 1;
	}

// sort the suffix array
	qsort(0, length);

}
int compare(mmer a, mmer b){
	uint64 mer_a = sequence[(a.position * alpha_bits) / 8];
	if(a.position * alpha_bits
}


MemorySuffixArray::MemorySuffixArray(SuffixArray& sa){
	
	description = sa.description;
	alpha_bits = sa.alpha_bits;
	mer_size = sa.mer_size;
	mer_bytes = (mer_size * alpha_bits) / 8;
	if((mer_size * alpha_bits) % 8 != 0)
		mer_bytes++;
	length = sa.length;

	// allocate memory
	trans = new gnSeqC[UINT8_MAX];
	if(trans == NULL){
		DebugMsg("Error allocating translation table.");
		return;
	}
	memcpy(trans, sa.trans, UINT8_MAX);
	
	orientation = new uint8[length/8 + 1];
	if(orientation == NULL){
		DebugMsg("Error allocating orientation table.");
		return;
	}

	sequence = new uint8[(length * alpha_bits) / 8 + 1];
	if(sequence == NULL){
		DebugMsg("Error allocating sequence memory.");
		return;
	}
	
	// now copy the actual data in.

}

MemorySuffixArray::~MemorySuffixArray(){
	if(trans != NULL)
		delete trans;
	if(orientation != NULL)
		delete trans;
	if(sequence != NULL)
		delete trans;
}

MemorySuffixArray* MemorySuffixArray::Clone(){
	return new MemorySuffixArray(*this);
}
	
boolean MemorySuffixArray::Read(vector<bmer>& readVector, gnSeqI size, gnSeqI offset){
	if(offset >= length)
		return false;
	gnSeqI realend = size + offset < length ? size + offset : length;

	for(uint32 i=offset; i < realend; i++){
		readVector.push_back((*this)[i]);
	}
}

bmer MemorySuffixArray::operator[](gnSeqI index){
	bmer merle;
	merle.position = position[index];
	merle.mer = new uint8[mer_bytes];
	uint32 seqpos = (merle.position * alpha_bits) / 8;
	uint32 bitpos = (merle.position * alpha_bits) % 8;
	memcpy(merle.mer, sequence + seqpos, mer_bytes);
	//shift everything if necessary
	if(bitpos != 0){
		for(uint32 i=0; i < mer_bytes - 1; i++){
			merle.mer[i] <<= bitpos;
			uint8 mask = merle.mer[i+1] >>> 8-bitpos;
			merle.mer[i] |= mask;
		}
		merle.mer[mer_bytes-1] <<= bitpos;
	}
}

gnSeqI MemorySuffixArray::Find(gnSeqC* query_seq){
	bmer query_mer;
	//convert query_seq to mer form
	query_mer.mer = query_seq;

	return bsearch(query_mer, 0, length);
}


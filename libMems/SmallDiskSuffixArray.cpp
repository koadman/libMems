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
#include "SmallDiskSuffixArray.h"

uint32 SmallDiskSuffixArray::FORMAT_VERSION = 1;

SmallDiskSuffixArray::SmallDiskSuffixArray(){
	sequence = NULL;
	header.version = FORMAT_VERSION;
	format_version = FORMAT_VERSION;
}
SmallDiskSuffixArray::~SmallDiskSuffixArray(){
	if(sequence != NULL)
		delete[] sequence;
}

SmallDiskSuffixArray::SmallDiskSuffixArray(const string& fname, const uint8* table, const uint32 alpha_bits){
	
	header.alphabet_bits = alpha_bits;
	memcpy(header.translation_table, table, UINT8_MAX);
	filename = fname;
	header.version = FORMAT_VERSION;
	format_version = FORMAT_VERSION;
	sequence = NULL;
}

SmallDiskSuffixArray::SmallDiskSuffixArray(const string& fname, const SuffixArray& sa){
	sa, fname;
}

SmallDiskSuffixArray& SmallDiskSuffixArray::operator=(const SmallDiskSuffixArray& sa){
	DiskSuffixArray::operator=(sa);
	if(header.length > 0 && sa.sequence != NULL){
		// copy raw sequence data
		gnSeqI binary_seq_len = (header.length * header.alphabet_bits) / 32;
		binary_seq_len++;	//fix for access violations.
		sequence = new uint32[binary_seq_len];
		if(sequence == NULL){
			DebugMsg("SmallDnaSar::LoadFile: Error allocating memory.");
			return *this;
		}
		memcpy(sequence, sa.sequence, binary_seq_len * sizeof(uint32));
	}

	return *this;
}

SmallDiskSuffixArray* SmallDiskSuffixArray::Clone() const{
	SmallDiskSuffixArray *bdsa = new SmallDiskSuffixArray();
	(*bdsa) = *this;
	return bdsa;
}

boolean SmallDiskSuffixArray::LoadFile(const string& fname){
	filename = fname;
	sarfile.open(fname.c_str(), ios::binary | ios::in );
	if(!sarfile.is_open()){
		sarfile.clear();
		DebugMsg("SmallDiskSuffixArray::LoadFile: Unable to open file.\n");
		return false;
	}
	// read the header
	sarfile.read((char*)&header, sizeof(struct SarHeader));
	if(sarfile.gcount() < sizeof(struct SarHeader)){
		DebugMsg("SmallDiskSuffixArray::LoadFile: Unable to read file.");
		sarfile.clear();
		return false;
	}
	if(header.version != format_version){
		DebugMsg("SmallDiskSuffixArray::LoadFile: Unsupported file format.");
		return false;
	}

	mer_mask = UINT32_MAX;
	mer_mask <<= 32;
	mer_mask |= UINT32_MAX;
	mer_mask <<= (64 - header.alphabet_bits * header.mer_size);

	//header is ok.  read the sequence.
	gnSeqI seq_len = header.length;
	if(header.circular)
		seq_len += header.mer_size - 1;
	gnSeqI binary_seq_len = (seq_len * header.alphabet_bits) / 32;
	binary_seq_len++;	//fix for access violations.
	if(sequence != NULL)
		delete[] sequence;
	sequence = new uint32[binary_seq_len];
	if(sequence == NULL){
		DebugMsg("SmallDiskSuffixArray::LoadFile: Error allocating memory.");
		return false;
	}
	sarfile.read((char*)sequence, binary_seq_len*sizeof(uint32));
	if(sarfile.gcount() < binary_seq_len*sizeof(uint32)){
		DebugMsg("SmallDiskSuffixArray::LoadFile: Error reading sequence data.");
		sarfile.clear();
		return false;
	}

	sarray_start_offset = sarfile.tellg();
	sarfile.seekg(sarray_start_offset + sizeof(gnSeqI) * header.length);
	if(!sarfile.good()){
		DebugMsg("SmallDiskSuffixArray::LoadFile: Premature end of file.");
		sarfile.clear();
		return false;
	}
	filename = fname;

	return true;
}

uint64 SmallDiskSuffixArray::GetNeededMemory(gnSeqI len){
	uint64 neededmem = (len * header.alphabet_bits) / 8;
	neededmem += len;
	neededmem += sizeof(bmer) * len;
	return neededmem;
}

void SmallDiskSuffixArray::FillSar(gnSeqC* seq_buf, gnSeqI seq_len, boolean circular, vector<bmer>& s_array){
	DiskSuffixArray::FillNormalSar(seq_buf, seq_len, circular, s_array);
}

void SmallDiskSuffixArray::WriteSequence(gnSeqC* seq_buf, gnSeqI seq_len){
	gnSeqI binary_seq_len = (seq_len * header.alphabet_bits) / 32;
	if((seq_len * header.alphabet_bits) % 32 != 0)
		binary_seq_len++;

	sequence = new uint32[binary_seq_len];
	translate32(sequence, seq_buf, seq_len);
	sarfile.write((char*)sequence, binary_seq_len*sizeof(uint32));
}

void SmallDiskSuffixArray::WriteSar(vector<bmer>& s_array){
	uint32 pos_size = sizeof(gnSeqI);
	gnSeqI sar_len = s_array.size();
	for(gnSeqI suffixI=0; suffixI < sar_len; suffixI++)
		sarfile.write((char*)&(s_array[suffixI].position), pos_size);
}

//Merges the supplied suffix arrays into this one, overwriting the existing sar.
//KNOWN BUG:  The first suffix array must have (length * alphabet_bits) / word_bits == 0
//for Merge to work properly.
boolean SmallDiskSuffixArray::Merge(SuffixArray& sa, SuffixArray& sa2){
	SarHeader sa_head = sa.GetHeader();
	SarHeader sa_head2 = sa2.GetHeader();
	
	//allocate some memory
	const uint32 SEQ_BUFFER_SIZE = 200000;
	uint32* seq_buf = new uint32[SEQ_BUFFER_SIZE];
	if(seq_buf == NULL){
		DebugMsg("SmallDiskSuffixArray::Merge: Insufficient memory available.\n");
		return false;
	}
	//do some sanity checks on the sars we're merging.
	if(sa_head.alphabet_bits != sa_head2.alphabet_bits ||
	  sa_head.version != sa_head2.version ||
	  memcmp(sa_head.translation_table, sa_head2.translation_table, UINT8_MAX)){
		DebugMsg("SmallDiskSuffixArray::Merge: Incompatible suffix arrays.\n");
		delete[] seq_buf;
		return false;
	}
	
	
	if(!OpenForWriting()){
		delete[] seq_buf;
		return false;
	}

	header = sa_head;
	//take the smaller mer_size
	header.mer_size = sa_head.mer_size < sa_head2.mer_size ? sa_head.mer_size : sa_head2.mer_size;
	header.unique_mers = NO_UNIQUE_COUNT;
	header.length += sa_head2.length;
	mer_mask = UINT32_MAX;
	mer_mask <<= 32;
	mer_mask |= UINT32_MAX;
	mer_mask <<= (64 - header.alphabet_bits * header.mer_size);

	//write the header
	sarfile.write((char*)&header, sizeof(struct SarHeader));
	if(!sarfile.good()){
		DebugMsg("SmallDiskSuffixArray::Merge: Error writing suffix array header to disk.");
		sarfile.clear();
		sarfile.close();
		sarfile.open(filename.c_str(), ios::binary | ios::in );
		delete[] seq_buf;
		return false;
	}

	//copy sequence data into memory.
	uint32 binary_seq_len = (header.length * header.alphabet_bits) / 32;
	if((header.length * header.alphabet_bits) % 32 > 0)
		binary_seq_len++;

	//The +1 is to avoid access violations when copying in the
	//binary sequence before shifting.
	sequence = new uint32[binary_seq_len+1];
	sa.GetBSequence(sequence, sa_head.length, 0);

	uint32 bseq_len1 = (sa_head.length * sa_head.alphabet_bits) / 32;
	uint32 bseq_remainder = (sa_head.length * sa_head.alphabet_bits) % 32;
	if(bseq_remainder > 0){
		sa2.GetBSequence(&(sequence[bseq_len1]), sa_head2.length, 0);
		//mask off the end of the first sequence
		uint32 end_mask = 0xFFFFFFFF;
		end_mask <<= bseq_remainder;
		sequence[bseq_len1] &= end_mask;

		//shift the second sequence over.
		for(uint32 i=bseq_len1; i < binary_seq_len; i++){
			uint32 tmp = sequence[i+1];
			tmp >>= 32 - bseq_remainder;
			sequence[i] |= tmp;
			sequence[i+1] <<= bseq_remainder;
		}
	}else
		sa2.GetBSequence(&(sequence[bseq_len1]), sa_head2.length, 0);
	
	//write the sequence
	sarfile.write((char*)sequence, binary_seq_len * sizeof(uint32));
	sarray_start_offset = sarfile.tellg();

	//merge and write the suffix arrays
	vector<bmer> array1, array2;
	uint32 SAR_BUFFER_SIZE = SEQ_BUFFER_SIZE/2;  //actual size is this number * 13 bytes
	uint32 len1 = sa_head.length;
	uint32 len2 = sa_head2.length;
	uint32 k=0, l=0;
	while(k < len1 && l < len2){
		sa.Read(array1, SAR_BUFFER_SIZE, k);
		sa2.Read(array2, SAR_BUFFER_SIZE, l);
		//mergesort them
		uint32 m = 0, n = 0;
		while(m < array1.size() && n < array2.size()){
			if(array1[m].mer <= array2[n].mer){
				seq_buf[m+n] = array1[m].position;
				m++;
			}else{
				seq_buf[m+n] = array2[n].position + sa_head.length;
				n++;
			}
		}
		for(;m < array1.size(); m++)
			seq_buf[m+n] = array1[m].position;
		for(;n < array2.size(); n++)
			seq_buf[m+n] = array2[n].position;
		sarfile.write((char*)seq_buf, (m + n)*sizeof(uint32));
		k += SAR_BUFFER_SIZE;
		l += SAR_BUFFER_SIZE;
	}
	delete[] seq_buf;

	if(!sarfile.good()){
		DebugMsg("SmallDiskSuffixArray::Merge: Error writing position array.");
		sarfile.clear();
		sarfile.close();
		sarfile.open(filename.c_str(), ios::binary | ios::in );
		return false;
	}
	// reopen the suffix array file read-only
	sarfile.close();
	sarfile.open(filename.c_str(), ios::binary | ios::in );
	if(!sarfile.is_open()){
		DebugMsg("SmallDiskSuffixArray::Merge: Error opening suffix array file.");
		sarfile.clear();
		return false;
	}
	return true;
}

boolean SmallDiskSuffixArray::Read(vector<bmer>& readVector, gnSeqI size, const gnSeqI offset){
	if(!sarfile.is_open()){
		DebugMsg("SmallDiskSuffixArray::Read: Error sar file not open.\n");
		return false;
	}
	gnSeqI total_len = header.length;
	if(offset >= total_len)
		return false;
	gnSeqI readlen = offset + size < total_len ? size : total_len - offset;

	uint64 seekpos = sarray_start_offset;
	//now seekpos contains the index of the first byte of the suffix
	//array.  seek to where we want to read
	seekpos += offset * sizeof(gnSeqI);
	sarfile.seekg(seekpos);
	
	//allocate memory to read
	readVector.clear();
	gnSeqI* pos_array = new gnSeqI[readlen];
	if(pos_array == NULL){
		DebugMsg("SmallDiskSuffixArray::Read: Out of memory.\n");
		return false;
	}
	readVector.reserve(readlen);

	//read
	sarfile.read((char*)pos_array, sizeof(gnSeqI) * readlen);
	gnSeqI success_count = sarfile.gcount() / sizeof(gnSeqI);

	//copy data to the vector
	for(uint32 j=0; j < success_count; j++){
		bmer tmp_mer;
		tmp_mer.position = pos_array[j];
		tmp_mer.mer = GetMer(pos_array[j]);
		readVector.push_back(tmp_mer);
	}
	delete[] pos_array;

	if( success_count < readlen){
		DebugMsg("SmallDiskSuffixArray::Read: Error reading data.\n");
		sarfile.clear();
		return false;
	}
	return true;
}

uint64 SmallDiskSuffixArray::GetMer(gnSeqI offset){
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
	return mer_a;
}

//potential buffer overflows here.  make dest extra big.
boolean SmallDiskSuffixArray::GetBSequence(uint32* dest, const gnSeqI len, const gnSeqI offset){
	//first determine the byte offset of the sequence within the file.
	if(offset >= header.length){
		DebugMsg("SmallDiskSuffixArray::GetBSequence: Offset out of range.");
		return false;
	}
	uint64 startpos = (offset * header.alphabet_bits) / 32;
	int begin_remainder = (offset * header.alphabet_bits) % 32;
	uint64 readlen = offset + len < header.length ? len : header.length - offset;

	gnSeqI word_read_len = (readlen * header.alphabet_bits) / 32;
	int end_remainder = (readlen * header.alphabet_bits) % 32;
	if(begin_remainder + (readlen * header.alphabet_bits) > 32
	   && end_remainder > 0)
		word_read_len++;
	if(begin_remainder > 0)
		word_read_len++;
	
	//now do the actual read
	memcpy((char*)dest, (char*)sequence + (startpos * 4), word_read_len * 4);
	
	//now shift if needed
	ShiftWords(dest, word_read_len, -begin_remainder);
	
	//now mask if needed
	if(end_remainder > begin_remainder){
		uint32 mask = 0xFFFFFFFF;
		mask <<= 32 - (end_remainder - begin_remainder);
		dest[word_read_len-1] &= mask;
	}else if(end_remainder < begin_remainder){
		uint32 mask = 0xFFFFFFFF;
		mask <<= (begin_remainder - end_remainder);
		dest[word_read_len-2] &= mask;
	}
	return true;
}

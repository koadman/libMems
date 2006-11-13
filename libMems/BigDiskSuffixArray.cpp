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
#include "BigDiskSuffixArray.h"

uint32 BigDiskSuffixArray::FORMAT_VERSION = 0;

BigDiskSuffixArray::BigDiskSuffixArray(){
	header.version = FORMAT_VERSION;
	format_version = FORMAT_VERSION;
}

BigDiskSuffixArray::BigDiskSuffixArray(const string& fname, const uint8* table, const uint32 alpha_bits){
	header.alphabet_bits = alpha_bits;
	memcpy(header.translation_table, table, UINT8_MAX);
	filename = fname;
	header.version = FORMAT_VERSION;
	format_version = FORMAT_VERSION;
}

BigDiskSuffixArray::BigDiskSuffixArray(const string& fname, const SuffixArray& sa){
}

BigDiskSuffixArray* BigDiskSuffixArray::Clone() const{
	BigDiskSuffixArray *bdsa = new BigDiskSuffixArray();
	(*bdsa) = *this;
	return bdsa;
}

boolean BigDiskSuffixArray::LoadFile(const string& fname){
	sarfile.open(fname.c_str(), ios::in | ios::binary );
	filename = fname;
	if(!sarfile.is_open()){
		DebugMsg("BigDiskSuffixArray::LoadFile: Unable to open file.\n");
		return false;
	}
	//file already exists.  try to use it.
	sarfile.read((char*)&header, sizeof(struct SarHeader));
	if(sarfile.gcount() < sizeof(struct SarHeader)){
		DebugMsg("BigDiskSuffixArray::LoadFile: Unable to read file.");
		sarfile.clear();
		sarfile.close();
		return false;
	}
	if(header.version != format_version){
		DebugMsg("BigDiskSuffixArray::LoadFile: Unsupported file format.");
		sarfile.close();
		return false;
	}
	mer_mask = UINT32_MAX;
	mer_mask <<= 32;
	mer_mask |= UINT32_MAX;
	mer_mask <<= (64 - header.alphabet_bits * header.mer_size);

	sarray_start_offset = sarfile.tellg();
	uint64 total_len = header.length;
	sarray_start_offset += (total_len * header.alphabet_bits) / 8;
	if((total_len * header.alphabet_bits) % 8)
		sarray_start_offset++;
	sarfile.seekg(sarray_start_offset + sizeof(bmer) * header.length);
	if(!sarfile.good()){
		DebugMsg("BigDiskSuffixArray::LoadFile: Premature end of file.");
		sarfile.clear();
		sarfile.close();
		return false;
	}
	return true;
}

uint64 BigDiskSuffixArray::GetNeededMemory(gnSeqI len){
	uint64 neededmem = 0;
	neededmem += len;
	neededmem += sizeof(bmer) * len;
	return neededmem;
}

void BigDiskSuffixArray::FillSar(gnSeqC* seq_buf, gnSeqI seq_len, boolean circular, vector<bmer>& s_array){
	DiskSuffixArray::FillNormalSar(seq_buf, seq_len, circular, s_array);
}


uint64  BigDiskSuffixArray::GetMer(gnSeqI position){
	ErrorMsg("BigDiskSuffixArray::GetMer: Not yet implemented.\n");
	return 0;
}

//Merges the supplied suffix arrays into this one, overwriting the existing sar.
boolean BigDiskSuffixArray::Merge(SuffixArray& sa, SuffixArray& sa2){
	SarHeader sa_head = sa.GetHeader();
	SarHeader sa_head2 = sa2.GetHeader();
	
	//allocate some memory
	const uint32 SEQ_BUFFER_SIZE = 200000;
	uint32* seq_buf = new uint32[SEQ_BUFFER_SIZE];
	if(seq_buf == NULL){
		DebugMsg("BigDiskSuffixArray::Merge: Insufficient memory available.\n");
		return false;
	}
	bmer* sar_buf = new bmer[SEQ_BUFFER_SIZE];
	if(sar_buf == NULL){
		DebugMsg("BigDiskSuffixArray::Merge: Insufficient memory available.\n");
		delete[] seq_buf;
		return false;
	}
	//do some sanity checks on the sars we're merging.
	if(sa_head.alphabet_bits != sa_head2.alphabet_bits ||
	  sa_head.version != sa_head2.version ||
	  memcmp(sa_head.translation_table, sa_head2.translation_table, UINT8_MAX)){
		DebugMsg("BigDiskSuffixArray::Merge: Incompatible suffix arrays.\n");
		delete[] seq_buf;
		delete[] sar_buf;
		return false;
	}
	
	
	// Open sarfile for writing
	if(!OpenForWriting()){
		delete[] seq_buf;
		delete[] sar_buf;
		return false;
	}
	header = sa_head;
	//take the smaller mer_size
	header.mer_size = sa_head.mer_size < sa_head2.mer_size ? sa_head.mer_size : sa_head2.mer_size;
	header.unique_mers = NO_UNIQUE_COUNT;
	header.length += sa_head2.length;

	//write the header
	sarfile.write((char*)&header, sizeof(struct SarHeader));
	if(!sarfile.good()){
		DebugMsg("BigDiskSuffixArray::Merge: Error writing suffix array header to disk.\n");
		sarfile.close();
		sarfile.open(filename.c_str(), ios::binary | ios::in );
		delete[] seq_buf;
		delete[] sar_buf;
		return false;
	}
	
	//write the first sequence
	gnSeqI seqI = 0;
	while(seqI < sa_head.length){
		gnSeqI read_size = sa_head.length - seqI < SEQ_BUFFER_SIZE ? sa_head.length - seqI : SEQ_BUFFER_SIZE;
		sa.GetBSequence(seq_buf, read_size, seqI);
		gnSeqI binary_seq_len = (read_size * header.alphabet_bits) / 32;
		if(binary_seq_len == 0){
			DebugMsg("BigDiskSuffixArray::Merge: Alignment bug encountered, file will be corrupt.\n");
			break;
		}
		sarfile.write((char*)seq_buf, binary_seq_len * sizeof(uint32));
		seqI += (binary_seq_len * 32) / header.alphabet_bits;
	}
	//write the second sequence
	seqI = 0;
	while(seqI < sa_head2.length){
		uint32 read_size = sa_head2.length - seqI < SEQ_BUFFER_SIZE ? sa_head2.length - seqI : SEQ_BUFFER_SIZE;
		sa2.GetBSequence(seq_buf, read_size, seqI);
		gnSeqI binary_seq_len = (read_size * header.alphabet_bits) / 32;
		if(binary_seq_len == 0)
			binary_seq_len++;
		sarfile.write((char*)seq_buf, binary_seq_len * sizeof(uint32));
		seqI += (binary_seq_len * 32) / header.alphabet_bits;
	}
	delete[] seq_buf;

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
				sar_buf[m+n] = array1[m];
				m++;
			}else{
				sar_buf[m+n] = array2[n];
				sar_buf[m+n].position += sa_head.length;
				n++;
			}
		}
		for(;m < array1.size(); m++)
			sar_buf[m+n] = array1[m];
		for(;n < array2.size(); n++)
			sar_buf[m+n] = array2[n];
		sarfile.write((char*)sar_buf, (m + n)*sizeof(bmer));
		if(sarfile.gcount() != (m + n)*sizeof(bmer)){
			DebugMsg("BigDiskSuffixArray::Merge: Error writing bmer array.");
			sarfile.close();
			sarfile.open(filename.c_str(), ios::binary | ios::in );
			delete[] sar_buf;
			return false;
		}
		k += SAR_BUFFER_SIZE;
		l += SAR_BUFFER_SIZE;
	}
	delete[] sar_buf;

	// reopen the suffix array file read-only
	sarfile.close();
	sarfile.open(filename.c_str(), ios::binary | ios::in );
	if(!sarfile.is_open()){
		DebugMsg("BigDiskSuffixArray::Merge: Error opening suffix array file.");
		return false;
	}
	return true;
}

boolean BigDiskSuffixArray::Read(vector<bmer>& readVector, gnSeqI size, const gnSeqI offset){
	gnSeqI total_len = header.length;
	if(offset >= total_len)
		return false;
	gnSeqI readlen = offset + size < total_len ? size : total_len - offset;
	uint64 seekpos = sarray_start_offset;
	//now seekpos contains the index of the first byte of the suffix
	//array.  seek to where we want to read
	uint64 beemer_size = sizeof(gnSeqI) + sizeof(uint64);
	seekpos += offset * beemer_size;
	sarfile.seekg(seekpos);
	
	//allocate memory to read
	readVector.clear();
	uint64 read_size = readlen * beemer_size;
	char* mer_array = new char[read_size];
	readVector.reserve(readlen);

	//read
	sarfile.read((char*)mer_array, read_size);
	gnSeqI success_count = sarfile.gcount() / beemer_size;

	//copy data to the vector
	for(uint32 j=0; j < success_count; j++){
		bmer tmp_mer;
		tmp_mer.position = ((gnSeqI*)(mer_array + j * beemer_size))[0];
		tmp_mer.mer = ((uint64*)(mer_array + j * beemer_size + sizeof(gnSeqI)))[0] &= mer_mask;
		readVector.push_back(tmp_mer);
	}
	delete[] mer_array;

	if( success_count < readlen)
		return false;
	return true;
}

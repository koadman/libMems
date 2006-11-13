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

#include "DiskSuffixArray.h"
#include "gnFilter.h"
#include <algorithm>

uint64 DiskSuffixArray::MEMORY_MINIMUM = DEFAULT_MEMORY_MINIMUM;

DiskSuffixArray& DiskSuffixArray::operator=(const DiskSuffixArray& sa){
	header = sa.header;
	filename = sa.filename;
	sarray_start_offset = sa.sarray_start_offset;
	mer_mask = sa.mer_mask;
	format_version = sa.format_version;
//	if(sa.sarfile.is_open()){
		sarfile.open(filename.c_str(), ios::binary | ios::in );
		if(!sarfile.is_open()){
			DebugMsg("DiskSuffixArray::=: Unable to open suffix array file.");
			sarfile.clear();
			return *this;
		}
//	}
	return *this;
}

boolean DiskSuffixArray::OpenForWriting(){
	// Open sarfile for writing
	boolean was_open = sarfile.is_open();
	if(was_open)
		sarfile.close();
	sarfile.open(filename.c_str(), ios::binary | ios::in | ios::out | ios::trunc );
	if(!sarfile.is_open() || !sarfile.good()){
		DebugMsg("DiskSuffixArray::OpenForWriting: Unable to open file for writing.");
		sarfile.clear();
		if(was_open)
			sarfile.open(filename.c_str(), ios::binary | ios::in );
		return false;
	}
	return true;
}

boolean DiskSuffixArray::WriteHeader(){
	if(!sarfile.is_open()){
		DebugMsg("DiskSuffixArray::WriteHeader: File is not valid.\n");
		return false;
	}
	boolean success = true;
	// Open sarfile for writing and write new header.
	sarfile.close();
	sarfile.open(filename.c_str(), ios::binary | ios::in | ios::out);
	if(!sarfile.is_open()){
		DebugMsg("DiskSuffixArray::WriteHeader: Unable to open file for writing.\n");
		success = false;
	}else{
		sarfile.write((char*)&header, sizeof(struct SarHeader));
		if(!sarfile.good()){
			DebugMsg("DiskSuffixArray::WriteHeader: Error writing header to disk.\n");
			success = false;
		}
	}
	// reopen the suffix array file read-only
	sarfile.close();
	sarfile.open(filename.c_str(), ios::binary | ios::in );
	if(!sarfile.is_open()){
		DebugMsg("DiskSuffixArray::WriteHeader: Error opening suffix array file.\n");
		success = false;
	}
	return success;
}

uint32 DiskSuffixArray::UniqueMerCount(){
	uint32 tmp_count = header.unique_mers;
	SuffixArray::UniqueMerCount();
	if(tmp_count != header.unique_mers)
		WriteHeader();
	return header.unique_mers;
}

//change the description in memory and on disk
void DiskSuffixArray::SetDescription(const string& d){
	strncpy(header.description, d.c_str(), DESCRIPTION_SIZE-1);
	WriteHeader();
} 

void DiskSuffixArray::SetID(const sarID_t d){
	header.id = d;
	WriteHeader();
}
//read binary sequence data from the file
//only works with byte aligned offsets for now.

boolean DiskSuffixArray::GetBSequence(uint32* dest, const gnSeqI len, const gnSeqI offset){
	//first determine the byte offset of the sequence within the file.
	uint64 startpos = sizeof(SarHeader);
	if(offset >= header.length){
		DebugMsg("DiskSuffixArray::GetBSequence: Offset out of range.");
		return false;
	}
	if((offset * header.alphabet_bits) % 8 != 0){
		DebugMsg("DiskSuffixArray::GetBSequence: Offset not byte aligned.");
		return false;
	}
	startpos += (offset * header.alphabet_bits) / 8;

	uint64 readlen = offset + len < header.length ? len : header.length - offset;

	gnSeqI byte_read_len = (readlen * header.alphabet_bits) / 8;
	if((readlen * header.alphabet_bits) % 8 > 0)
		byte_read_len++;
	
	//now do the actual read
	sarfile.seekg(startpos);
	sarfile.read((char*)dest, byte_read_len);
	if(sarfile.gcount() == byte_read_len)
		return true;
	return false;
}

boolean DiskSuffixArray::Create(const gnSequence& seq, const int8 id, const uint32 mersize){
	vector<bmer> s_array;
	gnSeqC* seq_buf;

	if(CalculateMaxMerSize() == 0){
		DebugMsg("DiskSuffixArray::Create: Alphabet size is too large\n");
		return false;
	}
	if(mersize > CalculateMaxMerSize()){
		DebugMsg("DiskSuffixArray::Create: Mer size is too large\n");
		return false;
	}
	gnSeqI mer_size = mersize;
	if(mer_size == 0){
		DebugMsg("DiskSuffixArray::Create: Cant have 0 mer size\n");
		return false;
	}
	if(!OpenForWriting())
		return false;
	//determine sequence and sar length and read in sequence
	gnSeqI seq_len = seq.length();
	gnSeqI sar_len = seq_len;
	if(seq.isCircular())
		seq_len += mer_size - 1;
	
	seq_buf = new char[seq_len];
	seq.ToArray(seq_buf, seq_len);
	if(seq.isCircular())
		seq.ToArray(seq_buf+sar_len, mer_size-1);

	// set header information
	header.version = format_version;
	header.description[0] = 0;
	header.length = sar_len;
	header.mer_size = mer_size;
	header.id = id;

	SetMerMaskSize(mersize);

	FillSar(seq_buf, seq_len, seq.isCircular(), s_array);
	sort(s_array.begin(), s_array.end(), &bmer_lessthan);
	
	/* now write out the file header */
	sarfile.write((char*)&header, sizeof(struct SarHeader));

	if(!sarfile.good()){
		DebugMsg("BigDiskSuffixArray::Create: Error writing suffix array header to disk.\n");
		delete[] seq_buf;
		sarfile.clear();
		return false;
	}

	/* write out the actual sequence */
	WriteSequence(seq_buf, seq_len);
	delete[] seq_buf;
	sarray_start_offset = sarfile.tellg();

	/* write out the suffix array */
	WriteSar(s_array);
	
	sarfile.flush();
	if(!sarfile.good()){
		DebugMsg("BigDiskSuffixArray::Create: Error writing suffix array to disk.\n");
		sarfile.clear();
		return false;
	}
	// reopen the suffix array file read-only
	sarfile.close();
	sarfile.open(filename.c_str(), ios::binary | ios::in );
	if(!sarfile.is_open()){
		DebugMsg("BigDiskSuffixArray::Create: Error opening suffix array file.\n");
		return false;
	}
	return true;
}
uint32 DiskSuffixArray::CalculateMaxMerSize() const{
	return 64 / header.alphabet_bits;
}

void DiskSuffixArray::WriteSequence(gnSeqC* seq_buf, gnSeqI seq_len){
	gnSeqI binary_seq_len = (seq_len * header.alphabet_bits) / 32;
	if((seq_len * header.alphabet_bits) % 32 != 0)
		binary_seq_len++;

	uint32* sequence = new uint32[binary_seq_len];
	translate32(sequence, seq_buf, seq_len);
	sarfile.write((char*)sequence, binary_seq_len*sizeof(uint32));
	delete[] sequence;
}

void DiskSuffixArray::WriteSar(vector<bmer>& s_array){
	uint32 pos_size = sizeof(gnSeqI);
	uint32 mer_size = sizeof(uint64);
	gnSeqI sar_len = s_array.size();
	for(gnSeqI suffixI=0; suffixI < sar_len; suffixI++){
		sarfile.write((char*)&(s_array[suffixI].position), pos_size);
		sarfile.write((char*)&(s_array[suffixI].mer), mer_size);
	}
}

void DiskSuffixArray::FillNormalSar(gnSeqC* seq_buf, gnSeqI seq_len, boolean circular, vector<bmer>& s_array){
	uint32 alpha_bits = header.alphabet_bits;
	uint32 mer_size = header.mer_size;
	gnSeqI sar_len = seq_len;
	if(circular)
		sar_len -= ( header.mer_size - 1);
	s_array.reserve(sar_len);

	bmer cur_suffix;
	cur_suffix.mer = 0;
	cur_suffix.position = 0;

	/* now fill in the suffix array with the forward sequence*/
	for(gnSeqI i=0; i < mer_size; i++){
		cur_suffix.mer <<= alpha_bits;
		cur_suffix.mer |= header.translation_table[seq_buf[i]];
	}
	uint8 dead_bits = 64 - (mer_size * alpha_bits);
	cur_suffix.mer <<= dead_bits;

	s_array.push_back(cur_suffix);

	//fill s_array with mers
	for(gnSeqI seqI = 1; seqI < seq_len; seqI++){//already added the
													//first one
		cur_suffix.position++;
		cur_suffix.mer <<= alpha_bits;
		uint64 new_mer = header.translation_table[seq_buf[seqI+(mer_size-1)]];
		new_mer <<= dead_bits;
		cur_suffix.mer |= new_mer;
		s_array.push_back(cur_suffix);
	}
}

void DiskSuffixArray::FillDnaSar(gnSeqC* seq_buf, gnSeqI seq_len, boolean circular, vector<bmer>& s_array){
	/* now fill in the suffix array with the forward sequence*/
	uint32 alpha_bits = header.alphabet_bits;
	uint32 mer_size = header.mer_size;
	gnSeqI sar_len = seq_len;
	if(circular)
		sar_len -= ( header.mer_size - 1);
	s_array.reserve(sar_len);

	gnSeqI binary_seq_len = (seq_len * header.alphabet_bits) / 32;
	if((seq_len * header.alphabet_bits) % 32 != 0)
		binary_seq_len++;
	
	uint32* seq = new uint32[binary_seq_len];
	translate32(seq, seq_buf, seq_len);

	uint32 dead_bits = 64 - (mer_size * alpha_bits);
	uint64 create_mask = UINT32_MAX;
	create_mask <<= 32;
	create_mask |= UINT32_MAX;
	create_mask <<= dead_bits;

	bmer cur_suffix, rcur_suffix;
	cur_suffix.mer = seq[0];
	cur_suffix.mer <<= 32;
	cur_suffix.mer |= seq[1];
	cur_suffix.mer &= create_mask;
	cur_suffix.position = 0;
	rcur_suffix.mer = 0;
	rcur_suffix.position = 0;
	
	//find the reverse complement of cur_suffix.mer and return it if it's
	//smaller
	uint64 mer_b = 0;
	mer_b = ~cur_suffix.mer;
	uint32 masq = 0xffffffff;
	masq >>= 32 - alpha_bits;
	for(uint32 i = 0; i < 64; i += alpha_bits){
		rcur_suffix.mer |= mer_b & masq;
		mer_b >>= alpha_bits;
		rcur_suffix.mer <<= alpha_bits;
	}
	rcur_suffix.mer <<= dead_bits - alpha_bits;
	rcur_suffix.mer |= 1;

	//add the first mer
	if(cur_suffix.mer < rcur_suffix.mer)
		s_array.push_back(cur_suffix);
	else
		s_array.push_back(rcur_suffix);

	//fill s_array with mers
	gnSeqI endI;
	if(circular)
		endI = sar_len + mer_size - 1;
	else
		endI = sar_len;
	uint32 rdead_bits = 64 - alpha_bits - dead_bits;
	uint64 tmp_rseq = 0;
	uint32 seqI = (mer_size * alpha_bits) / 32;
	int32 cur_bit = 32 - alpha_bits - ((mer_size * alpha_bits) % 32);
	uint32 cur_seq = seq[seqI];
	uint32 tmp_seq;
	uint32 alpha_mask = 0xFFFFFFFF;
	alpha_mask >>= 32 - alpha_bits;
	uint32 revalpha_mask = alpha_mask;
	revalpha_mask <<= dead_bits;

	//which is slower? a memory operation or a conditional?
	//probably a memory operation.
	for(gnSeqI cur_pos = mer_size; cur_pos < endI; cur_pos++){//already added the
													//first one
		//increment positions
		cur_suffix.position++;
		rcur_suffix.position++;
		
		//extract the next character
		tmp_seq = cur_seq;
		tmp_seq >>= cur_bit;
		tmp_seq &= alpha_mask;
		tmp_seq <<= dead_bits;
		
		//add it to the forward mer
		cur_suffix.mer <<= alpha_bits;
		cur_suffix.mer |= tmp_seq;

		//do the reverse complement mer
		tmp_seq = ~tmp_seq;
		tmp_seq &= revalpha_mask;
		tmp_rseq = tmp_seq;
		tmp_rseq <<= rdead_bits;
		rcur_suffix.mer >>= alpha_bits;
		rcur_suffix.mer |= tmp_rseq;
		rcur_suffix.mer &= create_mask;
		rcur_suffix.mer |= 1;

		if(cur_suffix.mer < rcur_suffix.mer)
			s_array.push_back(cur_suffix);
		else
			s_array.push_back(rcur_suffix);

		cur_bit -= alpha_bits;
		if(cur_bit < 0){
			cur_bit += alpha_bits;
			cur_seq <<= 16;		//trade bitwise ops for conditional
			cur_seq <<= 16 - (cur_bit);
			seqI++;
			tmp_seq = seq[seqI];
			tmp_seq >>= cur_bit;
			cur_seq |= tmp_seq;
			cur_bit += 32 - alpha_bits;
		}
	}
}

struct bmer DiskSuffixArray::operator[](const gnSeqI index){
	bmer merle;
	vector<bmer> the_mers;
	Read(the_mers, 1, index);
	the_mers.push_back(merle);
	return the_mers[0];
}

boolean DiskSuffixArray::BigCreate(const gnSequence& seq, const sarID_t id, const uint32 mersize){
	unsigned long freemem = wxGetFreeMemory();	//get the amount of free memory.
	unsigned long neededmem = GetNeededMemory(seq.length());
	if(neededmem >= freemem && neededmem > MEMORY_MINIMUM){ // divide and conquer
		boolean success = true;
		gnSeqI midpoint = seq.length() / 2;
		midpoint = (midpoint * header.alphabet_bits) / 32;
		midpoint = (midpoint / header.alphabet_bits) * 32;
		gnSequence seqA = seq.subseq(1, midpoint);
		gnSequence seqB = seq.subseq(1 + midpoint, seq.length() - midpoint);
		seqA.setCircular(false);
		seqB.setCircular(false);
		cout << "Splitting " << seq.length() << " to " << seqA.length() << " and " << seqB.length() << "\n";
		//create the first sar
		wxString temp_sarfile = wxGetTempFileName("bdsa_split");
//		wxRemoveFile(temp_sarfile);
		DiskSuffixArray* temp_sar = this->Clone();
		temp_sar->filename = temp_sarfile.c_str();
		success = success && temp_sar->BigCreate(seqA, id, mersize);
		//create the second sar
		wxString temp_sarfile2 = wxGetTempFileName("bdsa_split");
//		wxRemoveFile(temp_sarfile2);
		DiskSuffixArray* temp_sar2 = this->Clone();
		temp_sar2->filename = temp_sarfile2.c_str();
		success = success && temp_sar2->BigCreate(seqB, id, mersize);
		//merge them to this file
		cout << "Merging " << seqA.length() << " and " << seqB.length() << "\n";
		success = success && Merge(*temp_sar, *temp_sar2);
		//free up RAM
		delete temp_sar;
		delete temp_sar2;
		//erase the temp files.
		wxRemoveFile(temp_sarfile);
		wxRemoveFile(temp_sarfile2);
		return success;
	}else{
		return Create(seq, id, mersize);
	}
	return false;

}

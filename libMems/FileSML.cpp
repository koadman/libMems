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

#include "FileSML.h"
#include "gnFilter.h"
#include <algorithm>
#include <cmath>


FileSML& FileSML::operator=(const FileSML& sa)
{
 	SortedMerList::operator=( sa );
 	filename = sa.filename;
	sarray_start_offset = sa.sarray_start_offset;
	sarfile.open(filename.c_str(), ios::binary | ios::in );
	if(!sarfile.is_open()){
		DebugMsg("FileSML::=: Unable to open suffix array file.\n");
		sarfile.clear();
		return *this;
	}
	return *this;
}

void FileSML::LoadFile(const string& fname){
	filename = fname;
	sarfile.open(fname.c_str(), ios::binary | ios::in );
	if(!sarfile.is_open()){
		sarfile.clear();
		Throw_gnExMsg( FileNotOpened(), "Unable to open file.\n");
	}
	// read the header
	sarfile.read((char*)&header, sizeof(struct SMLHeader));
	if(sarfile.gcount() < sizeof(struct SMLHeader)){
		sarfile.clear();
		Throw_gnExMsg( FileUnreadable(), "Unable to read file.");
	}
	if(header.version != FormatVersion()){
		Throw_gnExMsg( FileUnreadable(), "Unsupported file format.");
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

	sarfile.read((char*)sequence, binary_seq_len*sizeof(uint32));
	if(sarfile.gcount() < binary_seq_len*sizeof(uint32)){
		sarfile.clear();
		Throw_gnExMsg( FileUnreadable(), "Error reading sequence data.");
	}

	sarray_start_offset = sarfile.tellg();
	sarfile.seekg(sarray_start_offset + sizeof(gnSeqI) * header.length);
	if(!sarfile.good()){
		sarfile.clear();
		Throw_gnExMsg( FileUnreadable(), "Premature end of file.");
	}
	filename = fname;

}

void FileSML::OpenForWriting(){
	// Open smlfile for writing
	boolean was_open = sarfile.is_open();
	if(was_open)
		sarfile.close();
	sarfile.open(filename.c_str(), ios::binary | ios::in | ios::out | ios::trunc );
	if(!sarfile.is_open() || !sarfile.good()){
		sarfile.clear();
		if(was_open)
			sarfile.open(filename.c_str(), ios::binary | ios::in );
		Throw_gnExMsg(FileNotOpened(), "Unable to open file for writing.");
	}
}

boolean FileSML::WriteHeader(){
	if(!sarfile.is_open()){
		Throw_gnExMsg(IOStreamFailed(), "File is not valid.");
	}
	boolean success = true;
	char* errormsg;
	// Open sarfile for writing and write new header.
	OpenForWriting();
	sarfile.write((char*)&header, sizeof(struct SMLHeader));
	if(!sarfile.good()){
		errormsg = "Error writing header to disk.";
		success = false;
	}

	// reopen the sorted mer list file read-only
	sarfile.close();
	sarfile.open(filename.c_str(), ios::binary | ios::in );
	if(!sarfile.is_open()){
		errormsg = "Error opening sorted mer list file.";
		success = false;
	}
	if(!success)
		Throw_gnExMsg(IOStreamFailed(), errormsg);
	return success;
}

uint32 FileSML::UniqueMerCount(){
	uint32 tmp_count = header.unique_mers;
	SortedMerList::UniqueMerCount();
	if(tmp_count != header.unique_mers)
		WriteHeader();
	return header.unique_mers;
}

//change the description in memory and on disk
void FileSML::SetDescription(const string& d){
	strncpy(header.description, d.c_str(), DESCRIPTION_SIZE-1);
	WriteHeader();
} 

void FileSML::SetID(const sarID_t d){
	header.id = d;
	WriteHeader();
}

void FileSML::Create(const gnSequence& seq, const uint32 mersize){
	try{
		SortedMerList::Create( seq, mersize );
		OpenForWriting();
	
		vector<bmer> sml_array;
		FillSML(seq, sml_array);

//	RadixSort(s_array);
	sort(sml_array.begin(), sml_array.end(), &bmer_lessthan);
	
	/* now write out the file header */
	sarfile.write((char*)&header, sizeof(struct SMLHeader));

	if(!sarfile.good()){
		sarfile.clear();
		Throw_gnExMsg( IOStreamFailed(), "Error writing sorted mer list header to disk.\n");
	}

	/* write out the actual sequence */
	sarfile.write((char*)sequence, binary_seq_len*sizeof(uint32));
	sarray_start_offset = sarfile.tellg();

	/* write out the sorted mer list */
	for(gnSeqI suffixI=0; suffixI < sml_array.size(); suffixI++)
		sarfile.write((char*)&(sml_array[suffixI].position), sizeof(gnSeqI));
	
	sarfile.flush();
	if(!sarfile.good()){
		sarfile.clear();
		Throw_gnExMsg( IOStreamFailed(), "Error writing sorted mer list to disk.\n");
	}
	// reopen the sorted mer list file read-only
	sarfile.close();
	sarfile.open(filename.c_str(), ios::binary | ios::in );
	if(!sarfile.is_open())
		Throw_gnExMsg( FileNotOpened(), "FileSML::Create: Error opening sorted mer list file.\n");

	}catch( bad_alloc& ba ){
		// If there wasn't enough memory then divide and conquer
//		if( sequence != NULL )
//			delete[] sequence;
		BigCreate( seq, 1, mersize );
		return;
	}catch( exception& e ){
//		if( sequence != NULL )
//			delete[] sequence;
		cerr << e.what() << endl;
		BigCreate( seq, 1, mersize );
		return;
	}catch( gnException& gne ){
//		if( sequence != NULL )
//			delete[] sequence;
		cerr << gne << endl;
		BigCreate( seq, 1, mersize );
		return;
	}
}

bmer FileSML::operator[](const gnSeqI index){
	vector<bmer> the_mers;
	Read(the_mers, 1, index);
	return the_mers[0];
}


boolean FileSML::Read(vector<bmer>& readVector, gnSeqI size, const gnSeqI offset){
	if(!sarfile.is_open()){
		DebugMsg("FileSML::Read: Error sar file not open.\n");
		return false;
	}
	gnSeqI total_len = header.length;
	if(offset >= total_len)
		return false;
	gnSeqI readlen = offset + size < total_len ? size : total_len - offset;
	
	wxMutexLocker wxml(*file_mutex);
	
	uint64 seekpos = sarray_start_offset;
	//now seekpos contains the index of the first byte of the sorted mer list
	//seek to where we want to read
	seekpos += offset * sizeof(gnSeqI);
	sarfile.seekg(seekpos);
	
	//allocate memory to read
	readVector.clear();
	Array<gnSeqI> pos_array( readlen );
	readVector.reserve(readlen);

	//read
	sarfile.read((char*)pos_array.data, sizeof(gnSeqI) * readlen);
	gnSeqI success_count = sarfile.gcount() / sizeof(gnSeqI);

	//copy data to the vector
	for(uint32 j=0; j < success_count; j++){
		bmer tmp_mer;
		tmp_mer.position = pos_array.data[j];
		tmp_mer.mer = GetMer(pos_array.data[j]);
		readVector.push_back(tmp_mer);
	}

	if( success_count < readlen){
		// couldn't read the entire request for some reason -- usually EOF.
		if(!sarfile.eof())
			ErrorMsg("Error reading from file.\n");
		sarfile.clear();
		return false;
	}
	return true;
}

void FileSML::BigCreate(const gnSequence& seq, const uint32 split_levels, const uint32 mersize){
//	unsigned long freemem = wxGetFreeMemory();	//get the amount of free memory.
//	unsigned long neededmem = GetNeededMemory(seq.length());
//	if(neededmem >= freemem && neededmem > MEMORY_MINIMUM){ // divide and conquer
	if(split_levels > 0){	// split_levels defines the number of times to divide and conquer
		uint64 midpoint = seq.length() / 2;
		midpoint = (midpoint * header.alphabet_bits) / 32;
		midpoint = (midpoint / header.alphabet_bits) * 32;
		gnSequence seqA = seq.subseq(1, midpoint);
		gnSequence seqB = seq.subseq(1 + midpoint, seq.length() - midpoint);
		seqA.setCircular(false);
		seqB.setCircular(false);
		cout << "Splitting " << seq.length() << " to " << seqA.length() << " and " << seqB.length() << "\n";

		//create the first sar
		wxString temp_sarfile = wxGetTempFileName("bdsa_split");
		FileSML* temp_sar = this->Clone();
		temp_sar->filename = temp_sarfile.c_str();
		temp_sar->BigCreate(seqA, split_levels - 1, mersize);

		//create the second sar
		wxString temp_sarfile2 = wxGetTempFileName("bdsa_split");
		FileSML* temp_sar2 = this->Clone();
		temp_sar2->filename = temp_sarfile2.c_str();
		temp_sar2->BigCreate(seqB, split_levels - 1, mersize);

		//merge them to this file
		cout << "Merging " << seqA.length() << " and " << seqB.length() << "\n";
		Merge(*temp_sar, *temp_sar2);
		//free up RAM
		delete temp_sar;
		delete temp_sar2;
		//erase the temp files.
		wxRemoveFile(temp_sarfile);
		wxRemoveFile(temp_sarfile2);
	}else{
		Create(seq, mersize);
	}
}

void FileSML::RadixSort(vector<bmer>& s_array){
	vector<bmer> *source_buckets;
	vector<bmer> *tmp_buckets;
	vector<bmer> *buckets;
	uint32 radix_size = 11;
	uint64 radix_mask = 0xFFFFFFFF;
	radix_mask <<= 32;
	radix_mask |= 0xFFFFFFFF;
	radix_mask >>= 64 - radix_size;
	
	uint32 bucket_count = (uint32) pow((double)2, (double)radix_size);
	uint32 cur_shift_bits = 0;
	buckets = new vector<bmer>[bucket_count];
	source_buckets = new vector<bmer>[bucket_count];
	uint64 cur_bucket;
	for(uint32 merI = 0; merI < s_array.size(); merI++){
		cur_bucket = s_array[merI].mer & radix_mask;
		source_buckets[cur_bucket].push_back(s_array[merI]);
	}
	s_array.clear();
	cur_shift_bits += radix_size;
	radix_mask <<= radix_size;
	while(cur_shift_bits < 64){
		for(uint32 bucketI = 0; bucketI < bucket_count; bucketI++){
			for(uint32 merI = 0; merI < source_buckets[bucketI].size(); merI++){
				cur_bucket = source_buckets[bucketI][merI].mer & radix_mask;
				cur_bucket >>= cur_shift_bits;
				buckets[cur_bucket].push_back(source_buckets[bucketI][merI]);
			}
			source_buckets[bucketI].clear();
		}
		cur_shift_bits += radix_size;
		radix_mask <<= radix_size;
		tmp_buckets = source_buckets;
		source_buckets = buckets;
		buckets = tmp_buckets;
	}
	s_array.clear();
	for(uint32 bucketI = 0; bucketI < bucket_count; bucketI++){
		for(uint32 merI = 0; merI < source_buckets[bucketI].size(); merI++){
			s_array.push_back(source_buckets[bucketI][merI]);
		}
		source_buckets[bucketI].clear();
	}
	delete[] source_buckets;
	delete[] buckets;
}

//Merges the supplied sorted mer lists into this one, overwriting the existing sml.
//KNOWN BUG:  The first sorted mer list must have (length * alphabet_bits) / word_bits == 0
//for Merge to work properly.
void FileSML::Merge(SortedMerList& sa, SortedMerList& sa2){
STACK_TRACE_START
	SMLHeader sa_head = sa.GetHeader();
	SMLHeader sa_head2 = sa2.GetHeader();
	
	//basic copying
	header = sa_head;
	//take the smaller mer_size
	if(sa_head.mer_size < sa_head2.mer_size){
		header.mer_size = sa_head.mer_size;
		mer_mask = sa.GetMerMask();
	}else{
		header.mer_size = sa_head2.mer_size;
		mer_mask = sa2.GetMerMask();
	}
	header.unique_mers = NO_UNIQUE_COUNT;
	header.length += sa_head2.length;

	//allocate some memory
	const uint32 SEQ_BUFFER_SIZE = 200000;
	Array<uint32> seq_buf ( SEQ_BUFFER_SIZE + header.mer_size );

	//do some sanity checks on the sars we're merging.
	if(sa_head.alphabet_bits != sa_head2.alphabet_bits ||
	  sa_head.version != sa_head2.version ||
	  memcmp(sa_head.translation_table, sa_head2.translation_table, UINT8_MAX)){
		Throw_gnExMsg(SMLMergeError(), "Incompatible sorted mer lists.");
	}
	
	OpenForWriting();

	//write the header
	sarfile.write((char*)&header, sizeof(struct SMLHeader));
	if(!sarfile.good()){
		sarfile.clear();
		sarfile.close();
		sarfile.open(filename.c_str(), ios::binary | ios::in );
		Throw_gnExMsg(IOStreamFailed(), "Error writing sorted mer list header to disk.");
	}

	//copy sequence data into memory.
	uint32 binary_seq_len = (header.length * header.alphabet_bits) / 32;
	if((header.length * header.alphabet_bits) % 32 > 0)
		binary_seq_len++;

	//The +1 is to avoid access violations when copying in the
	//binary sequence before shifting.
	if( sequence != NULL )
		delete[] sequence;
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

	//get new mers in the middle
	vector<bmer> middle_mers;
	bmer mid_mer;
	for(uint32 midI = sa_head.length - header.mer_size + 1; midI < sa_head.length; midI++){
		mid_mer.position = midI;
		mid_mer.mer = GetMer(midI);
		middle_mers.push_back(mid_mer);
	}
	sort(middle_mers.begin(), middle_mers.end(), &bmer_lessthan);
	//put a special mer at the end which will never go into the sorted mer list
	//since every possible mer is less than it.
	mid_mer.mer = 0xFFFFFFFF;
	mid_mer.mer <<= 32;
	mid_mer.mer |= 0xFFFFFFFF;
	mid_mer.position = GNSEQI_END;
	middle_mers.push_back(mid_mer);
	//merge and write the sorted mer lists
	vector<bmer> array1, array2;
	uint32 SAR_BUFFER_SIZE = SEQ_BUFFER_SIZE/2;  //actual size is this number * 13 bytes
	uint32 k=0, l=0, midI=0;
	uint32 m = 0, n = 0;
	gnSeqI bufferI=0;
	do{
		//mergesort them
		while(m < array1.size() && n < array2.size()){
			if(array1[m].mer <= array2[n].mer){
				if(array1[m].mer <= middle_mers[midI].mer){
					seq_buf.data[bufferI] = array1[m].position;
					m++;
					bufferI++;
				}else{
					seq_buf.data[bufferI] = middle_mers[midI].position;
					midI++;
					bufferI++;
				}
			}else if(array2[n].mer <= middle_mers[midI].mer){
				seq_buf.data[bufferI] = array2[n].position + sa_head.length;
				n++;
				bufferI++;
			}else{
				seq_buf.data[bufferI] = middle_mers[midI].position;
				midI++;
				bufferI++;
			}
			if(bufferI == SEQ_BUFFER_SIZE){
				sarfile.write((char*)seq_buf.data, bufferI * sizeof(uint32));
				bufferI = 0;
			}
		}
		if(m == array1.size()){
			sa.Read(array1, SAR_BUFFER_SIZE, k);
			k += array1.size();
			m = 0;
		}
		if(n == array2.size()){
			sa2.Read(array2, SAR_BUFFER_SIZE, l);
			l += array2.size();
			n = 0;
		}
	}while(array1.size() != 0 && array2.size() != 0);
	if(bufferI > 0)
		sarfile.write((char*)seq_buf.data, (bufferI)*sizeof(uint32));
	//consolidate the remaining mers to a known vector
	vector<bmer> remaining_mers;
	for(;m < array1.size(); m++)
		remaining_mers.push_back(array1[m]);
	for(;n < array2.size(); n++){
		remaining_mers.push_back(array2[n]);
		remaining_mers[remaining_mers.size()-1].position += sa_head.length;
	}
	for(;midI < middle_mers.size() - 1; midI++)
		remaining_mers.push_back(middle_mers[midI]);
	//merge them with the remaining middle_mers
	sort(remaining_mers.begin(), remaining_mers.end(), &bmer_lessthan);
	uint32 remI = 0;
	for(;remI < remaining_mers.size(); remI++)
		seq_buf.data[remI] = remaining_mers[remI].position;
	if(remI > 0)
		sarfile.write((char*)seq_buf.data, (remI)*sizeof(uint32));

	if(!sarfile.good()){
		sarfile.clear();
		sarfile.close();
		sarfile.open(filename.c_str(), ios::binary | ios::in );
		Throw_gnExMsg(IOStreamFailed(), "Error writing position array.");
	}
	// reopen the sorted mer list file read-only
	sarfile.close();
	sarfile.open(filename.c_str(), ios::binary | ios::in );
	if(!sarfile.is_open()){
		sarfile.clear();
		Throw_gnExMsg(FileNotOpened(), "Error opening sorted mer list file.");
	}
STACK_TRACE_END
}

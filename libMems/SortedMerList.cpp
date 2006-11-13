#include "SortedMerList.h"

const uint8* const SortedMerList::BasicDNATable(){
	static const uint8* const bdt = SortedMerList::CreateBasicDNATable();
	return bdt;
}

const uint8* const SortedMerList::ProteinTable(){
	static const uint8* const bdt = SortedMerList::CreateProteinTable();
	return bdt;
}

const uint8* const SortedMerList::CreateBasicDNATable(){
	uint8* bdt = new uint8[UINT8_MAX];
	memset(bdt, 0, UINT8_MAX);
	bdt['c'] = 1;
	bdt['C'] = 1;
	bdt['b'] = 1;
	bdt['B'] = 1;
	bdt['y'] = 1;
	bdt['Y'] = 1;
	bdt['g'] = 2;
	bdt['G'] = 2;
	bdt['s'] = 2;
	bdt['S'] = 2;
	bdt['k'] = 2;
	bdt['K'] = 2;
	bdt['t'] = 3;
	bdt['T'] = 3;
	return bdt;
}

const uint8* const SortedMerList::CreateProteinTable(){
	uint8* pt = new uint8[UINT8_MAX];
	memset(pt, 0, UINT8_MAX);
	pt['A'] = 0;
	pt['R'] = 1;
	pt['N'] = 2;
	pt['D'] = 3;
	pt['C'] = 4;
	pt['Q'] = 5;
	pt['E'] = 6;
	pt['G'] = 7;
	pt['H'] = 8;
	pt['I'] = 9;
	pt['L'] = 10;
	pt['K'] = 11;
	pt['M'] = 12;
	pt['F'] = 13;
	pt['P'] = 14;
	pt['S'] = 15;
	pt['T'] = 16;
	pt['W'] = 17;
	pt['Y'] = 18;
	pt['V'] = 19;
	
	pt['a'] = 0;
	pt['r'] = 1;
	pt['n'] = 2;
	pt['d'] = 3;
	pt['c'] = 4;
	pt['q'] = 5;
	pt['e'] = 6;
	pt['g'] = 7;
	pt['h'] = 8;
	pt['i'] = 9;
	pt['l'] = 10;
	pt['k'] = 11;
	pt['m'] = 12;
	pt['f'] = 13;
	pt['p'] = 14;
	pt['s'] = 15;
	pt['t'] = 16;
	pt['w'] = 17;
	pt['y'] = 18;
	pt['v'] = 19;
	return pt;
}

SortedMerList::SortedMerList(){
	//default to BasicDNA settings
	header.length = 0;
	header.alphabet_bits = 2;
	header.unique_mers = NO_UNIQUE_COUNT;
	memcpy(header.translation_table, BasicDNATable(), UINT8_MAX);
	header.description[0] = 0;
	header.mer_size = DNA_MER_SIZE;
	header.id = 0;
	header.circular = false;
	mask_size = DNA_MER_SIZE;
	// init sequence data to null
	sequence = NULL;
	binary_seq_len = 0;
}

SortedMerList::SortedMerList( const SortedMerList& sa ){
	sequence = NULL;
	*this = sa;
}

SortedMerList& SortedMerList::operator=(const SortedMerList& sa)
{
	header = sa.header;
	mer_mask = sa.mer_mask;
	mask_size = sa.mask_size;
	binary_seq_len = sa.binary_seq_len;

	// copy binary sequence data
	if( sa.sequence != NULL ){
		if( sequence != NULL )
			delete[] sequence;
		sequence = new uint32[binary_seq_len];
		memcpy(sequence, sa.sequence, sizeof(uint32) * binary_seq_len);
	}else
		sequence = NULL;

	return *this;
}

SortedMerList::~SortedMerList(){
	if( sequence != NULL )
		delete[] sequence;
}

uint32 SortedMerList::CalculateMaxMerSize() const{
	bmer tmp;
	return (sizeof(tmp.mer) * 8) / header.alphabet_bits;
}

boolean SortedMerList::FindMer(const uint64 query_mer, gnSeqI& result){
	bmer merle;
	merle.mer = query_mer;
	result = bsearch(merle, 0, Length() - 1);
	return ((*this)[result].mer == merle.mer);
}

boolean SortedMerList::Find(const string& query_seq, gnSeqI& result) {
	struct bmer merle;
	merle.mer = 0;

	//check the length to make sure it is small enough
	gnSeqI len = query_seq.length() * header.alphabet_bits < 64 ? 
		query_seq.length() : 64 / header.alphabet_bits;
		
	translate((uint8*)&merle.mer, query_seq.c_str(), len);
	result = bsearch(merle, 0, Length() - 1);
	return ((*this)[result].mer == merle.mer);
}

void SortedMerList::FindAll(const string& query_seq, vector<gnSeqI> result) {
	struct bmer merle;
	merle.mer = 0;

	//check the length to make sure it is small enough
	gnSeqI len = query_seq.length() * header.alphabet_bits < 64 ? 
		query_seq.length() : 64 / header.alphabet_bits;
		
	translate((uint8*)&merle.mer, query_seq.c_str(), len);
	
	//find the first match then start filling forward.
	gnSeqI matchI = 0;
	gnSeqI curend = Length() - 1;
	bmer matchmer;
	matchI = bsearch(merle, 0, curend);

	//first seek backwards
	int64 cur_matchI = matchI;
	matchmer = (*this)[matchI];
	while(cur_matchI >= 0 && matchmer.mer == merle.mer){
		cur_matchI--;
		matchmer = (*this)[cur_matchI];
	}
	int64 first_matchI = cur_matchI+1;

	//now seek forwards
	cur_matchI = matchI+1;
	matchmer = (*this)[cur_matchI];
	while(cur_matchI < GNSEQI_END && matchmer.mer == merle.mer){
		cur_matchI++;
		matchmer = (*this)[cur_matchI];
	}
	//fill the result array
	for(matchI = first_matchI; matchI < cur_matchI; matchI++)
		result.push_back(matchI);
}

string SortedMerList::Description() const{
	return header.description;
}

void SortedMerList::SetDescription(const string& d){
	strncpy(header.description, d.c_str(), DESCRIPTION_SIZE-1);
}

uint32 SortedMerList::MerSize() const{
	return header.mer_size;
}

boolean SortedMerList::IsCircular() const{
	return header.circular;
}

uint64 SortedMerList::GetMerMask() const{
	return mer_mask;
}

uint32 SortedMerList::GetMerMaskSize() const{
	return mask_size;
}

void SortedMerList::SetMerMaskSize(uint32 mer_size){
	if(mer_size > header.mer_size)
		mask_size = header.mer_size;
	else
		mask_size = mer_size;

	// calculate the mer mask
	mer_mask = UINT32_MAX;
	mer_mask <<= 32;
	mer_mask |= UINT32_MAX;
	mer_mask <<= (64 - header.alphabet_bits * header.mer_size);
}

gnSeqI SortedMerList::Length() const{
	return header.length;
}

sarID_t SortedMerList::GetID() const{
	return header.id;
}
void SortedMerList::SetID(const sarID_t d){
	header.id = d;
}

void SortedMerList::SetSequence(gnSeqC* seq_buf, gnSeqI seq_len){
	binary_seq_len = (seq_len * header.alphabet_bits) / 32;
	if((seq_len * header.alphabet_bits) % 32 != 0)
		binary_seq_len++;

	if( sequence != NULL )
		delete[] sequence;
	sequence = new uint32[binary_seq_len];
	translate32(sequence, seq_buf, seq_len);
}

uint64 SortedMerList::GetMer(gnSeqI position){
	//check this for access violations.
	uint64 mer_a;
	gnSeqI mer_word, mer_bit;
	uint32 merle;
	//get mer_a
	mer_a = 0;
	mer_word = (position * header.alphabet_bits) / 32;
	mer_bit = (position * header.alphabet_bits) % 32;
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
void SortedMerList::GetBSequence(uint32* dest, const gnSeqI len, const gnSeqI offset){
	//first determine the byte offset of the sequence within the file.
	if(offset >= header.length){
		Throw_gnEx( IndexOutOfBounds() );
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
}

gnSeqI SortedMerList::bsearch(const struct bmer& query_mer, const gnSeqI start, const gnSeqI end) {

	gnSeqI middle = (start + end) / 2;
	struct bmer midmer = (*this)[middle];
	if(midmer.mer == query_mer.mer)
		return middle;
	else if((midmer.mer < query_mer.mer) && (middle < end))
		return bsearch(query_mer, middle + 1, end);
	else if((midmer.mer > query_mer.mer) && (start < middle))
		return bsearch(query_mer, start, middle - 1);
	
	//if we get here then the mer was not found.
	//return where it would be if it existed.
	return middle;
}

//translate the character sequence to binary form based on the
//translation table.
void SortedMerList::translate(uint8* dest, const gnSeqC* src, const gnSeqI len) const{
	uint8 start_bit = 0;
	gnSeqI cur_byte = 0;
	uint32 alpha_bits = header.alphabet_bits;
	dest[cur_byte] = 0;
	for(uint32 i=0; i < len; i++){
		uint8 tmp = header.translation_table[src[i]];
		if(start_bit + alpha_bits <= 8){
			tmp <<= 8 - start_bit - alpha_bits;
			dest[cur_byte] |= tmp;
		}else{
			uint8 over_bits = (start_bit + alpha_bits) % 8;
			uint8 tmp2 = tmp;
			tmp2 <<= 8 - over_bits;
			tmp >>= over_bits;
			dest[cur_byte] |= tmp;
			dest[cur_byte+1] |= tmp2;
		}
		start_bit += alpha_bits;
		if(start_bit >= 8){
			start_bit %= 8;
			cur_byte++;
			dest[cur_byte] = 0;
		}
	}
}

void SortedMerList::translate32(uint32* dest, const gnSeqC* src, const gnSeqI len) const{
	if( len == 0 )
		return;
	uint8 start_bit = 0;
	gnSeqI cur_word = 0;
	uint32 alpha_bits = header.alphabet_bits;
	dest[cur_word] = 0;
	for(uint32 i=0; i < len; i++){
		uint32 tmp = header.translation_table[src[i]];
		if(start_bit + alpha_bits <= 32){
			tmp <<= 32 - start_bit - alpha_bits;
			dest[cur_word] |= tmp;
			start_bit += alpha_bits;
			if(start_bit >= 32 && i < len - 1){
				start_bit %= 32;
				cur_word++;
				dest[cur_word] = 0;
			}
		}else{
			uint8 over_bits = (start_bit + alpha_bits) % 32;
			uint32 tmp2 = tmp;
			tmp2 <<= 32 - over_bits;
			tmp >>= over_bits;
			dest[cur_word] |= tmp;
			cur_word++;
			dest[cur_word] = 0;
			dest[cur_word] |= tmp2;
			start_bit = over_bits;
		}
	}
}
SMLHeader SortedMerList::GetHeader() const{
	return header;
}

gnSeqI SortedMerList::UniqueMerCount(){
	if(header.unique_mers != NO_UNIQUE_COUNT)
		return header.unique_mers;

	uint32 MER_BUFFER_SIZE = 16384;  //not quite arbitrary (2^14)
	uint32 cur_pos = 0;
	vector<bmer> mer_vector;
	bmer prev_mer;
	uint32 m_unique = 0;
	while(cur_pos < header.length){
		if(!Read(mer_vector, MER_BUFFER_SIZE, cur_pos)){
			DebugMsg("SortedMerList::UniqueMerCount: Error reading bmer vector.");
			return NO_UNIQUE_COUNT;
		}
		uint32 mer_count = mer_vector.size();
		if(mer_count == 0)
			break;
		if(cur_pos > 0 && prev_mer.mer != mer_vector[0].mer)
			m_unique++;
		
		//count them up.
		uint32 i = 0;
		for(uint32 j = 1; j < mer_count; j++){
			if(mer_vector[i].mer != mer_vector[j].mer)
				m_unique++;
			i++;
		}
		prev_mer = mer_vector[i];
	}
	m_unique++;
	header.unique_mers = m_unique;
	return header.unique_mers;
}

//will not handle more than 8GB sequence on 32-bit systems
void SortedMerList::ShiftWords(unsigned int* data, uint32 length, int32 bits){
	int32 word_bits = 8 * sizeof(unsigned int);
	if(bits > 0 && bits < word_bits){
		//shift everything right starting at the end
		data[length - 1] >>= bits;
		for(int i=length-2; i >= 0; i--){
			uint32 tmp = data[i];
			tmp <<= word_bits - bits;
			data[i+1] |= tmp;
			data[i] >>= bits;
		}
	}else if(bits < 0 && bits > (-1)*word_bits){
		bits *= -1;
		//shift everything left
		data[0] <<= bits;
		for(uint32 i=0; i < length; i++){
			uint32 tmp = data[i+1];
			tmp >>= word_bits - bits;
			data[i] |= tmp;
			data[i+1] <<= bits;
		}
	}
}

void SortedMerList::FillSML(gnSeqC* seq_buf, gnSeqI seq_len, boolean circular, vector<bmer>& sml_array){
	uint32 alpha_bits = header.alphabet_bits;
	uint32 mer_size = header.mer_size;
	gnSeqI sar_len = seq_len;
	if(!circular)
		sar_len -= header.mer_size - 1;
	sml_array.reserve(sar_len);

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

	sml_array.push_back(cur_suffix);

	//fill sml_array with mers
	for(gnSeqI seqI = 1; seqI < sar_len; seqI++){//already added the
													//first one
		cur_suffix.position++;
		cur_suffix.mer <<= alpha_bits;
		uint64 new_mer = header.translation_table[seq_buf[seqI+(mer_size-1)]];
		new_mer <<= dead_bits;
		cur_suffix.mer |= new_mer;
		sml_array.push_back(cur_suffix);
	}
}

void SortedMerList::FillSML(const gnSequence& seq, vector<bmer>& sml_array){
	gnSeqI seq_len = seq.length();
	Array<gnSeqC> seq_buf( seq_len );
	seq.ToArray(seq_buf.data, seq_len);
	FillSML(seq_buf.data, seq_len, seq.isCircular(), sml_array);
}

void SortedMerList::FillSML(gnSeqI seq_len, vector<gnSeqI>& pos_array){
	pos_array.clear();
	pos_array.reserve( seq_len );
	for(gnSeqI seqI = 0; seqI < seq_len; seqI++ )
		pos_array.push_back(seqI);
}

uint64 SortedMerList::GetDnaMer(gnSeqI offset){
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
	
	// for debugging
	if( mer_c < mer_a )
		return mer_c;
	return mer_a < mer_c ? mer_a : mer_c;
}

void SortedMerList::FillDnaSML(const gnSequence& seq, vector<bmer>& sml_array){
	/* now fill in the suffix array with the forward sequence*/
	uint32 alpha_bits = header.alphabet_bits;
	uint32 mer_size = header.mer_size;
	gnSeqI sar_len = seq.length();
	if( sar_len < header.mer_size )
		return;	// can't have an sml if there ain't enough sequence
	if( !seq.isCircular() )
		sar_len -= ( header.mer_size - 1);
	sml_array.reserve(sar_len);

	uint32 dead_bits = 64 - (mer_size * alpha_bits);
	uint64 create_mask = UINT32_MAX;
	create_mask <<= 32;
	create_mask |= UINT32_MAX;
	create_mask <<= dead_bits;

	bmer cur_suffix, rcur_suffix;
	cur_suffix.mer = sequence[0];
	cur_suffix.mer <<= 32;
	cur_suffix.mer |= sequence[1];
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
		sml_array.push_back(cur_suffix);
	else
		sml_array.push_back(rcur_suffix);

	//fill sml_array with mers
	gnSeqI 	endI = sar_len + mer_size;
	if(seq.isCircular())
		endI += mer_size;

	uint32 rdead_bits = 64 - alpha_bits - dead_bits;
	uint64 tmp_rseq = 0;
	uint32 seqI = (mer_size * alpha_bits) / 32;
	int32 cur_bit = 32 - alpha_bits - ((mer_size * alpha_bits) % 32);
	uint32 cur_seq = sequence[seqI];
	uint64 tmp_seq;
	uint32 alpha_mask = 0xFFFFFFFF;
	alpha_mask >>= 32 - alpha_bits;
	uint64 revalpha_mask = alpha_mask;
	revalpha_mask <<= dead_bits;

	//which is slower? a memory operation or a conditional?
	//probably a memory operation.
	for(gnSeqI cur_pos = mer_size + 1; cur_pos < endI; cur_pos++){//already added the
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
			sml_array.push_back(cur_suffix);
		else
			sml_array.push_back(rcur_suffix);

		cur_bit -= alpha_bits;
		if(cur_bit < 0){
			cur_bit += alpha_bits;
			cur_seq <<= 16;		//trade bitwise ops for conditional
			cur_seq <<= 16 - (cur_bit);
			seqI++;
			tmp_seq = sequence[seqI];
			tmp_seq >>= cur_bit;
			cur_seq |= tmp_seq;
			cur_bit += 32 - alpha_bits;
		}
	}
}

void SortedMerList::Create(const gnSequence& seq, const uint32 mersize){
	
	if(CalculateMaxMerSize() == 0)
		Throw_gnExMsg( SMLCreateError(), "Alphabet size is too large" );
	if(mersize > CalculateMaxMerSize())
		Throw_gnExMsg( SMLCreateError(), "Mer size is too large" );

	gnSeqI mer_size = mersize;
	if(mer_size == 0)
		Throw_gnExMsg( SMLCreateError(), "Can't have 0 mer size" );

	//determine sequence and sar length and read in sequence
	gnSeqI seq_len = seq.length();
	gnSeqI sar_len = seq_len;
	if(!seq.isCircular()){
		sar_len -= mer_size;
		header.circular = false;
	}else
		header.circular = true;
	// use the nifty Array class as a wrapper for the buffer to ensure correct deallocation
	gnSeqI buf_len = seq.isCircular() ? seq_len + mer_size : seq_len;
	Array<gnSeqC> seq_buf( buf_len );
	seq.ToArray(seq_buf.data, seq_len);
	if( seq.isCircular() )
		seq.ToArray(seq_buf.data + seq_len, mer_size-1);

	// set header information
	header.length = seq_len;
	header.mer_size = mer_size;

	SetMerMaskSize( mersize );
	SetSequence( seq_buf.data, buf_len );
}

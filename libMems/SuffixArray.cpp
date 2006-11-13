#include "SuffixArray.h"

uint8* SuffixArray::BasicDNATable = SuffixArray::CreateBasicDNATable();
uint8* SuffixArray::ProteinTable = SuffixArray::CreateProteinTable();

uint8* SuffixArray::CreateBasicDNATable(){
	BasicDNATable = new uint8[UINT8_MAX];
	memset(BasicDNATable, 0, UINT8_MAX);
	BasicDNATable['c'] = 1;
	BasicDNATable['C'] = 1;
	BasicDNATable['b'] = 1;
	BasicDNATable['B'] = 1;
	BasicDNATable['y'] = 1;
	BasicDNATable['Y'] = 1;
	BasicDNATable['g'] = 2;
	BasicDNATable['G'] = 2;
	BasicDNATable['s'] = 2;
	BasicDNATable['S'] = 2;
	BasicDNATable['k'] = 2;
	BasicDNATable['K'] = 2;
	BasicDNATable['t'] = 3;
	BasicDNATable['T'] = 3;
	return BasicDNATable;
}
uint8* SuffixArray::CreateProteinTable(){
	ProteinTable = new uint8[UINT8_MAX];
	memset(ProteinTable, 0, UINT8_MAX);
	ProteinTable['A'] = 0;
	ProteinTable['R'] = 1;
	ProteinTable['N'] = 2;
	ProteinTable['D'] = 3;
	ProteinTable['C'] = 4;
	ProteinTable['Q'] = 5;
	ProteinTable['E'] = 6;
	ProteinTable['G'] = 7;
	ProteinTable['H'] = 8;
	ProteinTable['I'] = 9;
	ProteinTable['L'] = 10;
	ProteinTable['K'] = 11;
	ProteinTable['M'] = 12;
	ProteinTable['F'] = 13;
	ProteinTable['P'] = 14;
	ProteinTable['S'] = 15;
	ProteinTable['T'] = 16;
	ProteinTable['W'] = 17;
	ProteinTable['Y'] = 18;
	ProteinTable['V'] = 19;
	
	ProteinTable['a'] = 0;
	ProteinTable['r'] = 1;
	ProteinTable['n'] = 2;
	ProteinTable['d'] = 3;
	ProteinTable['c'] = 4;
	ProteinTable['q'] = 5;
	ProteinTable['e'] = 6;
	ProteinTable['g'] = 7;
	ProteinTable['h'] = 8;
	ProteinTable['i'] = 9;
	ProteinTable['l'] = 10;
	ProteinTable['k'] = 11;
	ProteinTable['m'] = 12;
	ProteinTable['f'] = 13;
	ProteinTable['p'] = 14;
	ProteinTable['s'] = 15;
	ProteinTable['t'] = 16;
	ProteinTable['w'] = 17;
	ProteinTable['y'] = 18;
	ProteinTable['v'] = 19;
	return ProteinTable;
}

SuffixArray::SuffixArray(){
	//default to BasicDNA settings
	header.length = 0;
	header.alphabet_bits = 2;
	header.unique_mers = NO_UNIQUE_COUNT;
	memcpy(header.translation_table, BasicDNATable, UINT8_MAX);
	header.description[0] = 0;
	header.mer_size = DNA_MER_SIZE;
	header.id = 0;
	header.circular = false;
	mask_size = DNA_MER_SIZE;
}

gnSeqI SuffixArray::Find(const string& query_seq) {
	struct bmer merle;
	merle.mer = 0;

	//check the length to make sure it is small enough
	gnSeqI len = query_seq.length() * header.alphabet_bits < 64 ? 
		query_seq.length() : 64 / header.alphabet_bits;
		
	translate((uint8*)&merle.mer, query_seq.c_str(), len);
	return bsearch(merle, 0, Length() - 1);
}

void SuffixArray::FindAll(const string& query_seq, vector<gnSeqI> result) {
	struct bmer merle;
	merle.mer = 0;

	//check the length to make sure it is small enough
	gnSeqI len = query_seq.length() * header.alphabet_bits < 64 ? 
		query_seq.length() : 64 / header.alphabet_bits;
		
	translate((uint8*)&merle.mer, query_seq.c_str(), len);
	
	//find the first match then start filling forward.
	gnSeqI matchI = 0;
	gnSeqI curend = Length() - 1;
	gnSeqI prevmatch = GNSEQI_END;
	do{
		matchI = bsearch(merle, 0, curend);
		if(matchI != GNSEQI_END)
			prevmatch = matchI;
		curend = matchI - 1;
	}while(matchI != GNSEQI_END && matchI != 0);

	//check whether anything was found
	if(prevmatch == GNSEQI_END)
		return;

	//fill the result array
	while((*this)[prevmatch].mer == merle.mer)
		result.push_back(prevmatch++);
}

string SuffixArray::Description() const{
	return header.description;
}

void SuffixArray::SetDescription(const string& d){
	strncpy(header.description, d.c_str(), DESCRIPTION_SIZE-1);
}

uint32 SuffixArray::MerSize() const{
	return header.mer_size;
}

boolean SuffixArray::IsCircular() const{
	return header.circular;
}

uint64 SuffixArray::GetMerMask() const{
	return mer_mask;
}

uint32 SuffixArray::GetMerMaskSize() const{
	return mask_size;
}
void SuffixArray::SetMerMaskSize(uint32 mer_size){
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

gnSeqI SuffixArray::Length() const{
	return header.length;
}

sarID_t SuffixArray::GetID() const{
	return header.id;
}
void SuffixArray::SetID(const sarID_t d){
	header.id = d;
}

gnSeqI SuffixArray::bsearch(const struct bmer& query_mer, const gnSeqI start, const gnSeqI end) {

	gnSeqI middle = (start + end) / 2;
	struct bmer midmer = (*this)[middle];
	if(midmer.mer == query_mer.mer)
		return middle;
	else if((midmer.mer < query_mer.mer) && (middle < end))
		return bsearch(query_mer, middle + 1, end);
	else if((midmer.mer > query_mer.mer) && (start < middle))
		return bsearch(query_mer, start, middle - 1);
	
	//if we get here then the mer was not found.
	return GNSEQI_END;
}

//translate the character sequence to binary form based on the
//translation table.
void SuffixArray::translate(uint8* dest, const gnSeqC* src, const gnSeqI len) const{
	uint8 start_bit = 0;
	gnSeqI cur_byte = 0;
	uint32 alpha_bits = header.alphabet_bits;
	dest[cur_byte] = 0;
	for(uint32 i=0; i < len; i++){
		if(i==426269)
			cout << "asdf";
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

void SuffixArray::translate32(uint32* dest, const gnSeqC* src, const gnSeqI len) const{
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
			if(start_bit >= 32){
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
SarHeader SuffixArray::GetHeader() const{
	return header;
}

uint32 SuffixArray::UniqueMerCount(){
	if(header.unique_mers != NO_UNIQUE_COUNT)
		return header.unique_mers;

	uint32 MER_BUFFER_SIZE = 16384;  //not quite arbitrary (2^14)
	uint32 cur_pos = 0;
	vector<bmer> mer_vector;
	bmer prev_mer;
	uint32 m_unique = 0;
	while(cur_pos < header.length){
		if(!Read(mer_vector, MER_BUFFER_SIZE, cur_pos)){
			DebugMsg("SuffixArray::UniqueMerCount: Error reading bmer vector.");
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
void SuffixArray::ShiftWords(unsigned int* data, uint32 length, int32 bits){
	int32 word_bits = 8*sizeof(unsigned int);
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

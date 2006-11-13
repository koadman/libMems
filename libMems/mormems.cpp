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
#include "MimHash.h"
#include "MemScorer.h"
#include "SmallDnaSar.h"
#include <vector>

boolean LoadSequence(gnSequence& file_sequence, SmallDnaSar& file_sar, MemHash& mh, string& test_seq1, string& test_sar1, uint32 requested_mer_size){
	// Load the sequence and tell the user if it loaded successfully
	if(file_sequence.LoadSource(test_seq1))
	{
		cout << "Sequence loaded successfully.\n";
	}else{
		cout << "Error loading file.\n";
		return false;
	}
	cout << test_seq1 << " " << file_sequence.length() << " base pairs.\n";

	// define a SmallDnaSar to store a suffix array
	boolean success = file_sar.LoadFile(test_sar1);
	boolean recreate = false;
	if(success && (file_sar.MerSize() != requested_mer_size)){
		cout << "Mer size mismatch.  A new suffix array will be created.\n";
		recreate = true;
	}
	if(!success || recreate){
		int new_id;
		if(!recreate){
			cout << "Please designate a unique numeric id for this suffix array.\n";
			cin >> new_id;
		}else
			new_id = file_sar.GetID();

		cout << "Creating suffix array\n";
		long start_time = wxGetLocalTime();
		file_sar.Create(file_sequence, (sarID_t)new_id, requested_mer_size);
		long end_time = wxGetLocalTime();
 		cout << "Create time was: " << end_time - start_time << " seconds.\n";
	}else
		cout << "Suffix Array loaded successfully\n";
	mh.AddSuffixArray(&file_sar, &file_sequence);
	return true;
}

boolean PrintStats(vector<MemHashEntry>& mh){

}

boolean AlignSequences(MimHash& mh, uint32 requested_mer_size, uint32 requested_gap_tolerance){

	cout << "Now I'm going to try to create a mem hash.\n";
	long start_time = wxGetLocalTime();
	mh.Create(requested_mer_size);
	long end_time = wxGetLocalTime();
	cout << "Mem hash time was: " << end_time - start_time << " seconds.\n";
	cout << mh.MemCount() << " mems were hashed.\n";
	cout << mh.MemCollisionCount() << " mers collided with existing mems\n";

	vector<MemHashEntry> mem_list;
	start_time = wxGetLocalTime();
	mh.EliminateInclusions();
	end_time = wxGetLocalTime();
	cout << "Eliminate Inclusions time was: " << end_time - start_time << " seconds.\n";

	mh.GetMemList(mem_list);
	ofstream mem_out("mem_list.txt");
	for(uint32 i=0; i < mem_list.size(); i++){
		mem_out << mem_list[i].Length();
		for(uint32 j=0; j < mh.SequenceCount(); j++)
			mem_out << '\t' << mem_list[i][j];
		mem_out << '\n';
	}
	mem_out.close();
	vector<MimHashEntry> mim_list;
	mh.GetMimList(mim_list);

	ofstream distrib_out("distrib_list.txt");
	mh.PrintDistribution(distrib_out);
	distrib_out.close();
	
	start_time = wxGetLocalTime();
	mh.ResolveMismatches(requested_gap_tolerance);
	end_time = wxGetLocalTime();
	cout << "Inexact match gap jumping took " << end_time - start_time << " seconds.\n";
	mh.GetMimList(mim_list);
	cout << "There are " << mim_list.size() << " mims in the list.\n";
	ofstream mim_out("mim_list.txt");
	for(uint32 i=0; i < mim_list.size(); i++){
		mim_out << mim_list[i];
		mim_out << '\n';
	}
	mim_out.close();
	return true;
}
uint32 GetUserSeqCount(){
	uint32 my_seq_count;
	cout << "How many sequences would you like to align? ";
	cin >> my_seq_count;
	return my_seq_count;
}
int main( void )
{
	string chlamydia_seqs[4] = {"cpneuJ.gbk",
								"cpneuA.gbk",
								"cpneu.gbk",
								"ctra.gbk"};
	string chlamydia_sars[4] = {"cpneuJ.sar",
								"cpneuA.sar",
								"cpneu.sar",
								"ctra.sar"};
	string diarrhea_seqs[5] = {	"diarrhea\\ecolim52.fas",
								"diarrhea\\ec-japan.fas",
								"diarrhea\\EDL933.fas",
								"diarrhea\\o157sakai.fas",
								"diarrhea\\cft_v16f.fas" };
	string diarrhea_sars[5] = {	"diarrhea\\ecolim52.sar",
								"diarrhea\\ec-japan.sar",
								"diarrhea\\EDL933.sar",
								"diarrhea\\o157sakai.sar",
								"diarrhea\\cft_v16f.sar" };
	string mindy_seqs[11] = {"D:\\Development\\mormems\\bin\\locus44\\oz44-2DEC_7A.seq",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2DEC_5C.seq",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2DEC_6B.seq",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2DEC_9E.seq",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2ECOR_14.seq",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2ECOR_37.seq",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2ECOR_47.seq",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2EDL_933.seq",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2MG1655.seq",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2MN_297.seq",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2WI_6.seq"};

	string mindy_sars[11] = {"D:\\Development\\mormems\\bin\\locus44\\oz44-2DEC_7A.sar",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2DEC_5C.sar",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2DEC_6B.sar",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2DEC_9E.sar",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2ECOR_14.sar",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2ECOR_37.sar",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2ECOR_47.sar",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2EDL_933.sar",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2MG1655.sar",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2MN_297.sar",
						"D:\\Development\\mormems\\bin\\locus44\\oz44-2WI_6.sar"};
	
	// define a string to store the sequence file name
	string filename, menu_choice;
	// define a gnSequence to store the sequence
	vector<gnSequence> sequences;
	vector<SuffixArray*> suffixArrays;
	SmallDnaSar sar, sar2, sar3, sar4;
	MimHash mh;
	MemScorer ms;
	uint32 requested_table_size = DEFAULT_MEM_TABLE_SIZE;
	uint32 requested_mer_size = DNA_MER_SIZE;
	uint32 requested_gap_tolerance = DNA_MER_SIZE - 1;
	uint32 new_rep_tol;
	
	uint32 seq_count;
	while(true){
		cout << "\nAlignment tester main menu:\n";
		cout << "1) Align Mel's sequences\n";
		cout << "2) Align your favorite STD\n";
		cout << "3) Align diarrhea\n";
		cout << "4) Pick your own sequences to align\n";
		cout << "===========================\n";
		cout << "5) Set the mem hash table size  [" << requested_table_size << "]\n";
		cout << "6) Set the default mer size  [" << requested_mer_size << "]\n";
		cout << "7) Set the mer repeat tolerance  [" << mh.GetRepeatTolerance() << "]\n";
		cout << "8) Set the mer repeat enumeration tolerance  [" << mh.GetEnumerationTolerance() << "]\n";
		cout << "9) Set the inexact match gap jumping tolerance  [" << requested_gap_tolerance << "]\n";
		cin >> menu_choice;
		switch(menu_choice[0]){
			case '1':	//align mel's sequences
				seq_count = GetUserSeqCount();
				for(uint32 i=0; i < seq_count; i++){
					suffixArrays.push_back(new SmallDnaSar());
					sequences.push_back(gnSequence());
					if(!LoadSequence(sequences[i], *(SmallDnaSar*)suffixArrays[i], mh, mindy_seqs[i], mindy_sars[i], requested_mer_size))
						return -1;
					ms.AddSuffixArray(suffixArrays[i], &sequences[i]);
				}
				AlignSequences(mh, requested_mer_size, requested_gap_tolerance);
				break;
			case '2':
				seq_count = GetUserSeqCount();
				for(uint32 i=0; i < seq_count; i++){
					suffixArrays.push_back(new SmallDnaSar());
					sequences.push_back(gnSequence());
					if(!LoadSequence(sequences[i], *(SmallDnaSar*)suffixArrays[i], mh, chlamydia_seqs[i], chlamydia_sars[i], requested_mer_size))
						return -1;
					ms.AddSuffixArray(suffixArrays[i], &sequences[i]);
				}
				AlignSequences(mh, requested_mer_size, requested_gap_tolerance);
				break;
			case '3':
				seq_count = GetUserSeqCount();
				for(uint32 i=0; i < seq_count; i++){
					suffixArrays.push_back(new SmallDnaSar());
					sequences.push_back(gnSequence());
					if(!LoadSequence(sequences[i], *(SmallDnaSar*)suffixArrays[i], mh, diarrhea_seqs[i], diarrhea_sars[i], requested_mer_size))
						return -1;
					ms.AddSuffixArray(suffixArrays[i], &sequences[i]);
				}
				AlignSequences(mh, requested_mer_size, requested_gap_tolerance);
				break;
			case '4':	//choose your own seq files
				seq_count = GetUserSeqCount();
				for(uint32 i=0; i < seq_count; i++){
					// Get the name of the sequence to load
					string sar_filename;
					cout << "Enter the name of the sequence file to load:\n" ;
					cin >> filename;
					cout << "Enter the name of the suffix array file to load or create:\n" ;
					cin >> sar_filename;
					sequences.push_back(gnSequence());
					suffixArrays.push_back(new SmallDnaSar());
					if(!LoadSequence(sequences[sequences.size() - 1], *(SmallDnaSar*)suffixArrays[suffixArrays.size() - 1], mh, filename, sar_filename, requested_mer_size))
						return -1;
					ms.AddSuffixArray(suffixArrays[suffixArrays.size() - 1], &sequences[sequences.size() - 1]);
				}
				AlignSequences(mh, requested_mer_size, requested_gap_tolerance);
				break;
			case '5':
				cout << "How many mem hash table buckets would you like?  More buckets " <<
						"requires more memory but will run faster.\n";
				cin >> requested_table_size;
				if(requested_table_size < 5000)
					cout << "Warning! " << requested_table_size << " is very small!  Try 40000 or more.\n";
				mh = MimHash(requested_table_size);
				break;
			case '6':
				cout << "The mer size should be odd.\n";
				cout << "The current mer size is " << requested_mer_size << "\n";
				cout << "What would you like to change it to? ";
				cin >> requested_mer_size;
				cout << "Ok, new mer size is " << requested_mer_size << ".\n";
				break;
			case '7':
				cout << "The mer repeat tolerance defines a maximum number of times a mer may repeat ";
				cout << "itself within a single genome.  Mers which repeat themselves more than this ";
				cout << "will be thrown out entirely.\n";
				cout << "The current mer repeat tolerance is " << mh.GetRepeatTolerance() << ".\n";
				cout << "What would you like to change it to? ";
				cin >> new_rep_tol;
				mh.SetRepeatTolerance(new_rep_tol);
				cout << "Ok, new repeat tolerance is " << mh.GetRepeatTolerance() << ".\n";
				break;
			case '8':
				cout << "The mer repeat enumeration tolerance specifies the maximum number of mers from ";
				cout << "any one genome which used to enumerate all possible pairings of repeated mers";
				cout << "across genomes.\n";
				cout << "The current mer repeat enumeration tolerance is " << mh.GetEnumerationTolerance() << ".\n";
				cout << "What would you like to change it to? ";
				cin >> new_rep_tol;
				mh.SetEnumerationTolerance(new_rep_tol);
				cout << "Ok, new repeat tolerance is " << mh.GetEnumerationTolerance() << ".\n";
				break;
			case '9':
				cout << "The inexact match gap jumping tolerance specifies the maximum number of ";
				cout << "base pairs separating two matching regions in the same sequence which will ";
				cout << "be acceptable to link the two matching regions together.\n";
				cout << "The current gap tolerance is: " << requested_gap_tolerance << "\n";
				cout << "What would you like to change it to? ";
				cin >> requested_gap_tolerance;
				cout << "Ok, new repeat tolerance is " << requested_gap_tolerance << ".\n";
				break;
		}
	}
	cout << "\n";
	cout << "Scoring sar similarity heuristically..\n";
//	cout << "Heuristic hit score is: " << ms.HeuristicHitsScore() << "\n";

	cout << "Scoring sar similarity exactly..\n";
//	cout << "Exact hit score is: " << ms.PercentHitsScore() << "\n";

/*	gnSeqC *curBase = new gnSeqC[seq_count];
	for(uint32 i=0; i < mem_list.size(); i++){
		for(uint32 seqI = 0; seqI < seq_count; seqI++){
			curBase[seqI] = sequences[seqI][mem_list[i].start[seqI]];
			if(curBase[seqI]
		}
	}
	delete[] curBase;
*/	cin >> filename;
	return 0;
}
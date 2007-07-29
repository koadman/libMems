/*******************************************************************************
 * $Id: parallel_test.cpp,v 1.5 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

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
#include "libMems/MimHash.h"
#include "libMems/MemScorer.h"
#include "libMems/DNAFileSML.h"
#include <vector>
#include "libMems/RepeatHash.h"
#include "libMems/ParallelMemHash.h"

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

	// define a SmallDnaSar to store a sorted mer list
	boolean success = file_sar.LoadFile(test_sar1);
	boolean recreate = false;
	if(success && (file_sar.MerSize() != requested_mer_size)){
		cout << "Mer size mismatch.  A new sorted mer list will be created.\n";
		recreate = true;
	}
	if(!success || recreate){
		int new_id;
		if(!recreate){
			cout << "Please designate a unique numeric id for this sorted mer list.\n";
			cin >> new_id;
		}else
			new_id = file_sar.GetID();

		cout << "Creating sorted mer list\n";
		long start_time = wxGetLocalTime();
		file_sar.Create(file_sequence, (sarID_t)new_id, requested_mer_size);
		long end_time = wxGetLocalTime();
 		cout << "Create time was: " << end_time - start_time << " seconds.\n";
	}else
		cout << "Sorted mer list loaded successfully\n";
	mh.AddSequence(&file_sar, &file_sequence);
	return true;
}

boolean FindReps(RepeatHash& rephash){
	cout << "Now I'm going to try to create a mem hash.\n";
	long start_time = wxGetLocalTime();
	rephash.CreateMatches();
	long end_time = wxGetLocalTime();
	cout << "Repeat hash time was: " << end_time - start_time << " seconds.\n";
	cout << rephash.MemCount() << " mems were hashed.\n";
	cout << rephash.MemCollisionCount() << " mers collided with existing mems\n";

	vector<Match> mem_list;
	rephash.GetMemList(mem_list);
	ofstream mem_out("repeat_list.txt");
	for(uint32 i=0; i < mem_list.size(); i++){
		mem_out << mem_list[i];
		mem_out << '\n';
	}
	mem_out.close();
	return true;
}

boolean AlignSequences(ParallelMemHash& mh, uint32 requested_mer_size, uint32 requested_gap_tolerance){

	cout << "Now I'm going to try to create a mem hash.\n";
	long start_time = wxGetLocalTime();
	mh.CreateMems(3);
	long end_time = wxGetLocalTime();
	cout << "Mem hash time was: " << end_time - start_time << " seconds.\n";
	cout << mh.MemCount() << " mems were hashed.\n";
	cout << mh.MemCollisionCount() << " mers collided with existing mems\n";

	vector<Match> mem_list;
	start_time = wxGetLocalTime();
//	mh.EliminateInclusions();
	end_time = wxGetLocalTime();
	cout << "Eliminate Inclusions time was: " << end_time - start_time << " seconds.\n";

	mh.GetMemList(mem_list);
	ofstream mem_out("mem_list.txt");
	for(uint32 i=0; i < mem_list.size(); i++){
		mem_out << mem_list[i].Length();
		for(uint32 j=0; j < mem_list[i].SeqCount(); j++)
			mem_out << '\t' << mem_list[i][j];
		mem_out << '\n';
	}
	mem_out.close();
	vector<MimHashEntry> mim_list;

	ofstream distrib_out("distrib_list.txt");
	mh.PrintDistribution(distrib_out);
	distrib_out.close();

/*	mh.GetMimList(mim_list);
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
*/	return true;
}

uint32 GetUserSeqCount(){
	uint32 my_seq_count;
	cout << "How many sequences would you like to align? ";
	cin >> my_seq_count;
	return my_seq_count;
}


void Clear(MemHash& mh, vector<gnSequence>& sequences, vector<SuffixArray*>& suffixArrays){
	mh.Clear();
	sequences.clear();
	for(uint32 sarI = 0; sarI < suffixArrays.size(); sarI++)
		delete suffixArrays[sarI];
	suffixArrays.clear();
}

class DerivedApp : public wxApp
{
public:
  virtual bool OnInit();
};

IMPLEMENT_APP(DerivedApp)

bool DerivedApp::OnInit()
{

	string filename;
try{
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
	string menu_choice;
	// define a gnSequence to store the sequence
	vector<gnSequence> sequences;
	vector<SuffixArray*> suffixArrays;
	SmallDnaSar sar, sar2, sar3, sar4;
	MimHash mh;
	ParallelMemHash pmh;
	MemScorer ms;
	uint32 requested_table_size = DEFAULT_MEM_TABLE_SIZE;
	uint32 requested_mer_size = DNA_MER_SIZE;
	uint32 requested_gap_tolerance = DNA_MER_SIZE - 1;
	uint32 new_rep_tol;
	
	uint32 seq_count;

	seq_count = GetUserSeqCount();
	for(uint32 i=0; i < seq_count; i++){
		suffixArrays.push_back(new SmallDnaSar());
		sequences.push_back(gnSequence());
		if(!LoadSequence(sequences[i], *(SmallDnaSar*)suffixArrays[i], pmh, chlamydia_seqs[i], chlamydia_sars[i], requested_mer_size))
			return -1;
		ms.AddSequence(suffixArrays[i], &sequences[i]);
	}
//	ms.SubsequenceIdentityScore();
//	cout << "Identity Matrix: \n";
//	ms.PrintIDMatrix(cout);
//	cout << "Hit Matrix: \n";
//	ms.PrintHitMatrix(cout);
//	cout << "Detailed List: \n";
//	ms.PrintDetailList(cout);
	AlignSequences(pmh, requested_mer_size, requested_gap_tolerance);
	Clear(mh, sequences, suffixArrays);

}catch(gnException& e){
	cout << e;
}
	cin >> filename;
  ExitMainLoop();
  return TRUE;
}

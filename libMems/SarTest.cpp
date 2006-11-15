/*******************************************************************************
 * $Id: SarTest.cpp,v 1.8 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
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
#include "BigDiskSuffixArray.h"
#include "SmallDiskSuffixArray.h"
#include "SmallDnaSar.h"
#include "BigDnaSar.h"
#include "libMems/MemHash.h"
#include "MemorySuffixArray.h"

int main( void )
{
/*	string filename, output_filename;
	ifstream mem_file;
	ofstream out_file;
	uint32 count = 0;
	while(!mem_file.is_open()){
		cout << "Enter the name of the mem_list to load:\n" ;
		cin >> filename;
		mem_file.open(filename.c_str());
	}
	while(!out_file.is_open()){
		cout << "Enter the name of the output file:\n" ;
		cin >> output_filename;
		out_file.open(output_filename.c_str());
	}
	while(mem_file.good()){
		int32 len, g1, g2, g3, g4, g5;
		mem_file >> len;
		mem_file >> g1;
		mem_file >> g2;
		mem_file >> g3;
		mem_file >> g4;
		mem_file >> g5;
		if(g1 != 0 && g2 != 0 && g1 > 3000000){
			out_file << len << '\t';
			out_file << g1 << '\t';
			out_file << g2 << '\t';
			out_file << g3 << '\t';
			out_file << g4 << '\t';
			out_file << g5;
			out_file << "\n";
			count++;
		}
	}
	cout << "Found " << count << " mems\n";
	cin >> filename;
*/
	// define a string to store the sequence file name
	string filename;
try{

	string test_seq1 = "cpneuJ.gbk";
	string test_sar1 = "cpneuJsmall.sar";
	string test_sar2 = "cpneuJbig.fsar";
	boolean user_input = false;

	// define a gnSequence to store the sequence
	gnSequence file_sequence, file_sequence_rc;
	
	// Get the name of the sequence to load
	cout << "Enter the name of the sequence file to load:\n" ;
	if(user_input)
		cin >> test_seq1;
	else
		cout << test_seq1 << "\n";
	// Load the sequence and tell the user if it loaded successfully
	if(file_sequence.LoadSource(test_seq1))
	{
		cout << "Sequence loaded successfully.\n";
	}else{
		cout << "Error loading file.\n";
		return -1;
	}
	cout << test_seq1 << " " << file_sequence.length() << " base pairs.\n";

	cout << "Enter the name of the suffix array to create:\n" ;
	if(user_input)
		cin >> test_sar1;
	else
		cout << test_sar1 << "\n";

	// define a BigDiskSuffixArray to store a suffix array
	SmallDnaSar file_sar(test_sar1);
	long start_time = wxGetLocalTime();
//	file_sar.LoadFile(test_sar1);
	file_sar.Create(file_sequence, 2);
	long end_time = wxGetLocalTime();
	cout << "Create time was: " << end_time - start_time << " seconds.\n";
	file_sar.SetDescription("This is a simple suffix array.");

	MemorySuffixArray file_sar2;
	start_time = wxGetLocalTime();
//	file_sar2.LoadFile(test_sar2);
	file_sar2.Create(file_sequence, 2);
	end_time = wxGetLocalTime();
	cout << "Create time was: " << end_time - start_time << " seconds.\n";
	file_sar2.SetDescription("This is a simple suffix array.");

/*	BigDnaSar file_sar2;
	start_time = wxGetLocalTime();
	file_sar2.LoadFile(test_sar2);
//	file_sar2.Create(file_sequence, 2);
	end_time = wxGetLocalTime();
	cout << "Create time was: " << end_time - start_time << " seconds.\n";
	file_sar2.SetDescription("This is a simple suffix array.");
*/
	if(file_sequence.length() != file_sar.Length())
		cout << "Porknobbin hobjob\n";
	else
		cout << "asdf\n";
		
	uint32 lennard = file_sar2.Length();
	uint32 i = 0;
	for(; i < lennard; i++){
		bmer a = file_sar[i];
		bmer b = file_sar2[i];
		if(a.mer != b.mer ){
//		if(a.mer != b.mer || a.position != b.position){
			cout << "Suffix arrays are not equal";
			break;
		}
	}
	if(i == lennard)
		cout << "Suffix arrays are equal under god.\n";
}catch( gnException& gne ){
	cout << gne;
}
	
	cin >> filename;
	return 0;
}

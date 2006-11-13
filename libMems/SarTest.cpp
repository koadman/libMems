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
#include "MemHash.h"

int main( void )
{
	string test_seq1 = "cpneuJ.fas";
	string test_sar1 = "cpneuJsmall.fsar";
	string test_sar2 = "cpneuJbig.fsar";
	boolean user_input = false;

	// define a string to store the sequence file name
	string filename;
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
	file_sar.BigCreate(file_sequence, 2);
	long end_time = wxGetLocalTime();
	cout << "Create time was: " << end_time - start_time << " seconds.\n";
	file_sar.SetDescription("This is a simple suffix array.");

	BigDnaSar file_sar2(test_sar2);
	start_time = wxGetLocalTime();
	file_sar2.Create(file_sequence, 2);
	end_time = wxGetLocalTime();
	cout << "Create time was: " << end_time - start_time << " seconds.\n";
	file_sar2.SetDescription("This is a simple suffix array.");
	
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
	cin >> filename;
	return 0;
}
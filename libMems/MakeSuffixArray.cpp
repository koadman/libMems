/*******************************************************************************
 * $Id: MakeSuffixArray.cpp,v 1.6 2004/03/01 02:40:08 darling Exp $
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

#include "libGenome/gnSequence.h"
#include "libMems/DNAFileSML.h"

void print_usage(char* pname){
	cout << "Usage: " << pname << " <sequence file> <sar file> <mer size> [split count]";
}

int main(int argc, char* argv[]){
	
	if(argc != 4 && argc != 5){
		print_usage(argv[0]);
		return -1;
	}
	
	string seq_filename(argv[1]);
	string sar_filename(argv[2]);
	
	gnSequence file_sequence;
	DNAFileSML file_sar;
	uint32 mer_size = atoi(argv[3]);
	uint32 split_size = 0;
	if(argc == 5)
		split_size = atoi(argv[4]);

	// Load the sequence and tell the user if it loaded successfully
	if(file_sequence.LoadSource(seq_filename))
	{
		cout << "Sequence loaded successfully.\n";
	}else{
		cout << "Error loading file.\n";
		return false;
	}
	cout << seq_filename << " " << file_sequence.length() << " base pairs.\n";

	// define a SmallDnaSar to store a sorted mer list
	boolean success = true;
	try{
		file_sar.LoadFile(sar_filename);
	}catch(...){
		success = false;
	}
	boolean recreate = false;
	if(success && (file_sar.SeedLength() != mer_size)){
		cout << "Mer size mismatch.  A new sorted mer list will be created.\n";
		recreate = true;
	}
	if(!success || recreate){
		cout << "Creating sorted mer list\n";
		long start_time = wxGetLocalTime();
		file_sar.Create(file_sequence, mer_size);
		long end_time = wxGetLocalTime();
 		cout << "Create time was: " << end_time - start_time << " seconds.\n";
	}else
		cout << "Sorted mer list loaded successfully\n";

}

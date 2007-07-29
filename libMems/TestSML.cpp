/*******************************************************************************
 * $Id: TestSML.cpp,v 1.6 2004/03/01 02:40:08 darling Exp $
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

#include "wx/timer.h"
#include "libGenome/gnSequence.h"
#include "libMems/DNAFileSML.h"
#include "libMems/DNAMemorySML.h"
#include "libMems/MemHash.h"

int main( int argc, char* argv[] )
{
	// define a string to store the sequence file name
	string filename;
	try{
		
		if( argc != 5 ){
			cout << "Usage: TestSML <seq1> <sml1> <seq2> <sml2>\n";
			return -1;
		}
		string test_seq1 = argv[1];
		string test_sar1 = argv[2];
		string test_seq2 = argv[3];
		string test_sar2 = argv[4];

		// define a gnSequence to store the sequence
		gnSequence file_sequence, file_sequence_rc;
		
		// Load the sequence and tell the user if it loaded successfully
/*		if(file_sequence.LoadSource(test_seq1))
		{
			cout << "Sequence loaded successfully.\n";
		}else{
			cout << "Error loading file.\n";
			return -1;
		}
		cout << test_seq1 << " " << file_sequence.length() << " base pairs.\n";
*/
		// define a DNAFileSML to store a sorted mer list
		DNAFileSML file_sar(test_sar1);
		long start_time = wxGetLocalTime();
		file_sar.LoadFile(test_sar1);
//		file_sar.Create(file_sequence, 31);
		long end_time = wxGetLocalTime();
		cout << "Load time was: " << end_time - start_time << " seconds.\n";
		file_sar.SetDescription("This is a simple sorted mer list.");

//		DNAMemorySML file_sar2;
		DNAFileSML file_sar2(test_sar2);
		start_time = wxGetLocalTime();
		file_sar2.LoadFile(test_sar2);
//		file_sar2.Create(file_sequence);
		end_time = wxGetLocalTime();
		cout << "Load time was: " << end_time - start_time << " seconds.\n";
		file_sar2.SetDescription("This is a simple sorted mer list.");
		
		SMLHeader header = file_sar.GetHeader();
		cout << "Header for " << test_sar1 << endl;
		cout << "Format Version: " << header.version << endl;
		cout << "alphabet bits: " << header.alphabet_bits << endl;
		cout << "seed weight: " << header.seed_weight << endl;
		cout << "length: " << header.length << endl;
		cout << "unique_mers: " << header.unique_mers << endl;
		cout << "word_size: " << header.word_size << endl;
		cout << "little endian: " << (int)header.little_endian << endl;
		cout << "id: " << (int)header.id << endl;
		cout << "circular: " << (int)header.circular << endl;
		cout << endl;

                header = file_sar2.GetHeader();
                cout << "Header for " << test_sar2 << endl;
                cout << "Format Version: " << header.version << endl;
                cout << "alphabet bits: " << header.alphabet_bits << endl;
                cout << "seed weight: " << header.seed_weight << endl;
                cout << "length: " << header.length << endl;
                cout << "unique_mers: " << header.unique_mers << endl;
                cout << "word_size: " << header.word_size << endl;
                cout << "little endian: " << (int)header.little_endian << endl;
                cout << "id: " << (int)header.id << endl;
                cout << "circular: " << (int)header.circular << endl;

//		if(file_sequence.length() != file_sar.Length())
//			cout << "";
//		else
//			cout << "asdf\n";
			
		uint32 lennard = file_sar2.Length();
//		if(!file_sequence.isCircular())
//			lennard -= file_sar2.MerSize();
		uint32 i = 0;
		uint64 mer_mask = file_sar.GetMerMask();
		bmer a, b, c;
		boolean a_decrease = false;
		boolean b_decrease = false;
		for(i = 0; i < lennard; i++){
			b = file_sar2[i];
			if( i > 0 && b.mer & mer_mask < c.mer & mer_mask ){
				cout << i - 1 << '\t' << i << " are decreasing\n";
				b_decrease = true;
			}
			c = b;
		}
		if( !b_decrease ){
			cout << test_sar2 << " is non-decreasing\n";
		}
		for(i = 0; i < lennard; i++){
			a = file_sar[i];
			b = file_sar2[i];
			if( (a.mer & mer_mask) != (b.mer & mer_mask) ){
	//		if(a.mer != b.mer || a.position != b.position){
				cout << "Sorted Mer Lists are not equal";
				break;
			}
		}
                if(i == lennard)
                        cout << "Sorted Mer Lists are equal under god.\n";

		for(i = 0; i < lennard; i++){
			a = file_sar[i];
                        if( i > 0 && a.mer & mer_mask < c.mer & mer_mask ){
                                cout << i - 1 << '\t' << i << " are decreasing\n";
                                a_decrease = true;
                        }
			c = a;
		}
                if( !a_decrease ){
                        cout << test_sar1 << " is non-decreasing\n";
                }

	}catch( gnException& gne ){
		cout << gne;
	}
	
	return 0;
}


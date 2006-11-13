#ifndef _MatchList_h_
#define _MatchList_h_
#include <iostream>
#include <list>
#include "SortedMerList.h"
#include "gn/gnSequence.h"
#include "Match.h"
#include "MemList.h"

class MatchList {
public:
	MatchList(){};
	MatchList( const MatchList& ml );
	MatchList& operator=( const MatchList& ml );
	
	/**
	 * Attempts to load the sequences and smls designated by the
	 * elements of the sml_filename and seq_filename vectors.  This
	 * method will create the sorted mer lists if they do not exist.
	 * The gnSequence and SortedMerList objects are created on the heap
	 * and are not deallocated when this class is destroyed.  They should
	 * be manually destroyed when no longer in use.
	 * @throws InvalidArgument will be thrown if the number of entries in seq_filename does not match sml_filename
	 */
	void LoadSequences( uint mer_size, ostream* log_stream );

	/**
	 * Reads a MatchList from an input stream
	 * Sequence and SML file names are read into the seq_filename
	 * and sml_filename vectors, but the actual files are not
	 * opened.  The calling function should load them after
	 * using this method.
	 * @param match_stream The input stream to read from
	 */
	void ReadList( istream& match_stream );

	/**
	 *  Writes a MatchList to the designated output stream
	 * @param match_stream The outptu stream to write to
	 */
	void WriteList( ostream& match_stream ) const;
	
	void MultiplicityFilter( unsigned mult );
	void ExactFilter( valarray<bool>& filter_spec );
	void IntersectFilter( valarray<bool>& filter_spec );
	/**
	 * Keeps only matches that are not subsets and whose first matching
	 * sequence is the specified sequence.
	 */
	void UnlinkedFirstStartFilter( unsigned start_seq );
	
	void FromMemList( const list<MemHashEntry*>& mem_list );
	void ToMemList( list<MemHashEntry*>& mem_list ) const;

	list<Match*> match_list;
	vector<string> sml_filename;
	vector<string> seq_filename;
	vector<SortedMerList*> sml_table;
	vector<gnSequence*> seq_table;

protected:

};

CREATE_EXCEPTION( InvalidArgument );


#endif	//_MatchList_h_

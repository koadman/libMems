/*******************************************************************************
 * $Id: MatchList.h,v 1.10 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifndef _MatchList_h_
#define _MatchList_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <list>
#include "libMems/SortedMerList.h"
#include "libGenome/gnSequence.h"
#include "libMems/Match.h"

namespace mems {

class MatchList : public std::vector< Match* > {
public:
	MatchList(){};
	MatchList( const MatchList& ml );
	MatchList& operator=( const MatchList& ml );

	/**
	 * Attempts to load the sequences designated by the
	 * elements of the seq_filename vector.
	 * The genome::gnSequence objects are created on the heap
	 * and are not deallocated when this class is destroyed.  They should
	 * be manually destroyed when no longer in use.
	 */
	void LoadSequences( std::ostream* log_stream );

	/**
	 * Attempts to load SMLs designated by the
	 * elements of the sml_filename vector.  This
	 * method will create the sorted mer lists if they do not exist.
	 * The DNAFileSML objects are created on the heap
	 * and are not deallocated when this class is destroyed.  They should
	 * be manually destroyed when no longer in use.
	 * @param	seed_rank	The rank of the seed to use, 0-2 are ranked spaced seeds, 
	 *						other options include CODING_SEED and SOLID_SEED
	 */
	void LoadSMLs( uint mer_size, std::ostream* log_stream, int seed_rank = 0 );

	/**
	 * Loads sequences to align from a Multi-FastA file and constructs a SML
	 * for each sequence entry in the file.
	 * The genome::gnSequence and SortedMerList objects are created on the heap
	 * and are not deallocated when this class is destroyed.  They should
	 * be manually destroyed when no longer in use.
	 *
	 * @param mfa_filename	The name of the Multi-FastA file to read in.  Each 
	 *						sequence entry will be treated as a separate sequence to 
	 *						be aligned.
	 * @param mer_size		The seed size to use when constructing the sorted mer lists
	 * @param log_stream	An output stream to log messages to.  If NULL no logging is done
	 * @param load_smls		Specifies whether sorted mer lists should be created 
	 * 						for each sequence entry
	 */
	void LoadMFASequences( const std::string& mfa_filename, uint mer_size, std::ostream* log_stream, boolean load_smls = true, int seed_rank = 0 );

	/**
	 * Calculates a default search mer size for the given set of sequences
	 * @param seq_table		The vector of sequences to calculate a default mer size for
	 */
	static uint GetDefaultMerSize( const std::vector< genome::gnSequence* >& seq_table );
	
	/**
	 * Deletes the genome::gnSequence, SortedMerList, and Match objects associated
	 * with this MatchList.
	 */
	void Clear();

	/**
	 * Reads a MatchList from an input stream
	 * Sequence and SML file names are read into the seq_filename
	 * and sml_filename vectors, but the actual files are not
	 * opened.  The calling function should load them after
	 * using this method.
	 * @param match_stream The input stream to read from
	 */
	void ReadList( std::istream& match_stream );

	/**
	 *  Writes a MatchList to the designated output stream
	 * @param match_stream The output stream to write to
	 */
	void WriteList( std::ostream& match_stream ) const;
	
	/**
	 * Removes all matches that have a multiplicity lower than the specified level
	 * @param mult	The multiplicity filter threshold
	 */
	void MultiplicityFilter( unsigned mult );

	/**
	 * Removes all matches that shorter than the specified length
	 * @param length	The minimum length
	 */
	void LengthFilter( gnSeqI length );

	/**
	 * Removes matches that do not match in exactly the sequences specified in filter_spec
	 * @param filter_spec 	The specification of the exact filter, true designates that the
	 *						match must exist in that sequence.  filter_spec must contain
	 *						one boolean entry for every sequence.
	 */
//	void ExactFilter( valarray< bool >& filter_spec );
	/**
	 * Removes matches that do not intersect with the sequences specified in filter_spec
	 * @param filter_spec 	The specification of the intersection filter, true designates
	 *						match must exist in that sequence.  filter_spec must contain
	 *						one boolean entry for every sequence.
	 */
//	void IntersectFilter( valarray< bool >& filter_spec );

	/**
	 * Keeps only matches that are not subsets and whose first matching
	 * sequence is the specified sequence.
	 */
	void UnlinkedFirstStartFilter( unsigned start_seq );
	
	std::vector<std::string> sml_filename;		/**< The file names of the sorted mer list for each sequence, may be empty or null */
	std::vector<std::string> seq_filename;		/**< The file names of the sequence data, may be empty or null */
	std::vector<SortedMerList*> sml_table;	/**< The sorted mer list associated with each sequence, may be empty or null */
	std::vector<genome::gnSequence*> seq_table;		/**< The actual sequences associated with the matches stored in this list.  Should not be empty or null. */

protected:

};

CREATE_EXCEPTION( InvalidArgument );

}

#endif	//_MatchList_h_

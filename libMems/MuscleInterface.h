/*******************************************************************************
 * $Id: MuscleInterface.h,v 1.12 2004/04/19 23:10:50 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifndef _MuscleInterface_h_
#define _MuscleInterface_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/NumericMatrix.h"
#include "libGenome/gnFilter.h"
#include "libGenome/gnSequence.h"
#include "libMems/GappedAlignment.h"
#include "libMems/GappedAligner.h"

namespace mems {

extern bool debug_muscle;

class MuscleInterface : public GappedAligner {
public:
	~MuscleInterface()
	{
		ClearCommandLine();
	}
	/**
	 * Returns a reference to a usable MuscleInterface
	 */
	static MuscleInterface& getMuscleInterface();

	/**
	 * Parse the execution path from argv[0] and set the muscle
	 * path accordingly
	 */
	void ParseMusclePath( const char* argv0 );

	/** 
	 * Set the path to the muscle executable
	 * Defaults to "muscle"
	 */
	void SetMusclePath( const std::string& path );

	/** 
	 * Set the arguments to use when executing muscle 
	 * Defaults to -stable -quiet
	 */
	void SetMuscleArguments( const std::string& args );

	/**
	 * Attempts to perform a multiple alignment using Muscle between
	 * <code>r_begin</code> and <code>r_end</code>
	 */
	boolean Align( GappedAlignment& cr, Match* r_begin, Match* r_end, std::vector< genome::gnSequence* >& seq_table );

	bool Refine( GappedAlignment& ga, size_t windowsize = 0 );

	/**
	 * Given two gapped alignments in ga1 and ga2, align them and store the result in aln.  ga1 and
	 * ga2 must have equal sequence count and contain disjoint sets of sequences, e.g. for any given
	 * seqI, if ga1.LeftEnd(seqI) != NO_MATCH, then ga2.LeftEnd(seqI) == NO_MATCH 
	 */
	bool ProfileAlign( const GappedAlignment& ga1, const GappedAlignment& ga2, GappedAlignment& aln, bool anchored = true );

protected:
	std::string muscle_path;
	std::string muscle_arguments;
	char** muscle_cmdline;
	boolean CallMuscle( std::vector< std::string >& aln_matrix, const std::vector< std::string >& seq_table );
	void ClearCommandLine()
	{
		if( muscle_cmdline != NULL )
		{
			size_t cmdI = 0;
			while(muscle_cmdline[cmdI] != NULL)
			{
				delete[] muscle_cmdline[cmdI];
				cmdI++;
			}
			delete[] muscle_cmdline;
		}
	}

private:
	MuscleInterface( const MuscleInterface& ci ){ *this = ci; }
	MuscleInterface& operator=( const MuscleInterface& ci );
	MuscleInterface();
};

}

#endif // _MuscleInterface_h_

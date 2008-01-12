/*******************************************************************************
 * $Id: ParallelMemHash.cpp,v 1.32 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/ParallelMemHash.h"
#include <vector>

#ifdef _OPENMP


using namespace std;
using namespace genome;
namespace mems {

	ParallelMemHash::ParallelMemHash() : MemHash()

{
}

ParallelMemHash::ParallelMemHash(const ParallelMemHash& mh) : MemHash(mh)
{
	*this = mh;
}

ParallelMemHash& ParallelMemHash::operator=( const ParallelMemHash& mh ){
	thread_mem_table = mh.thread_mem_table;
	return *this;
}

ParallelMemHash* ParallelMemHash::Clone() const{
	return new ParallelMemHash(*this);
}

void ParallelMemHash::FindMatches( MatchList& ml ) 
{
	for( uint32 seqI = 0; seqI < ml.seq_table.size(); ++seqI ){
                if( !AddSequence( ml.sml_table[ seqI ], ml.seq_table[ seqI ] ) ){
                        ErrorMsg( "Error adding " + ml.seq_filename[seqI] + "\n");
                        return;
                }
        }

	size_t CHUNK_SIZE = 500000;
	// break up the SMLs into nice small chunks
	vector< vector< gnSeqI > > chunk_starts;
	vector< gnSeqI > chunk_lengths;

	// break up on the longest SML
	int max_length_sml = -1;
	size_t maxlen = 0;
	for( size_t i = 0; i < ml.sml_table.size(); i++ )
		if( ml.sml_table[i]->Length() > maxlen )
		{
			maxlen = ml.sml_table[i]->Length();
			max_length_sml = i;
		}

	chunk_starts.push_back( vector< gnSeqI >( seq_count, 0 ) );

	while( chunk_starts.back()[max_length_sml] + CHUNK_SIZE < ml.sml_table[max_length_sml]->Length() )
	{
		vector< gnSeqI > tmp( seq_count, 0 );
		GetBreakpoint(max_length_sml, chunk_starts.back()[max_length_sml] + CHUNK_SIZE, tmp);
		chunk_starts.push_back(tmp);
	}
	
	// now that it's all chunky, search in parallel
#pragma omp parallel for schedule(dynamic)
	for( int i = 0; i < chunk_starts.size(); i++ )
	{
		if(thread_mem_table.get().size() != mem_table.size())
			thread_mem_table.get().resize( mem_table.size() );

		vector< gnSeqI > chunk_lens(seq_count);
		if( i + 1 < chunk_starts.size() )
		{
			for( size_t j = 0; j < seq_count; j++ )
				chunk_lens[j] = chunk_starts[i+1][j] - chunk_starts[i][j];
		}else
			chunk_lens = vector< gnSeqI >( seq_count, GNSEQI_END );
		SearchRange( chunk_starts[i], chunk_lens );
		MergeTable();
	}
	GetMatchList( ml );	
}

void ParallelMemHash::MergeTable()
{
#pragma omp critical
	{
		size_t buckets = thread_mem_table.get().size();
		for( size_t bI = 0; bI < buckets; bI++ )
		{
			vector< MatchHashEntry* >& bucket = thread_mem_table.get()[bI];
			for( size_t mI = 0; mI < bucket.size(); mI++ )
				AddHashEntry(bucket[mI], mem_table);
		}
		thread_mem_table.get() = mem_table;
	}
}

/** This is a class for guard objects using OpenMP
*  It is adapted from the book
*  "Pattern-Oriented Software Architecture". */
class omp_guard {
public:
    /** Acquire the lock and store a pointer to it */
	omp_guard (omp_lock_t &lock){
      acquire ();
	}

	/** Set the lock explicitly */
	void acquire (){
		omp_set_lock (lock_);
		owner_ = true;
	}

	/** Release the lock explicitly (owner thread only!) */
	void release (){
		if (owner_) {
			owner_ = false;
			omp_unset_lock (lock_);
		};
	}

	/** Destruct guard object */
	~omp_guard (){ release (); }
 
private:
    omp_lock_t *lock_;  // pointer to our lock
    bool owner_;   // is this object the owner of the lock?
   
    // Disallow copies or assignment
    omp_guard (const omp_guard &);
    void operator= (const omp_guard &);
};


MatchHashEntry* ParallelMemHash::AddHashEntry(MatchHashEntry& mhe){
	// do the normal procedure, but use the thread-local mem table.
	return MemHash::AddHashEntry(mhe, thread_mem_table.get());
}


} // namespace mems

#endif // _OPENMP

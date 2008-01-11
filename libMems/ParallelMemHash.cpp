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
	omp_locks.resize(table_size);
}

ParallelMemHash::ParallelMemHash(const ParallelMemHash& mh) : MemHash(mh), allocator( SlotAllocator<MatchHashEntry>::GetSlotAllocator() )
{
	*this = mh;
}

ParallelMemHash& ParallelMemHash::operator=( const ParallelMemHash& mh ){
	omp_locks = mh.omp_locks;
	return *this;
}

ParallelMemHash* ParallelMemHash::Clone() const{
	return new ParallelMemHash(*this);
}


void ParallelMemHash::SetTableSize(uint32 new_table_size)
{
	MemHash::SetTableSize(new_table_size);
	omp_locks.resize(new_table_size);
}

void ParallelMemHash::FindMatches( MatchList& ml ) 
{
	size_t CHUNK_SIZE = 500000;
	// break up the SMLs into nice small chunks
	vector< vector< gnSeqI > > chunk_starts;
	vector< gnSeqI > chunk_lengths;

	// break up on the longest SML
	int max_length_sml = -1;
	size_t maxlen = 0;
	for( size_t i = 0; i < ml.sml_table.size(); i++ )
		if( ml.sml_table[i].Length() > maxlen )
		{
			maxlen = ml.sml_table[i].Length();
			max_length_sml = i;
		}

	chunk_starts.push_back( vector< gnSeqI >( seq_count, 0 ) );

	size_t cur_len = 0;
	while( cur_len < ml.sml_table[max_length_sml].Length() )
	{
		vector< gnSeqI > tmp( seq_count, 0 );
		GetBreakpoint(max_length_sml, chunk_starts.back()[max_length_sml] + break_size, tmp);
		chunk_starts.push_back(tmp);
	}
	
	// now that it's all chunky, search in parallel
#pragma omp parallel for shedule(dynamic)
	for( int i = 0; i < chunk_starts.size()-1; i++ )
	{
		vector< gnSeqI > chunk_lens;
		for( size_t j = 0; j < seq_count; j++ )
			chunk_lens[j] = chunk_starts[i+1] - chunk_starts[i];
		SearchRange( chunk_starts, chunk_lens );
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


// Tries to add a new mem to the mem hash table
// If the mem already exists in the table, a pointer to it
// is returned.  Otherwise mhe is added and a pointer to
// it is returned.
MatchHashEntry* ParallelMemHash::AddHashEntry(MatchHashEntry& mhe){
	//first compute which hash table bucket this is going into
	int64 offset = mhe.Offset();

	uint32 bucketI = ((offset % table_size) + table_size) % table_size;

	// lock the bucket
	omp_guard(omp_locks[bucketI]);

	// do the normal procedure
	return MemHash::AddHashEntry(mhe);
}


} // namespace mems

#endif // _OPENMP

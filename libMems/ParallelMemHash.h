#ifndef _ParallelMemHash_h_
#define _ParallelMemHash_h_

#include "MemHash.h"

/**
 * The default amount of sequence to preprocess for matches before spawning threads
 * to continue searching.  Preprocessing is done to avoid heavy lock contention among
 * threads.  This should be a number between 0 and 1. 
 */
const float DEFAULT_PREPROCESS_PERCENT = .01; 

/**
 * Implements a parallelized version of exact match extraction
 */
class ParallelMemHash : public MemHash{
public:
	ParallelMemHash();
	~ParallelMemHash();
	ParallelMemHash(const ParallelMemHash& mh);
	virtual ParallelMemHash* Clone() const;
	virtual void Clear();
	/**
	 * Search for matches using the specified number of threads.
	 * @param thread_count The number of threads to spawn.  
	 * Generally this should be set to the number of processors in the system.
	 */
	virtual void CreateMems(uint32 thread_count);

protected:
	float m_preprocess_percent;
};


#endif //_ParallelMemHash_h_

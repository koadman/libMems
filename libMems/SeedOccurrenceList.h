#ifndef __SeedOccurrenceList_h__
#define __SeedOccurrenceList_h__

#include <vector>
#include "libMems/SortedMerList.h"

namespace mems
{

class SeedOccurrenceList
{
public:
	typedef float32 frequency_type;

	SeedOccurrenceList(){}

	template< typename SMLType >
	void construct( SMLType& sml )
	{
		std::vector<mems::bmer> mer_vec;
		sml.Read( mer_vec, sml.SMLLength(), 0 );
		count.resize( sml.Length() );
		size_t seed_start = 0;
		size_t cur_seed_count = 1;
		uint64 mer_mask = sml.GetSeedMask();
		size_t seedI = 1;
		for( ; seedI < mer_vec.size(); ++seedI )
		{
			if( (mer_vec[seedI].mer & mer_mask) == (mer_vec[seedI-1].mer & mer_mask) )
			{
				++cur_seed_count;
				continue;
			}
			// set seed frequencies
			for( size_t i = seed_start; i < seedI; ++i )
				count[mer_vec[i].position] = (frequency_type)cur_seed_count;
			seed_start = seedI;
			cur_seed_count = 1;
		}
		// set seed frequencies for the last few
		for( size_t i = seed_start; i < seedI && i < mer_vec.size(); ++i )
			count[mer_vec[i].position] = (frequency_type)cur_seed_count;
		// hack: fudge the last few values on the end of the sequence
		for( ; seedI < count.size(); ++seedI )
			count[seedI]=1;

		smoothFrequencies( sml );

		// wipe out any stray zeros
		for( size_t i = 0; i < count.size(); ++i )
			if( count[i]== 0 )
				count[i] = 1;
	}


	frequency_type getFrequency( gnSeqI position )
	{
		return count[position];
	}

protected:
	/**
	 * converts position freqs to the average freq of all k-mers containing that position
	 */
	template< typename SMLType >
	void smoothFrequencies( const SMLType& sml )
	{
		size_t seed_length = sml.SeedLength();
		// hack: for beginning (seed_length) positions assume that previous
		// containing seeds were unique
		double sum = seed_length - 1 + count[0];
		std::vector<frequency_type> buf(seed_length, 1);
		buf[0] = count[0];
		for( size_t i = 1; i < count.size(); i++ )
		{
			count[i-1] = sum / seed_length;
			sum += count[i];
			size_t bufI = i % seed_length;
			sum -= buf[bufI];
			buf[bufI] = count[i];
		}
	}


	std::vector<frequency_type> count;	
};

}

#endif	// __SeedOccurrenceList_h__


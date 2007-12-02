/*******************************************************************************
 * $Id: GenericMatch.h,v 1.10 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _Match_h_
#define _Match_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnClone.h"
#include <iostream>
#include <set>
#include "libMems/UngappedLocalAlignment.h"
#include "libMems/SparseAbstractMatch.h"
#include "libMems/DenseAbstractMatch.h"
#include "libMems/HybridAbstractMatch.h"

namespace mems {

/**
 * The GenericMatch class stores the location of an <b>equal size</b> (inexact or exactly) 
 * matching region
 * between several sequences.  There are numerous functions in this
 * class which can be used to compare and manipulate this match.
 */
template< class BaseType >
class GenericMatch : public BaseType 
{

public:
	GenericMatch(){};
	/**
	 * Creates a new GenericMatch.
	 * @param seq_count The total number of sequences in the alignment
	 * @param mersize The size of the mers used in the sorted mer lists.
	 * @param m_type The type of mem to create, can either be a seed or already extended.
	 * @see MemType
	 */
	GenericMatch( const uint seq_count );

	// use compiler generated copy constructor, assignment operator, and destructor

	GenericMatch* Clone() const;
	GenericMatch* Copy() const;
	virtual void Free();
	
	/** comparison operator, compares two matches to see if they are the same */
	boolean operator==(const GenericMatch<BaseType>& mhe) const;
			
protected:
	void HomogenizeParity(std::vector<int64>& start1, std::vector<int64>& start2, const uint32 startI) const;

};

/** define the Match class as a specialized version of GenericMatch for backwards API compatibility */
//typedef GenericMatch< UngappedLocalAlignment< DenseAbstractMatch<8> > > Match;
//typedef GenericMatch< UngappedLocalAlignment< SparseAbstractMatch< boost::pool_allocator<gnSeqI>, boost::pool_allocator<uint> > > > Match;
//typedef GenericMatch< UngappedLocalAlignment< SparseAbstractMatch<> > > Match;
typedef GenericMatch< UngappedLocalAlignment< HybridAbstractMatch<> > > Match;

template< class BaseType >
GenericMatch<BaseType>* GenericMatch<BaseType>::Copy() const
{
	return m_allocateAndCopy( *this );
}
template< class BaseType >
void GenericMatch<BaseType>::Free()
{
	m_free(this);
}

template< class BaseType >
GenericMatch<BaseType>::GenericMatch(uint seq_count)
 : BaseType( seq_count )
{}


template< class BaseType >
GenericMatch<BaseType>* GenericMatch<BaseType>::Clone() const
{
	return new GenericMatch(*this);
}

template< class BaseType >
boolean GenericMatch<BaseType>::operator==(const GenericMatch<BaseType>& mhe) const
{
	return BaseType::operator ==(mhe);
}


template< class BaseType >
void GenericMatch<BaseType>::HomogenizeParity(std::vector<int64>& start1, std::vector<int64>& start2, const uint32 startI) const
{
	//make sure the parity of each genome is the forward when the start level is different 
	if( start1[startI] < 0 ){
		//invert m_start1...
		for(uint32 j = startI; j < start1.size(); j++)
			start1[j] *= -1;
	}
	if( start2[startI] < 0 ){
		//invert m_start2...
		for( uint32 j=startI; j < start2.size(); j++)
			start2[j] *= -1;
	}
}


	static uint seq_compare_start;


}

#endif // _Match_h_

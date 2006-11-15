/*******************************************************************************
 * $Id: GenericMatch.h,v 1.10 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
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

	//GenericMatch count
	
	/**
	 * Links a subset match to this mem and this mem to the subset.
	 * @param subset The subset match which this mem overlaps
	 */
	virtual void LinkSubset( GenericMatch<BaseType>* subset );
	/** Returns a copy of set of the matches which overlap this match */
	std::set< GenericMatch<BaseType>* > GetSubsets() const { return m_subsets; }
	/** Returns a copy of set of the matches which this mem overlaps */
	std::set< GenericMatch<BaseType>* > GetSupersets() const { return m_supersets; }
	
	/** Returns a reference to the set of matches which overlap this match */
	std::set< GenericMatch<BaseType>* >& Subsets() { return m_subsets; }
	/** Returns a reference to the set of matches which this match overlaps */
	std::set< GenericMatch<BaseType>* >& Supersets() { return m_supersets; }
		
	/**
	 *  Will return true if this match evenly overlaps mhe
	 *  An even overlap implies that for each genome which this match has 
	 *  coordinates, the difference in start positions must always be the
	 *  same and this match's coordinates completely contain mhe's in those
	 *  genomes.
	 *  mhe may have coordinates in genomes which this mem does not.
	 * @param mhe The mem to check for intersection.
	 * @return True if this mem intersects with mhe.
	 */
	virtual boolean EvenlyOverlaps(const GenericMatch<BaseType>& mhe) const;
	

	virtual boolean GetUniqueStart(const GenericMatch<BaseType>& mhe, GenericMatch<BaseType>& unique_mhe) const;
	virtual boolean GetUniqueEnd(const GenericMatch<BaseType>& mhe, GenericMatch<BaseType>& unique_mhe) const;

	/** Ensures all linked supersets are really supersets and are correctly linked */
	virtual void UpdateSupersets();
	/** Unlinks this mem from any supersets and subsets it may be linked to. */
	virtual void UnlinkSelf();

	virtual void AddSubset( GenericMatch<BaseType>* overlapper );
	virtual void AddSuperset( GenericMatch<BaseType>* underlapper );
protected:
	void HomogenizeParity(std::vector<int64>& start1, std::vector<int64>& start2, const uint32 startI) const;

	virtual void UnlinkSubset( GenericMatch<BaseType>* subset );
	virtual void UnlinkSuperset( GenericMatch<BaseType>* superset );
	std::set< GenericMatch<BaseType>* > m_subsets;
	std::set< GenericMatch<BaseType>* > m_supersets;

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
void GenericMatch<BaseType>::UnlinkSubset( GenericMatch* subset ) {
	m_subsets.erase( subset );
}

template< class BaseType >
void GenericMatch<BaseType>::UnlinkSuperset( GenericMatch<BaseType>* superset ) {
	m_supersets.erase( superset );
}

template< class BaseType >
void GenericMatch<BaseType>::UnlinkSelf() 
{
	// unlink from the supersets
	typename std::set< GenericMatch<BaseType>* >::iterator match_iter = this->m_supersets.begin();
	for(; match_iter != this->m_supersets.end(); ++match_iter )
		(*match_iter)->UnlinkSubset( this );
	this->m_supersets.clear();

	// unlink from the subsets
	match_iter = this->m_subsets.begin();
	for(; match_iter != this->m_subsets.end(); ++match_iter )
		(*match_iter)->UnlinkSuperset( this );
	this->m_subsets.clear();

}

template< class BaseType >
void GenericMatch<BaseType>::UpdateSupersets() {
	typename std::set< GenericMatch<BaseType>* >::iterator super_iter = this->m_supersets.begin();

	while(super_iter != this->m_supersets.end()) {			
		GenericMatch<BaseType>* supermem = *super_iter;
		if( this->EvenlyOverlaps(*supermem) ){
			supermem->m_subsets.insert( this );
			++super_iter;
		}else {
			typename std::set< GenericMatch<BaseType>* >::iterator to_del = super_iter;
			++super_iter;
			this->m_supersets.erase( to_del );
		}
	}
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
void GenericMatch<BaseType>::LinkSubset( GenericMatch<BaseType>* subset )
{
	this->m_subsets.insert( subset );
	subset->AddSuperset( this );
}

template< class BaseType >
void GenericMatch<BaseType>::AddSuperset( GenericMatch<BaseType>* overlapper )
{
	this->m_supersets.insert(overlapper);
}

template< class BaseType >
void GenericMatch<BaseType>::AddSubset( GenericMatch<BaseType>* underlapper )
{
	this->m_subsets.insert(underlapper);
}

template< class BaseType >
boolean GenericMatch<BaseType>::GetUniqueStart(const GenericMatch<BaseType>& mhe, GenericMatch<BaseType>& unique_mhe) const
{
	uint i;
	int64 diff_i;
	int64 diff;
	uint seq_count = mhe.SeqCount();
	std::vector<int64> m_start1, m_start2;
	for( uint seqI = 0; seqI < this->SeqCount(); seqI++ ){
		m_start1.push_back( this->Start( seqI ) );
		m_start2.push_back( mhe.Start( seqI ) );
	}

	//find the first matching sequence...
	for(i=0; i < seq_count; i++){
		if(mhe.Start(i) == NO_MATCH || this->Start(i) == NO_MATCH)
			continue;
		else
			break;
	}

	uint fm = i;

	//make sure the parity of each genome is forward
//	if(this->FirstStart() != mhe.FirstStart())
		HomogenizeParity(m_start1, m_start2, i);
	diff = m_start2[i] - m_start1[i];

	//check for overlao properties
	if(diff >= this->Length() || diff <= 0)
		return false;

	//everything is ok so far, check for alignment
	int64 diff_rc = (int64)mhe.Length() - (int64)this->Length() + diff;
	for(i++; i < seq_count; i++){
		//check for consistent alignment between all genomes
		//in the case of revcomp, diff_i must equal m.length - diff
		diff_i = m_start2[i] - m_start1[i];
		//it's ok if this mem has a sequence which mhe doesn't
		if(m_start2[i] == NO_MATCH || m_start1[i] == NO_MATCH)
			continue;
		else if(m_start2[i] < 0 && diff_rc == diff_i)
			continue;
		else if(diff != diff_i )
			return false;
	}
	//it was aligned and intersects.
	unique_mhe = *this;
	if( this->Start(fm) > 0 ){
		unique_mhe.CropEnd( this->Length() - diff );
	}else{
		unique_mhe.CropEnd( mhe.Length() + diff );
	}
	return true;
}

template< class BaseType >
boolean GenericMatch<BaseType>::GetUniqueEnd(const GenericMatch<BaseType>& mhe, GenericMatch<BaseType>& unique_mhe) const
{
	uint32 i;
	int64 diff_i;
	int64 diff;
	uint32 seq_count = mhe.SeqCount();
	std::vector<int64> m_start1, m_start2;
	for( uint seqI = 0; seqI < this->SeqCount(); seqI++ ){
		m_start1.push_back( this->Start( seqI ) );
		m_start2.push_back( mhe.Start( seqI ) );
	}
	uint fm;
	
	//find the first matching sequence...
	for(i=0; i < seq_count; i++){
		if(mhe.Start(i) == NO_MATCH || this->Start(i) == NO_MATCH)
			continue;
		else
			break;
	}
	fm = i;

	//make sure the parity of each genome is the same
//	if(FirstStart() != mhe.FirstStart())
		HomogenizeParity(m_start1, m_start2, i);
	diff = m_start2[i] - m_start1[i];

	//check for overlap properties
	if(this->Length() <= mhe.Length() + diff)
		return false;


	//everything is ok so far, check for alignment
	int64 diff_rc = (int64)mhe.Length() - (int64)this->Length() + diff;
	for(i++; i < seq_count; i++){
		//check for consistent alignment between all genomes
		//in the case of revcomp, diff_i must equal m.length - diff
		diff_i = m_start2[i] - this->Start(i);
		//it's ok if this mem has a sequence which mhe doesn't
		if(m_start2[i] == NO_MATCH || m_start1[i] == NO_MATCH)
			continue;
		else if(m_start2[i] < 0 && diff_rc == diff_i)
			continue;
		else if(diff != diff_i )
			return false;
	}

	//it was aligned and intersects.
	unique_mhe = *this;
	if( this->Start(fm) > 0 ){
		unique_mhe.CropStart(mhe.Length() + diff);
	}else{
		unique_mhe.CropStart( this->Length() - diff );
	}
	return true;
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



//check if this mem evenly overlaps mhe in every sequence for which
//this match is defined.
template< class BaseType >
boolean GenericMatch<BaseType>::EvenlyOverlaps(const GenericMatch<BaseType>& mhe) const
{
	uint i;
	int64 diff_i;
	int64 diff;
	std::vector<int64> m_start1, m_start2;
	for( uint seqI = 0; seqI < this->SeqCount(); seqI++ ){
		m_start1.push_back( this->Start(seqI) );
		m_start2.push_back( mhe.Start(seqI) );
	}
	uint seq_count = m_start2.size();

	//ensure the first sequence in this mem exists in both mems...
	i = this->FirstStart();
	if(m_start2[i] == NO_MATCH)
		return false;

	//make sure the parity of each genome is the same
	if(this->FirstStart() != mhe.FirstStart())
		HomogenizeParity(m_start1, m_start2, i);
	diff = m_start2[i] - m_start1[i];

	//check for even overlap properties
	if(diff >= this->Length() || -diff >= mhe.Length())
		return false;

	//everything is ok so far, check for alignment
	int64 diff_rc = (int64)mhe.Length() - (int64)this->Length() + diff;
	for(i++; i < seq_count; i++){
		//check for consistent alignment between all genomes
		//in the case of revcomp, diff_i must equal m.length - diff
		diff_i = m_start2[i] - m_start1[i];
		//it's ok if mhe has a sequence which this mem doesn't
		//but not vice versa
		if(m_start1[i] == NO_MATCH)
			continue;
		if(m_start2[i] == NO_MATCH)
			return false;
		else if(m_start2[i] < 0 && diff_rc == diff_i)
			continue;
		else if(diff != diff_i )
			return false;
	}
	//it was an even overlap
	return true;
}


	static uint seq_compare_start;


}

#endif // _Match_h_

/*******************************************************************************
 * $Id: GappedAlignment.h,v 1.12 2004/04/19 23:10:50 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __GappedAlignment_h__
#define __GappedAlignment_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnFilter.h"
#include "libGenome/gnSequence.h"
#include "libMems/SparseAbstractMatch.h"
#include "libMems/AbstractGappedAlignment.h"

namespace mems {

class GappedAlignment : public AbstractGappedAlignment< SparseAbstractMatch<> >
{
public:
	GappedAlignment();
	GappedAlignment( uint seq_count, gnSeqI align_length );
	
	GappedAlignment* Clone() const { return new GappedAlignment( *this ); }
	GappedAlignment* Copy() const;
	virtual void Free();
	
	void SetAlignment( const std::vector< std::string >& seq_align );

	/**
	 * Writes this GappedAlignment to the specified output stream (e.g. cout).
	 */
	friend std::ostream& operator<<(std::ostream& os, const GappedAlignment& ga); //write to source.

	/**
	 * Reads a GappedAlignment from the specified input stream (e.g. cin).
	 */
	friend std::istream& operator>>(std::istream& is, GappedAlignment& ga); //read from source

	// Inherited methods from AbstractMatch:
	virtual void Invert();
	virtual void CropStart(gnSeqI crop_amount);
	virtual void CropEnd(gnSeqI crop_amount);

	virtual void CropLeft(gnSeqI crop_amount, uint seqI);
	virtual void CropRight(gnSeqI crop_amount, uint seqI);

	void GetAlignment( std::vector< bitset_t >& align_matrix ) const;

	friend const std::vector<std::string>& GetAlignment( const GappedAlignment& ga, const std::vector< genome::gnSequence* >& seq_table );

	void GetColumn( gnSeqI col, std::vector<gnSeqI>& pos, std::vector<bool>& column ) const;

protected:
	std::vector< std::string > align_matrix;
};


inline
GappedAlignment* GappedAlignment::Copy() const
{
	return m_allocateAndCopy( *this );
}
inline
void GappedAlignment::Free()
{
	m_free(this);
}

inline
void GappedAlignment::Invert(){
	const genome::gnFilter* rc_filter = genome::gnFilter::DNAComplementFilter();
	for(uint startI = 0; startI < SeqCount(); startI++)
		rc_filter->ReverseFilter( align_matrix[ startI ] );
	AbstractGappedAlignment< SparseAbstractMatch<> >::Invert();
}

inline
void GappedAlignment::CropStart(gnSeqI crop_amount){
	if( crop_amount > AlignmentLength() )
		Throw_gnEx( genome::SeqIndexOutOfBounds() );
	for( uint i=0; i < SeqCount(); i++ ){
		gnSeqI char_count = 0;
		for( gnSeqI cropI = 0; cropI < crop_amount; cropI++ )
			if( align_matrix[i][cropI] != '-' )
				char_count++;
		align_matrix[ i ] = align_matrix[ i ].substr( crop_amount );
		if( Start(i) > 0 )
			SetStart(i, Start(i) + char_count);
		SetLength(Length(i)-char_count, i);
		if( Length(i) == 0 )
			SetLeftEnd(i, NO_MATCH);
	}

	SetAlignmentLength( AlignmentLength() - crop_amount );
}

inline
void GappedAlignment::CropEnd(gnSeqI crop_amount){
	if( crop_amount > AlignmentLength() )
		Throw_gnEx( genome::SeqIndexOutOfBounds() );
	SetAlignmentLength( AlignmentLength() - crop_amount );

	for( uint i=0; i < SeqCount(); i++ ){
		gnSeqI char_count = 0;
		for( gnSeqI cropI = align_matrix[i].length() - crop_amount; cropI < align_matrix[i].length(); cropI++ )
			if( align_matrix[i][cropI] != '-' )
				char_count++;
		align_matrix[ i ] = align_matrix[ i ].substr( 0, AlignmentLength() );
		if( Start(i) < 0 )
			SetStart(i, Start(i)-char_count);
		SetLength(Length(i)-char_count, i);
		if( Length(i) == 0 )
			SetLeftEnd(i, NO_MATCH);
	}
}

inline
void GappedAlignment::CropLeft(gnSeqI crop_amount, uint seqI)
{
	// count "crop_amount" characters into seqI and crop there
	size_t left_col = 0;
	if( Orientation(seqI) == AbstractMatch::forward )
	{
		for( ; crop_amount > 0 && left_col < align_matrix[seqI].size(); ++left_col )
			if( align_matrix[seqI][left_col] != '-' )
				--crop_amount;

		CropStart(left_col);
	}else{
		left_col = align_matrix[seqI].size();
		for( ; crop_amount > 0 && left_col > 0; --left_col )
			if( align_matrix[seqI][left_col-1] != '-' )
				--crop_amount;
		CropEnd(AlignmentLength()-left_col);
	}
}

inline
void GappedAlignment::CropRight(gnSeqI crop_amount, uint seqI)
{
	// TODO: remove the dependency on Invert() since it will be slow
	Invert();
	CropLeft(crop_amount, seqI);
	Invert();
}

inline
void GappedAlignment::GetAlignment( std::vector< bitset_t >& align_matrix ) const
{
	align_matrix = std::vector< bitset_t >( this->align_matrix.size(), bitset_t(this->AlignmentLength(), false) );
	for( size_t seqI = 0; seqI < this->align_matrix.size(); seqI++ )
	{
		if( LeftEnd(seqI) == NO_MATCH )
			continue;
		for( std::string::size_type charI = 0; charI < this->align_matrix[seqI].size(); charI++ )
			if( this->align_matrix[seqI][charI] != '-' )
				align_matrix[seqI].set(charI);
	}
}

const std::vector<std::string>& GetAlignment( const GappedAlignment& ga, const std::vector< genome::gnSequence* >& seq_table );
inline
const std::vector<std::string>& GetAlignment( const GappedAlignment& ga, const std::vector< genome::gnSequence* >& seq_table )
{
	return ga.align_matrix;
}

inline
void GappedAlignment::GetColumn( gnSeqI col, std::vector<gnSeqI>& pos, std::vector<bool>& column ) const
{
	pos = std::vector<gnSeqI>(SeqCount(), NO_MATCH);
	column = std::vector<bool>(SeqCount(), false);
	for( uint seqI = 0; seqI < SeqCount(); seqI++ )
	{
		if( align_matrix[seqI][col] != '-' )
			column[seqI] = true;

		gnSeqI count = 0;
		for( size_t colI = 0; colI <= col; colI++ )
			if( align_matrix[seqI][colI] != '-' )
				count++;

		if( count > 0 )
		{
			if( Orientation(seqI) == forward )
				pos[seqI] = LeftEnd(seqI) + count - 1;
			else if( Orientation(seqI) == reverse )
				pos[seqI] = RightEnd(seqI) - count + 1;
		}
	}
}

}

#endif // __GappedAlignment_h__


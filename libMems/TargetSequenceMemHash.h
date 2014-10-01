/*******************************************************************************
 * $Id: TargetSequenceMemHash.h,v 1.3 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _TargetSequenceMemHash_h_
#define _TargetSequenceMemHash_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MemHash.h"

namespace mems {

/**
 * Finds matches that meet a particular sequence mask, e.g. 0b11111 for 5-way matches
 * Doesn't filter anything unless a mask is set using SetMask().  The
 * filter can be cleared by calling SetMask(0)
 */
class TargetSequenceMemHash : public MemHash{
public:
	TargetSequenceMemHash();
	~TargetSequenceMemHash(){};
	TargetSequenceMemHash(const TargetSequenceMemHash& mh);
	TargetSequenceMemHash& operator=( const TargetSequenceMemHash& mh );
	virtual TargetSequenceMemHash* Clone() const;
	virtual void SetTarget( uint64 target ){ this->target = target; }
protected:
	/**
	 * Subsets must include the target sequence
	 */
	virtual void FindSubsets(const Match& mhe, std::vector<Match>& subset_matches);
	virtual boolean HashMatch(std::list<idmer>& match_list);
	uint64 target;
};

}

#endif //_TargetSequenceMemHash_h_

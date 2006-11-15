/*******************************************************************************
 * $Id: RepeatHashThread.h,v 1.3 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifndef _RepeatHashThread_h_
#define _RepeatHashThread_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MemHashThread.h"
#include "libMems/RepeatHash.h"

namespace mems {

class RepeatHashThread : public MemHashThread, public RepeatHash{
public:
	RepeatHashThread();
	~RepeatHashThread();
	RepeatHashThread(const RepeatHashThread& mh);
	virtual RepeatHashThread* Clone() const;
protected:
};

}

#endif //_RepeatHashThread_h_

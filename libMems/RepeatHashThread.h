/*******************************************************************************
 * $Id: RepeatHashThread.h,v 1.3 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
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

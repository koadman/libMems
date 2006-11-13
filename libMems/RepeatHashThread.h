#ifndef _RepeatHashThread_h_
#define _RepeatHashThread_h_

#include "MemHashThread.h"
#include "RepeatHash.h"

class RepeatHashThread : public MemHashThread, public RepeatHash{
public:
	RepeatHashThread();
	~RepeatHashThread();
	RepeatHashThread(const RepeatHashThread& mh);
	virtual RepeatHashThread* Clone() const;
protected:
};

#endif //_RepeatHashThread_h_

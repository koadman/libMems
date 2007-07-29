/*******************************************************************************
 * $Id: MemHashThread.h,v 1.6 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _MemHashThread_h_
#define _MemHashThread_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
    #pragma hdrstop
#endif

// for all others, include the necessary headers (this file is usually all you
// need because it includes almost all "standard" wxWindows headers)
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include "libMems/MemHash.h"

namespace mems {

class MemHashThread : public MemHash, public wxThread{
public:
	MemHashThread();
	~MemHashThread();
	MemHashThread(const MemHash& mf);
	MemHashThread(const MemHashThread& mh);
	virtual MemHashThread* Clone() const;
	virtual void Clear();
	virtual ExitCode Entry();
	virtual void SetMemTableMutex(wxMutex *mem_mutex);
	virtual void SetMemTable(std::set<Match*, MheCompare>* shared_table);
	virtual void SetSharedMemHash(MemHash *memhash);
	virtual void SetSearchSpace(std::vector<gnSeqI>& start_points, std::vector<gnSeqI>& search_length);
	virtual void OnExit();
protected:
	virtual Match* AddHashEntry(Match& mhe);
	wxMutex *s_mutexProtectingMemTable;
	MemHash* m_memhash;

	//keep track of how much of each genome is being searched.
	std::vector<gnSeqI> start_loc;
	std::vector<gnSeqI> search_len;
};

inline
void MemHashThread::SetMemTableMutex(wxMutex *mem_mutex){
	s_mutexProtectingMemTable = mem_mutex;
}

inline
void MemHashThread::SetMemTable(set<Match*, MheCompare>* shared_table){
	mem_table = shared_table;
}

inline
void MemHashThread::SetSharedMemHash(MemHash *memhash){
	m_memhash = memhash;
}

inline
void MemHashThread::SetSearchSpace(std::vector<gnSeqI>& start_points, std::vector<gnSeqI>& search_length){
	start_loc = start_points;
	search_len = search_length;
}

}

#endif //_MemHashThread_h_

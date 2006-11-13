#ifndef _MemHashThread_h_
#define _MemHashThread_h_

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

#include "MemHash.h"

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
	virtual void SetMemTable(set<MemHashEntry*, MheCompare>* shared_table);
	virtual void SetSharedMemHash(MemHash *memhash);
	virtual void SetSearchSpace(vector<gnSeqI>& start_points, vector<gnSeqI>& search_length);
	virtual void OnExit();
protected:
	virtual MemHashEntry* AddHashEntry(MemHashEntry& mhe);
	wxMutex *s_mutexProtectingMemTable;
	MemHash* m_memhash;

	//keep track of how much of each genome is being searched.
	vector<gnSeqI> start_loc;
	vector<gnSeqI> search_len;
};

inline
void MemHashThread::SetMemTableMutex(wxMutex *mem_mutex){
	s_mutexProtectingMemTable = mem_mutex;
}

inline
void MemHashThread::SetMemTable(set<MemHashEntry*, MheCompare>* shared_table){
	mem_table = shared_table;
}

inline
void MemHashThread::SetSharedMemHash(MemHash *memhash){
	m_memhash = memhash;
}

inline
void MemHashThread::SetSearchSpace(vector<gnSeqI>& start_points, vector<gnSeqI>& search_length){
	start_loc = start_points;
	search_len = search_length;
}

#endif //_MemHashThread_h_

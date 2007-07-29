/*******************************************************************************
 * $Id: ParallelMemHash.cpp,v 1.9 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

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

#include "libMems/ParallelMemHash.h"
#include "libMems/MemHashThread.h"

using namespace std;
namespace mems {

ParallelMemHash::ParallelMemHash(){
	m_preprocess_percent = DEFAULT_PREPROCESS_PERCENT;
}

ParallelMemHash::~ParallelMemHash(){

}

ParallelMemHash::ParallelMemHash(const ParallelMemHash& mh){
	m_preprocess_percent = mh.m_preprocess_percent;
}

ParallelMemHash* ParallelMemHash::Clone() const{
	return new ParallelMemHash(*this);
}

void ParallelMemHash::Clear(){

}

void ParallelMemHash::CreateMems(uint32 thread_count){
	vector<MemHashThread*> thread_table;
	//make a self pointer so the threads will know where to
	//store their data.
	//instantiate the mutex
	wxMutex* s_mutexProtectingMemTable = new wxMutex();
	
	//preprocess the first part of the sorted mer lists to avoid
	//heavy lock contention.
	vector<gnSeqI> start_points, search_len;
	for(uint32 seqI = 0; seqI < sar_table.size(); seqI++)
		start_points.push_back(0);
	gnSeqI breakpoint = (gnSeqI) (m_preprocess_percent * GetSar(0)->Length());
	GetBreakpoint(0, breakpoint, search_len);
	//Do the initial search
	MatchFinder::SearchRange(start_points, search_len);

	gnSeqI break_size = (GetSar(0)->Length() - breakpoint) / thread_count;
	//create threads for the remaining portions of the sorted mer lists
	for(uint32 threadI = 0; threadI < thread_count; threadI++){
		MemHashThread* mh_thread = new MemHashThread(*this);
		//delete the allocated mem_table and replace it with
		//the shared one.
		mh_thread->SetSharedMemHash(this);
		mh_thread->SetMemTable(mem_table);
		mh_thread->SetMemTableMutex(s_mutexProtectingMemTable);

		//set the search space for each genome
		for(uint32 seqI = 0; seqI < sar_table.size(); seqI++)
			start_points[seqI] += search_len[seqI];
		if(threadI < thread_count - 1){
			GetBreakpoint(0, start_points[0] + break_size, search_len);
			for(uint32 seqI = 0; seqI < sar_table.size(); seqI++)
				search_len[seqI] = search_len[seqI] - start_points[seqI];
		}else{
			for(uint32 seqI = 0; seqI < sar_table.size(); seqI++)
				search_len[seqI] = GNSEQI_END;
		}
		mh_thread->SetSearchSpace(start_points, search_len);

		//create the actual thread data
		if(mh_thread->Create() != wxTHREAD_NO_ERROR){
			//bad error.  clean up and bail out.
			ErrorMsg("ParallelMemHash::CreateMems: Error creating thread.\n");
		}
		mh_thread->SetPriority(WXTHREAD_MAX_PRIORITY);
		thread_table.push_back(mh_thread);
	}

	//Run each thread
	for(uint32 threadI = 0; threadI < thread_count; threadI++)
		thread_table[threadI]->Run();

	//wait for all the threads
	//only do this if they are joinable
	for(uint32 threadI = 0; threadI < thread_count; threadI++)
		thread_table[threadI]->Wait();
}

} // namespace mems

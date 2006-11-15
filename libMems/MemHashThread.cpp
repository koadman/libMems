/*******************************************************************************
 * $Id: MemHashThread.cpp,v 1.9 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MemHashThread.h"

using namespace std;
namespace mems {

MemHashThread::MemHashThread() : wxThread(wxTHREAD_JOINABLE){
	s_mutexProtectingMemTable = NULL;
	m_memhash = NULL;
	if(mem_table != NULL)
		delete[] mem_table;
}

MemHashThread::~MemHashThread(){

}

MemHashThread::MemHashThread(const MemHash& mh) : MemHash(mh), wxThread(wxTHREAD_JOINABLE){
	s_mutexProtectingMemTable = NULL;
	m_memhash = NULL;
	if(mem_table != NULL)
		delete[] mem_table;
}

MemHashThread::MemHashThread(const MemHashThread& mh) : MemHash(mh), wxThread(wxTHREAD_JOINABLE){
	s_mutexProtectingMemTable = mh.s_mutexProtectingMemTable;
	m_memhash = mh.m_memhash;
	start_loc = mh.start_loc;
	search_len = mh.search_len;
}

MemHashThread* MemHashThread::Clone() const{
	return new MemHashThread(*this);
}

void MemHashThread::Clear(){
	start_loc.clear();
	search_len.clear();
	m_memhash = NULL;
	s_mutexProtectingMemTable = NULL;
}

wxThread::ExitCode MemHashThread::Entry(){
	while( !MatchFinder::SearchRange(start_loc, search_len) )
		;
	return 0;
}

// This is a parallel implementation of AddHashEntry
// First it checks the hash table for existence of the new seed,
// if it is not found then it locks the mutex and tries to insert
// the seed again.  If successful, it extends the seed into a full
// match.
Match* MemHashThread::AddHashEntry(Match& mhe){
	//first compute which hash table bucket this is going into
	int64 offset = mhe.Offset();

	uint32 bucketI = ((offset % table_size) + table_size) % table_size;
	
	set<Match*, MheCompare>::iterator insert_he;
	insert_he = mem_table[bucketI].find(&mhe);
	if( insert_he != mem_table[bucketI].end()){
		// it already exists!
		m_collision_count++;
		return false;
	}

	//if we made it this far there were no collisions
	//lock the mutex and try to insert it again.

	wxMutexLocker wxml(*s_mutexProtectingMemTable);
	// let the MemHash class do the rest.
	return MemHash::AddHashEntry( mhe );

}

void MemHashThread::OnExit(){
	wxMutexLocker wxml(*s_mutexProtectingMemTable);
//	m_memhash->m_mem_count += m_mem_count;
//	for(uint32 bucketI=0; bucketI < table_size; bucketI++)
//		m_memhash->mem_table_count[bucketI] += mem_table_count[bucketI];
}

} // namespace mems

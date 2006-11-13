#ifndef _FileSML_h_
#define _FileSML_h_

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

#include "gnSequence.h"
#include "SortedMerList.h"
#include <fstream>
#include <vector>
#include <string>

//sequence database size will be
//base_count / 4 + base_count * 12 bytes

#define DEFAULT_MEMORY_MINIMUM 20971520  //~20 Megabytes

class FileSML : public SortedMerList
{
public:
	FileSML() : SortedMerList() {file_mutex = new wxMutex();};
	FileSML& operator=(const FileSML& sa);
	virtual FileSML* Clone() const = 0;
	/**
	 * Loads an existing sorted mer list from a file on disk.
	 * @param fname The name of the file to load
	 * @throws FileNotOpened thrown if the file could not be opened
	 * @throws FileUnreadable thrown if the file was corrupt or not a sorted mer list
	 */
	virtual void LoadFile(const string& fname);
	/**
	 * Creates large sorted mer lists which do not fit entirely in memory.
	 * BigCreate uses an external mergesort to create large sorted mer lists.
	 * It will divide the data a number of times specified by the split_levels
	 * parameter.  Each split is written to temp files on disk and merged.
	 * @param seq The sequence to create an SML for.
	 * @param split_levels The number of times to divide the sequence in half.
	 * @param mersize The size of the mers to sort on.
	 * @see FileSML::Create
	 */
	virtual void BigCreate(const gnSequence& seq, const uint32 split_levels, const uint32 mersize = DNA_MER_SIZE);
	virtual void Create(const gnSequence& seq, const uint32 mersize = DNA_MER_SIZE);
	virtual boolean Read(vector<bmer>& readVector, gnSeqI size, gnSeqI offset = 0);
	virtual void Merge(SortedMerList& sa, SortedMerList& sa2);

	virtual bmer operator[](const gnSeqI index);

	virtual uint32 UniqueMerCount();
	virtual void SetDescription(const string& d);
	virtual void SetID(const sarID_t d);
	
	virtual uint32 FormatVersion();
	static uint64 MemoryMinimum();
	virtual void RadixSort(vector<bmer>& s_array);
protected:
	/**
	 * Reopens the sarfile fstream in read/write mode
	 * @throws FileNotOpened thrown if the file could not be opened for writing
	 */
	virtual void OpenForWriting();
	/**
	 * Writes the SML header to disk
	 * @throws FileNotOpened thrown if the file could not be opened for writing
	 * @throws IOStreamFailed thrown if an error occurred writing the data
	 */
	virtual boolean WriteHeader();
	/**
	 * Calculates and returns the amount of memory needed to create a sorted
	 * mer list for a sequence of the specified length.
	 * @param len The length of the sequence
	 * @return The amount of memory needed in bytes.
	 */
	virtual uint64 GetNeededMemory(gnSeqI len) = 0;

	string filename;
	fstream sarfile;
	uint64 sarray_start_offset;
	wxMutex* file_mutex;
};

inline
uint32 FileSML::FormatVersion(){
	static uint32 f_version = 2;
	return f_version;
}

inline
uint64 FileSML::MemoryMinimum(){
	static uint32 m_minimum = DEFAULT_MEMORY_MINIMUM;
	return m_minimum;
}

#endif   //_FileSML_h_

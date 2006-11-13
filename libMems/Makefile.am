
if DEBUG
MEMLIB = libMemsd
D_CXXFLAGS = -Wall -g -O2 -D__GNDEBUG__
else
MEMLIB = libMems
P_CXXFLAGS = -O2
endif

CXXFLAGS = $(D_CXXFLAGS) $(P_CXXFLAGS) `@CONFIG_WXWINDOWSDIR@/bin/wx-config --cxxflags` 

LIBMEMS_SRC = \
MemScorer.cpp        RepeatHash.cpp        \
DNAFileSML.cpp       MatchFinder.cpp      MemSubsets.cpp \
DNAMemorySML.cpp     MemorySML.cpp        SortedMerList.cpp \
FileSML.cpp          MemHash.cpp          MemHashEntry.cpp \
MemHashThread.cpp    MemList.cpp \
ParallelMemHash.cpp \
Match.cpp            MatchList.cpp        Aligner.cpp 

if DEBUG
lib_LIBRARIES = libMemsd.a
libMemsd_a_SOURCES = $(LIBMEMS_SRC)
lib_LTLIBRARIES = libMemsd.la
libMemsd_la_SOURCES = $(LIBMEMS_SRC)
else
lib_LIBRARIES = libMems.a
libMems_a_SOURCES = $(LIBMEMS_SRC)
lib_LTLIBRARIES = libMems.la
libMems_la_SOURCES = $(LIBMEMS_SRC)
endif

#check_PROGRAMS = TestRevComp
#TestRevComp_SOURCES = TestRevComp.cpp
#TestRevComp_LDADD = $(top_builddir)/src/libGenome.a

EXTRA_PROGRAMS = Mormems MakeSuffixArray
Mormems_SOURCES = mormems.cpp
Mormems_INCLUDES = -I$(top_srcdir)/include/ -I@CONFIG_LIBGENOMEDIR@/include -I@CONFIG_LIBGENOMEDIR@/include/gn -I@CONFIG_WXWINDOWSDIR@/include/
Mormems_LDADD = $(top_builddir)/src/libMems.a -lGenome -lwx_gtk -L@CONFIG_LIBGENOMEDIR@/lib -L@CONFIG_WXWINDOWSDIR@/lib/

MakeSuffixArray_SOURCES = MakeSuffixArray.cpp
MakeSuffixArray_INCLUDES = -I$(top_srcdir)/include/ -I@CONFIG_LIBGENOMEDIR@/include -I@CONFIG_LIBGENOMEDIR@/include/gn -I@CONFIG_WXWINDOWSDIR@/include/
MakeSuffixArray_LDADD = $(top_builddir)/src/libMems.a -lGenome -lwx_gtk -L@CONFIG_LIBGENOMEDIR@/lib -L@CONFIG_WXWINDOWSDIR@/lib/

INCLUDES = -I$(top_srcdir)/include/ -I@CONFIG_LIBGENOMEDIR@/include -I@CONFIG_LIBGENOMEDIR@/include/gn  

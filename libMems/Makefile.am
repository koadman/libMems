SUBDIRS = dmSML

if DEBUG
D_CXXFLAGS = -Wall -g -D__GNDEBUG__
endif
AM_CXXFLAGS = $(D_CXXFLAGS)  

LIBMEMS_H = \
MemSubsets.h       RepeatHash.h      MatchHashEntry.h \
DNAFileSML.h     MemorySML.h         MatchProjectionAdapter.h \
DNAMemorySML.h  MatchFinder.h           SortedMerList.h IntervalList.h \
FileSML.h      gnAlignedSequences.h  Interval.h        ClustalInterface.h  \
MemHash.h      MemScorer.h          AbstractMatch.h    SlotAllocator.h \
Aligner.h   Match.h     MatchList.h Matrix.h NumericMatrix.h \
Islands.h   MaskedMemHash.h   SlotAllocator.h SeedMasks.h GappedAlignment.h \
MuscleInterface.h GappedAligner.h PhyloTree.h SparseAbstractMatch.h \
DenseAbstractMatch.h RepeatMatch.h UngappedLocalAlignment.h \
AbstractGappedAlignment.h CompactGappedAlignment.h HybridAbstractMatch.h \
twister.h SubstitutionMatrix.h RepeatHashCat.h RepeatMatchList.h

LIBMEMS_SRC = \
MemScorer.cpp        RepeatHash.cpp       ClustalInterface.cpp \
DNAFileSML.cpp       MatchFinder.cpp      MemSubsets.cpp \
DNAMemorySML.cpp     MemorySML.cpp        SortedMerList.cpp \
FileSML.cpp          MemHash.cpp          MatchHashEntry.cpp \
Interval.cpp	     IntervalList.cpp     twister.c \
gnAlignedSequences.cpp                     \
MatchList.cpp        Aligner.cpp \
Islands.cpp          MaskedMemHash.cpp    GappedAlignment.cpp \
MuscleInterface.cpp  PhyloTree.cpp         \
RepeatHashCat.cpp    RepeatMatchList.cpp  RepeatMatch.cpp 


library_includedir=$(includedir)/$(GENERIC_LIBRARY_NAME)-$(GENERIC_API_VERSION)/$(GENERIC_LIBRARY_NAME)

library_include_HEADERS = $(LIBMEMS_H) 

# build libraries with gcc (no suffix)
lib_LTLIBRARIES = libMems-1.5.la
libMems_1_5_la_SOURCES = $(LIBMEMS_SRC)

libMems_1_5_la_LDFLAGS= -version-info $(GENERIC_LIBRARY_VERSION) -release $(GENERIC_RELEASE) 
libMems_1_5_la_LIBADD= dmSML/libdmSML.la


#EXTRA_PROGRAMS = TestSML TestSMLstatic


#TestSMLstatic_SOURCES = TestSML.cpp
#TestSMLstatic_INCLUDES = -I$(top_srcdir)/include/ `wx-config --cxxflags`
#TestSMLstatic_LDFLAGS = -static $(top_builddir)/libMems/libMems.a `wx-config --static --libs`

INCLUDES = $(DEPS_CFLAGS) @BOOST_CPPFLAGS@

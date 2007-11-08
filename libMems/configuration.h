#ifndef __libMems_configuration_h__
#define __libMems_configuration_h__

#if defined(WIN32)||defined(WIN64)

// set the mems library name to include based on the configuration...

#if defined(WIN64)&&defined(NDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "mems64omp.lib")
#else if defined(WIN64)&&!defined(NDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "mems64fdomp.lib")
#else if defined(WIN32)&&defined(NDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "memsomp.lib")
#else if defined(WIN32)&&!defined(NDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "memsfdomp.lib")
#if defined(WIN64)&&defined(NDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "mems64.lib")
#else if defined(WIN64)&&!defined(NDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "mems64fd.lib")
#else if defined(WIN32)&&defined(NDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "mems.lib")
#else if defined(WIN32)&&!defined(NDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "memsfd.lib")
#endif


#endif

#endif // __libMems_configuration_h__
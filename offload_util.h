#ifndef H_OFFLOADUTIL
#define H_OFFLOADUTIL

// macro the offload declspec -- 1. easier to write, 2. ctags doesn't like declspec
#ifdef __INTEL_OFFLOAD
#define _OFFLOADABLE __declspec(target(mic))
#else
#define _OFFLOADABLE 
#endif

// TODO: this should be a runtime parameter
#define MIC_NCORES 59

// Intel LEO shorthands
#define ALLOC alloc_if(1)
#define FREE free_if(1)
#define RETAIN free_if(0)
#define REUSE alloc_if(0)

// global vars
#ifdef __cplusplus
extern C {
#endif
extern int offload_target;
#ifdef __cplusplus
}
#endif

#endif

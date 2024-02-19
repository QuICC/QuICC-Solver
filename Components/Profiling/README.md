# Profiler

## How to use
You need to set the following cmake variables
```bash
cmake <path/to/QuICC> \
-DQUICC_PROFILE=ON \                        # activate
-DQUICC_PROFILER_BACKEND=<native|likwid> \  # select profiler
-DQUICC_PROFILER_LEVEL=[0|1|...]            # choose granularity of regions profiled
```
In the code you need to initialize/finalize the profiler and mark the regions to be profiled
(see tests as example)
```bash
using namespace QuICC::Profiler;
Initialize();
...
RegionStart("myRegion");
...
RegionEnd("myRegion");
...
Finalize();
```
alternatively one can use RAII style regions
```bash
using namespace QuICC::Profiler;
Initialize();
...
{
    RegionFixture myFix("myRegion");
    ...
}
...
Finalize();
```
in both the above cases one can specify a depth level of activation (default=0)
```bash
...
RegionStart<2>("myRegion");
...
RegionEnd<2>("myRegion");
...
{
    RegionFixture<1> myFix("myRegion");
    ...
}
...
```
so that different granularities can be achieved.
The counter of collected values for all regions can be reset:
```bash
RegionResetAll();
``` 
## Existing backends

### native
A simple backend implemented direclty in QuICC that does not depend on external
libraries.
Timing is based on `std::chrono::steady_clock` and memory tracking is based on
the resident set size which is accessed through the UNIX `rusage` struct.

### likwid
Please refer to the likwid project [manual](https://github.com/RRZE-HPC/likwid).

#### MPI syncronization
To attribute uncore counters (for instance, memory traffic) to the correct region, one needs to syncronize MPI processes.
A MPI barrier can be added (before and/or after a region) at runtime by setting an enviroment variable as follows
```bash
export QUICC_PROFILER_MPI_BARRIER="after;before"
```

## How to implement new backends
- add backend to `QUICC_PROFILER_BACKEND` in `Debug.cmake`
- add implementation in `Components/Profiling/Profiler/Interface.hpp`;
    you need to implement the following methods
```bash
...
#elif defined QUICC_PROFILER_BACKEND_NEWAWESOMEBACKEND

namespace details {
    inline static void Initialize(){...}
    inline static void Finalize(){...}
    inline static void RegionStart(const std::string& name){...}
    inline static void RegionEnd(const std::string& name){...}
    inline static void RegionResetAll(){...}
}
...
```
- add test

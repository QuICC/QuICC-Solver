/**
 * @file ProfilerMacro.h
 * @brief Preprocessor macros to setup profiling timers if requested in configuration 
 */

#ifndef QUICC_DEBUG_PROFILERMACRO_H
#define QUICC_DEBUG_PROFILERMACRO_H

#ifdef QUICC_PROFILE
   #ifdef QUICC_PROFILE_PERCORE
   #include "QuICC/Io/Ascii/DirectAsciiWriter.hpp"
   #endif //QUICC_PROFILE_PERCORE
   #include "QuICC/Debug/Profiler/Tools.hpp"
   #ifdef QUICC_MPI
      // include MPI profiler
      #include "QuICC/Debug/Profiler/MpiProfiler.hpp"

      namespace QuICC {
         /// Typedef for a profiler
         typedef Debug::Profiler::MpiProfiler ProfilerMacro;
      }
   #else
      // include serial profiler
      #include "QuICC/Debug/Profiler/SerialProfiler.hpp"

      namespace QuICC {
         /// Typedef for a profiler
         typedef Debug::Profiler::SerialProfiler ProfilerMacro;
      }
   #endif // QUICC_MPI

   /// Define profiler initialisation macro function
   #define ProfilerMacro_init()  QuICC::ProfilerMacro::init()

   /// Reset profiler reset macro function
   #define ProfilerMacro_reset()  QuICC::ProfilerMacro::reset()

   // Create profiler ID variable
   #define ProfilerMacro_var(V,P)  auto V = P

   /// Define profiler start macro function
   #define ProfilerMacro_start(P)  QuICC::ProfilerMacro::start(P)
   #define ProfilerMacro_start_lvl1(P,N)  QuICC::ProfilerMacro::start(static_cast<Debug::Profiler::BreakPoint>(static_cast<int>(P)+N*Debug::Profiler::LVL1))
   #define ProfilerMacro_start_lvl2(P,N)  QuICC::ProfilerMacro::start(static_cast<Debug::Profiler::BreakPoint>(static_cast<int>(P)+N*Debug::Profiler::LVL2))
   #define ProfilerMacro_start_lvl3(P,N)  QuICC::ProfilerMacro::start(static_cast<Debug::Profiler::BreakPoint>(static_cast<int>(P)+N*Debug::Profiler::LVL3))
   #define ProfilerMacro_start_lvl2_lvl3(P,M,N)  QuICC::ProfilerMacro::start(static_cast<Debug::Profiler::BreakPoint>(static_cast<int>(P)+M*Debug::Profiler::LVL2+N*Debug::Profiler::LVL3))

   /// Define profiler stop macro function
   #define ProfilerMacro_stop(P)  QuICC::ProfilerMacro::stop(P)
   #define ProfilerMacro_stop_lvl1(P,N)  QuICC::ProfilerMacro::stop(static_cast<Debug::Profiler::BreakPoint>(static_cast<int>(P)+N*Debug::Profiler::LVL1))
   #define ProfilerMacro_stop_lvl2(P,N)  QuICC::ProfilerMacro::stop(static_cast<Debug::Profiler::BreakPoint>(static_cast<int>(P)+N*Debug::Profiler::LVL2))
   #define ProfilerMacro_stop_lvl3(P,N)  QuICC::ProfilerMacro::stop(static_cast<Debug::Profiler::BreakPoint>(static_cast<int>(P)+N*Debug::Profiler::LVL3))
   #define ProfilerMacro_stop_lvl2_lvl3(P,M,N)  QuICC::ProfilerMacro::stop(static_cast<Debug::Profiler::BreakPoint>(static_cast<int>(P)+M*Debug::Profiler::LVL2+N*Debug::Profiler::LVL3))

   #ifdef QUICC_PROFILE_PERCORE
      /// Define profiler printInfo macro function
      #define ProfilerMacro_printInfo()   Array ts;\
                                          ProfilerMacro::getTimings(ts);\
                                          std::stringstream sid;\
                                          sid << QuICCEnv().id();\
                                          Io::Ascii::DirectAsciiWriter prof_out("prof_out", "." + sid.str(), "Detailed profiling for rank " + sid.str(), "ASCII", "1.0");\
                                          prof_out.initDebug();\
                                          QuICC::Debug::Profiler::Tools::printInfo(&prof_out)\
                                          prof_out.finalizeDebug();
   #else
      /// Define profiler printInfo macro function
      #define ProfilerMacro_printInfo()  QuICC::Debug::Profiler::Tools::printInfo()
   #endif // QUICC_PROFILE_PERCORE

#else
   /// Define profiler initialisation macro function
   #define ProfilerMacro_init()

   /// Define profiler reset macro function
   #define ProfilerMacro_reset()

   /// Define empty profiler start macro function
   #define ProfilerMacro_start(P)  
   #define ProfilerMacro_start_lvl1(P,N)  
   #define ProfilerMacro_start_lvl2(P,N)  
   #define ProfilerMacro_start_lvl3(P,N)  
   #define ProfilerMacro_start_lvl2_lvl3(P,M,N)  

   /// Define empty profiler stop macro function
   #define ProfilerMacro_stop(P)  
   #define ProfilerMacro_stop_lvl1(P,N)  
   #define ProfilerMacro_stop_lvl2(P,N)  
   #define ProfilerMacro_stop_lvl3(P,N)  
   #define ProfilerMacro_stop_lvl2_lvl3(P,M,N)  

   /// Define empty profiler printInfo macro function
   #define ProfilerMacro_printInfo()  
#endif // QUICC_PROFILE

#endif // QUICC_DEBUG_PROFILERMACRO_H

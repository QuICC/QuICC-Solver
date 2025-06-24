/**
 * @file MergeTest.hpp
 * @brief Run Catch2 session with merged tests
 */

#ifndef QUICC_TESTSUITE_TRANSFORM_MERGETEST_HPP
#define QUICC_TESTSUITE_TRANSFORM_MERGETEST_HPP

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "QuICC/TestSuite/Transform/TestArgs.hpp"
#include "TestSuite/Io.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {
/// @brief namespace for TestSuite common utilities
namespace TestSuite {

namespace Transform {

   struct MergeInfo 
   {
      std::string file;
      std::size_t jid;
      std::size_t jN;

      MergeInfo();
      ~MergeInfo() = default;
   };

/// @brief Set standard argc/argv based on options vector
template <typename T> int runMergedTests(const MergeInfo& info, T& session, TestArgs& testArgs);

template <typename T>
   int runMergedTests(const MergeInfo& info, T& session, TestArgs& testArgs)
{
   int returnCode = -1;
   std::vector<std::string> commands;
   QuICC::TestSuite::readLines(commands, info.file, info.jid, info.jN);

   std::vector<std::string> failed;
   std::size_t counter = 0;
   for(const auto& command: commands)
   {
      // Reset session
      QuICC::Profiler::RegionResetAll();
      session.useConfigData({0});
      std::string testType = "";
      testArgs.clear();

      // Process next line
      std::vector<std::string> options;
      QuICC::TestSuite::splitLine(options, command, ';');

      std::cout << "[" << counter << "/" << commands.size() << "]:";
      for(auto&& o: options)
      {
         std::cout << ' ' << o;
      }
      std::cout << std::endl;

      int run_argc;
      std::vector<char *> run_argv;
      QuICC::TestSuite::getCommand(run_argc, run_argv, options);

      // Let Catch (using Clara) parse the command line
      int ret = session.applyCommandLine( run_argc, run_argv.data() );
      if( ret != 0 ) // Indicates a command line error
      {
         return ret;
      }

      // Set test from command line
      if(testType != "")
      {
         testArgs.setType(testType);
      }

      if(testArgs.params.size() > 0)
      {
         testArgs.useDefault = false;
      }

      ret = session.run();
      counter++;

      if(ret != 0)
      {
         failed.push_back(options.at(1));
      }
   }

   returnCode = failed.size();

   if(returnCode == 0)
   {
      std::cout << "All merged tests passed" << std::endl;
   }
   else
   {
      std::cout << "Following merged tests failed: " << std::endl;
      for(auto&& tag: failed)
      {
         std::cout << "   " << tag << std::endl;
      }
   }

   return returnCode;
}

} // namespace Transform
} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_TRANSFORM_MERGETEST_HPP

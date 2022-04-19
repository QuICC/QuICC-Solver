#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <sstream>
#include <functional>

// save
#if defined QUICC_PROFILE_LEVEL
    #define GLOBAL_QUICC_PROFILE_LEVEL QUICC_PROFILE_LEVEL
    #undef QUICC_PROFILE_LEVEL
#endif
#define QUICC_PROFILE_LEVEL 1

#define QUICC_PROFILER_BACKEND_TESTER
#include "Profiler/Interface.hpp"


// helper class to capture std::cout
struct CaptureStdOut{
    CaptureStdOut()
    {
        // Redirect std::cout to buffer
        prevcoutbuf = std::cout.rdbuf(buffer.rdbuf());
    }
    ~CaptureStdOut()
    {
        // Restore original buffer before exiting
        std::cout.rdbuf(prevcoutbuf);
    }

    int compare(const std::string& expected)
    {
        // Use the string value of buffer to compare against expected output
        std::string stdout = buffer.str();
        return stdout.compare(expected);
    }

    std::stringstream buffer;
    std::streambuf* prevcoutbuf;

};

// helper class to check the captured output
template <class T>
void captureAndCheck(T F, const std::string& ref)
{
    CaptureStdOut checker;
    F();
    INFO("Expected: \""+ref+'\"');
    INFO("Found: \""+checker.buffer.str()+'\"');
    REQUIRE(checker.compare(ref) == 0);
}


TEST_CASE("Open Close", "[OpenClose]")
{
    using namespace QuICC::Profiler;
    captureAndCheck(std::bind(Initialize), "initall\n");
    captureAndCheck(std::bind(Finalize), "endall\n");
}

TEST_CASE("Start Stop region with String", "[StartStopString]")
{
    using namespace QuICC::Profiler;
    std::string region = "Important";
    captureAndCheck(std::bind(RegionStart, region), "start: "+region+'\n');
    captureAndCheck(std::bind(RegionEnd, region), "end: "+region+'\n');
}

TEST_CASE("Start Stop region with Literal", "[StartStopLiteral]")
{
    using namespace QuICC::Profiler;
    captureAndCheck(std::bind(RegionStart, "Important"), std::string("start: ")+"Important"+'\n');
    captureAndCheck(std::bind(RegionEnd, "Important"), std::string("end: ")+"Important"+'\n');
}

TEST_CASE("RAII", "[RAII]")
{
    std::string name = "RAII";
    captureAndCheck([name](){QuICC::Profiler::RegionFixture aFix(name);},
        "start: "+name+'\n'+
        "end: "+name+'\n');
}

TEST_CASE("Nested Levels: on, on, off", "[Levels]")
{
    captureAndCheck([]()
    {
        using namespace QuICC::Profiler;
        RegionFixture<0> mainFix("Main");
        {
            RegionFixture<1> nestedFix("Nested");
            {
                RegionFixture<2> innerFix("Innermost");
            }
        }
    },
    "start: Main\nstart: Nested\nend: Nested\nend: Main\n");
}

// cleanup & restore
#undef QUICC_PROFILER_BACKEND_TESTER
#undef QUICC_PROFILE_LEVEL
#if defined GLOBAL_QUICC_PROFILE_LEVEL
    #define QUICC_PROFILE_LEVEL GLOBAL_QUICC_PROFILE_LEVEL
    #undef GLOBAL_QUICC_PROFILE_LEVEL
#endif

/**
 * @file Tracker.hpp
 * @brief Track time and memory usage for a region.
 * Print summary with min/max/avg across mpi ranks
 *
 */

#include "Tracker.hpp"

#include <array>
#include <iomanip>
#include <utility>

#ifdef QUICC_PROFILE_NATIVE_WRITER_HIGHFIVE
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#endif

#include "PeakRss.hpp"
#include "gitHash.hpp"

namespace QuICC {
namespace Profiler {


#ifdef QUICC_MPI
void Tracker::init(MPI_Comm comm)
{
    mComm = comm;
}
#else
void Tracker::init(){}
#endif

void Tracker::initRegion(const std::string& regionName)
{
    mRegions.insert({regionName, Tracker::tracking_t{
        std::array<double, collectSize>(), 0, 0, 0, SteadyTimer(), MemoryTracker()}});
}

void Tracker::start(const std::string& regionName)
{
    // search for region and register if doesn't exist
    auto reg = mRegions.find(regionName);
    if (reg == mRegions.end())
    {
        initRegion(regionName);
        // retry
        start(regionName);
    }
    else {
        // start memory and time tracker (time second!)
        std::get<tracking::memTracker>(reg->second).start();
        std::get<tracking::timTracker>(reg->second).start();
    }
}

void Tracker::stop(const std::string& regionName)
{
    // stop time and mem tracker (time first!)
    auto reg = mRegions.find(regionName);
    std::get<tracking::timTracker>(reg->second).stop();
    std::get<tracking::memTracker>(reg->second).stop();

    auto count = std::get<tracking::count>(reg->second);
    // time
    auto id = count % collectSize;
    std::get<tracking::time>(reg->second)[id] =
        std::get<tracking::timTracker>(reg->second).time();
    // update count
    count += 1;
    std::get<tracking::count>(reg->second) = count;
    // high watermark memory (MB)
    std::get<tracking::memory>(reg->second) =
        std::max(std::get<tracking::memory>(reg->second), PeakRss() * 1.0e-6);
    // memory delta (MB)
    std::get<tracking::memoryDelta>(reg->second) +=
        std::get<tracking::memTracker>(reg->second).diff() * 1.0e-6;

}

Tracker::tracking_t Tracker::get(const std::string& regionName)
{
    auto reg = mRegions.find(regionName);
    return reg->second;
}

void Tracker::print_hdf5()
{
#ifdef QUICC_PROFILE_NATIVE_WRITER_HIGHFIVE

    std::string outFile = "profile.hdf5";

    int rank{};
    int nRanks{1};

    #if QUICC_MPI
    MPI_Comm_rank(mComm, &rank);
    MPI_Comm_size(mComm, &nRanks);
    #endif

    // Write in a format compatible to conduit
    using namespace HighFive;

    File file(outFile, File::ReadWrite | File::Create | File::Truncate
    #if QUICC_MPI
    , MPIOFileDriver(mComm, MPI_INFO_NULL)
    #endif
    );

    // Write metadata
    auto gi = file.createGroup("info");
    std::vector<size_t> dimScalar(1);
    dimScalar[0] = std::size_t(1);
    DataSet dataRanks = gi.createDataSet<int>("ranks", DataSpace(dimScalar));
    dataRanks.write(nRanks);
    DataSet dataCommit = gi.createDataSet<char[10]>("git-commit", DataSpace(dimScalar));
    dataCommit.write(Framework::gitHash);

    // Write timings
    auto gT = file.createGroup("timings");
    for (auto reg = mRegions.begin(); reg != mRegions.end(); ++reg)
    {
        // compute actual sample size
        auto count = std::get<tracking::count>(reg->second);
        auto sampleSize = (count > collectSize) ? collectSize : count;

        // scalars
        auto gR = gT.createGroup(reg->first);
        gR.createDataSet<std::uint32_t>("count", DataSpace::From(count)).write(count);
        gR.createDataSet<std::uint32_t>("sampleSize", DataSpace::From(sampleSize))
            .write(sampleSize);
        gR.createDataSet<double>("memory", DataSpace::From(dimScalar))
            .write(&std::get<tracking::memory>(reg->second));
        gR.createDataSet<double>("memoryDelta", DataSpace::From(dimScalar))
            .write(&std::get<tracking::memoryDelta>(reg->second));

        // timings
        std::vector<size_t> dimSample = {nRanks*sampleSize};
        gR.createDataSet<double>("time", DataSpace(dimSample))
            .select({std::size_t(rank)*sampleSize}, {sampleSize})
            .write(std::get<tracking::time>(reg->second).data());

        // wait for everyone to be done
        MPI_Barrier(mComm);
    }
#endif
}

#ifdef QUICC_MPI
void Tracker::print(std::ostream& os)
{

    int rank, nRanks;
    MPI_Comm_rank(mComm, &rank);
    MPI_Comm_size(mComm, &nRanks);

    for (auto reg = mRegions.begin(); reg != mRegions.end(); ++reg)
    {
        // collect stats memory
        constexpr std::size_t stat_size = 3;
        std::array<double, stat_size> data{}, data_min{}, data_max{}, data_avg{};

        // avg/min/max time
        auto count = std::get<tracking::count>(reg->second);
        auto sampleSize = (count > collectSize) ? collectSize : count;
        double timeAvg{};
        double timeMin{1e+300};
        double timeMax{-1e+300};
        for (std::size_t i = 0; i < sampleSize; ++i)
        {
            auto t = std::get<tracking::time>(reg->second)[i];
            timeMin = std::min(t, timeMin);
            timeMax = std::max(t, timeMax);
            timeAvg += t;
        }
        timeAvg /= sampleSize;

        // abusing the fact that data types are mapped to the lower part of the enum
        data[tracking::memory] = std::get<tracking::memory>(reg->second);
        data[tracking::memoryDelta] = std::get<tracking::memoryDelta>(reg->second);
        data[tracking::time] = timeMin;
        MPI_Reduce(data.data(), data_min.data(), stat_size, MPI_DOUBLE, MPI_MIN, 0, mComm);
        data[tracking::time] = timeMax;
        MPI_Reduce(data.data(), data_max.data(), stat_size, MPI_DOUBLE, MPI_MAX, 0, mComm);
        data[tracking::time] = timeAvg;
        MPI_Reduce(data.data(), data_avg.data(), stat_size, MPI_DOUBLE, MPI_SUM, 0, mComm);
        for (auto& d : data_avg)
        {
            d /= nRanks;
        }

        if (rank == 0)
        {
            os << reg->first
                << '\t' << std::get<tracking::count>(reg->second) << " times\n"
                << "\t stats by rank                         min             max             avg\n"
                << "\t high watermark memory (MB):"
                << '\t' << std::setw(10) << data_min[tracking::memory]
                << '\t' << std::setw(10) << data_max[tracking::memory]
                << '\t' << std::setw(10) << data_avg[tracking::memory]
                << '\n'
                << "\t memory delta (MB)         :"
                << '\t' << std::setw(10) << data_min[tracking::memoryDelta]
                << '\t' << std::setw(10) << data_max[tracking::memoryDelta]
                << '\t' << std::setw(10) << data_avg[tracking::memoryDelta]
                << '\n'
                << "\t time (s)                  :"
                << '\t' << std::setw(10) << data_min[tracking::time]
                << '\t' << std::setw(10) << data_max[tracking::time]
                << '\t' << std::setw(10) << data_avg[tracking::time]
                << '\n';
        }
    }
}
#else
void Tracker::print(std::ostream& os)
{
    for (auto reg = mRegions.begin(); reg != mRegions.end(); ++reg)
    {
        // average time
        auto count = std::get<tracking::count>(reg->second);
        auto sampleSize = (count > collectSize) ? collectSize : count;
        double timeAvg{};
        for (std::size_t i = 0; i < sampleSize; ++i)
        {
            timeAvg += std::get<tracking::time>(reg->second)[i];
        }
        timeAvg /= sampleSize;

        os << reg->first
            << '\t' << std::get<tracking::count>(reg->second) << " times\n"
            << "\t high watermark memory (MB):"
            << '\t' << std::setw(10) << std::get<tracking::memory>(reg->second)
            << '\n'
            << "\t memory delta (MB)         :"
            << '\t' << std::setw(10) << std::get<tracking::memoryDelta>(reg->second)
            << '\n'
            << "\t time (s)                  :"
            << '\t' << std::setw(10) << timeAvg
            << '\n';
    }
}
#endif

// static members init
std::map<std::string, Tracker::tracking_t>
    Tracker::mRegions{};

#ifdef QUICC_MPI
MPI_Comm Tracker::mComm = MPI_COMM_NULL;
#endif

} // namespace Profiler
} // namespace QUICC

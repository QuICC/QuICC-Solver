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
#include "Framework/gitHash.hpp"
#if defined QUICC_HAS_CUDA_BACKEND
#include <cuda_runtime_api.h>
#endif

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
    else 
    {
        #if defined QUICC_HAS_CUDA_BACKEND
        // sync cuda kernel
        // WARNING
        // this does not allow for kernel overlap
        // to be substituted with cudaEventCreate/Record/Sync/Elapsed
        cudaDeviceSynchronize();
        #endif
        // start memory and time tracker (time second!)
        std::get<tracking::memTracker>(reg->second).start();
        std::get<tracking::timTracker>(reg->second).start();
    }
}

void Tracker::stop(const std::string& regionName)
{
    // stop time and mem tracker (time first!)
    auto reg = mRegions.find(regionName);
    #if defined QUICC_HAS_CUDA_BACKEND
    // sync cuda kernel
    // WARNING
    // this does not allow for kernel overlap
    // to be substituted with cudaEventCreate/Record/Sync/Elapsed
    cudaDeviceSynchronize();
    #endif
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
    if (reg == mRegions.end())
    {
        // not found return object set to zero
        return Tracker::tracking_t{};
    }
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
    DataSet dataCommit = gi.createDataSet<char[10]>("git-commit", DataSpace(dimScalar));

    if (rank == 0)
    {
        dataRanks.write(nRanks);
        dataCommit.write(Framework::gitHash);
    }

    // Write timings
    auto gT = file.createGroup("timings");
    for (auto reg = mRegions.begin(); reg != mRegions.end(); ++reg)
    {
        // compute actual sample size
        auto count = std::get<tracking::count>(reg->second);
        auto sampleSize = (count > collectSize) ? collectSize : count;

        // scalars
        auto gR = gT.createGroup(reg->first);
        DataSet dataCount = gR.createDataSet<std::uint32_t>("count", DataSpace::From(count));
        DataSet dataSampleSize = gR.createDataSet<std::uint32_t>("sampleSize", DataSpace::From(sampleSize));
        DataSet dataMemory = gR.createDataSet<double>("memory", DataSpace::From(dimScalar));
        DataSet dataMemoryDelta = gR.createDataSet<double>("memoryDelta", DataSpace::From(dimScalar));
        if (rank == 0)
        {
            dataCount.write(count);
            dataSampleSize.write(sampleSize);
            dataMemory.write(&std::get<tracking::memory>(reg->second));
            dataMemoryDelta.write(&std::get<tracking::memoryDelta>(reg->second));
        }

        // timings
        std::vector<size_t> dimSample = {nRanks*sampleSize};
        gR.createDataSet<double>("time", DataSpace(dimSample))
            .select({std::size_t(rank)*sampleSize}, {sampleSize})
            .write(std::get<tracking::time>(reg->second).data());

        #if QUICC_MPI
        // wait for everyone to be done
        MPI_Barrier(mComm);
        #endif
    }
#endif
}

void Tracker::print(std::ostream& os)
{

    int rank = 0;
    #ifdef QUICC_MPI
    int nRanks = 1;
    MPI_Comm_rank(mComm, &rank);
    MPI_Comm_size(mComm, &nRanks);
    #endif

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

        #ifdef QUICC_MPI
        // reduce across processes
        MPI_Reduce(data.data(), data_min.data(), stat_size, MPI_DOUBLE, MPI_MIN, 0, mComm);
        data[tracking::time] = timeMax;
        MPI_Reduce(data.data(), data_max.data(), stat_size, MPI_DOUBLE, MPI_MAX, 0, mComm);
        data[tracking::time] = timeAvg;
        MPI_Reduce(data.data(), data_avg.data(), stat_size, MPI_DOUBLE, MPI_SUM, 0, mComm);
        for (auto& d : data_avg)
        {
            d /= nRanks;
        }
        #else
        // serial case
        data_avg[tracking::memory] = data[tracking::memory];
        data_avg[tracking::memoryDelta] = data[tracking::memoryDelta];
        data_min[tracking::time] = timeMin;
        data_max[tracking::time] = timeMax;
        data_avg[tracking::time] = timeAvg;
        #endif

        if (rank == 0)
        {
            os << reg->first
                << '\t' << std::get<tracking::count>(reg->second) << " times\n"
                << "\t stats                                 min             max             avg\n"
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

// static members init
std::map<std::string, Tracker::tracking_t>
    Tracker::mRegions{};

#ifdef QUICC_MPI
MPI_Comm Tracker::mComm = MPI_COMM_NULL;
#endif

} // namespace Profiler
} // namespace QUICC

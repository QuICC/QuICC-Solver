/**
 * @file MpiUtils.hpp
 * @brief Methods for mpi enabled transform
 */
#pragma once

// External includes
//
#include <array>
#include <vector>
#include <mpi.h>

// Project includes
//
namespace QuICC {
namespace Transpose {

constexpr int dimSize = 3;

using point_t = std::array<int, dimSize>;

/// @brief Build send or recv displacement
/// @param absCooNew ending coordinates
/// @param absCooOld starting coordinates
/// @param comm mpi communicator spanning all involved ranks
/// @return send/recv displacement
std::vector<std::vector<int>> getDispls(std::vector<point_t>& absCooNew,
    const std::vector<point_t>& absCooOld, const MPI_Comm comm = MPI_COMM_WORLD) {

    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);

    std::vector<std::vector<int>> sendDispls(ranks);

    std::map<point_t, int> locOldIdx;
    for (std::size_t i = 0; i < absCooOld.size(); ++i) {
        auto&& p = absCooOld[i];
        locOldIdx[p] = i;
    }
    for (int r = 0; r < ranks; ++r) {
        // get new coo from other rank and check if it is here
        std::map<point_t, int> remNewIdx;
        // comm remote coo size
        int remAbsCooNewSize = absCooNew.size();
        MPI_Bcast(&remAbsCooNewSize, 1, MPI_INT, r, comm);
        if (r == rank) {
            MPI_Bcast(absCooNew.data(), absCooNew.size()*dimSize, MPI_INT, r, comm);
            // setup remote map
            for (std::size_t i = 0; i < absCooNew.size(); ++i) {
                auto&& p = absCooNew[i];
                remNewIdx[p] = i;
            }
        }
        else {
            // comm remote coordinates
            std::vector<point_t> remAbsCooNew(remAbsCooNewSize);
            MPI_Bcast(remAbsCooNew.data(), remAbsCooNew.size()*dimSize, MPI_INT, r, comm);
            // setup remote map
            for (std::size_t i = 0; i < remAbsCooNew.size(); ++i) {
                auto&& p = remAbsCooNew[i];
                remNewIdx[p] = i;
            }
        }

        // loop over loc coo to find match
        for (auto itLCoo = locOldIdx.begin(); itLCoo != locOldIdx.end();) {
            auto lCoo = (*itLCoo).first;
            if (auto itRCoo = remNewIdx.find(lCoo); itRCoo != remNewIdx.end()) {
                // std::cout << "rank " << r << ": "
                //     << lCoo[0] << ',' << lCoo[1] << " found " << r << '\n';
                sendDispls[r].push_back((*itLCoo).second);
                itLCoo = locOldIdx.erase(itLCoo);
                remNewIdx.erase(itRCoo);
            }
            else {
                ++itLCoo;
                // std::cout << "not found\n";
            }
        }
    }
    return sendDispls;
}

/// @brief Build set of communicating ranks from displacements
/// This is part of the alltoallw setup
/// @param sendDispls sending displacements
/// @param recvDispls receiving displacements
/// @return set of communicating ranks
std::vector<int> getReducedRanksSet(const std::vector<std::vector<int>>& sendDispls,
    const std::vector<std::vector<int>>& recvDispls) {
    std::set<int> redSet;
    // Save non-empty exchanges
    for (std::size_t i = 0; i < sendDispls.size(); ++i) {
        if (sendDispls[i].size() > 0) {
            redSet.insert(i);
        }
    }
    for (std::size_t i = 0; i < recvDispls.size(); ++i) {
        if (recvDispls[i].size() > 0) {
            redSet.insert(i);
        }
    }
    // Copy to vector
    std::vector<int> res;
    res.reserve(redSet.size());
    for (auto it = redSet.begin(); it != redSet.end(); ++it) {
        res.push_back(*it);
    }
    return res;
}

/// @brief Reduce diplacements to the to the reduced set
/// @param sendDispls
/// @param recvDispls
/// @param redSet
void redDisplsFromSet(std::vector<std::vector<int>>& sendDispls,
    std::vector<std::vector<int>>& recvDispls, std::vector<int>& redSet)
{
    auto redSize = redSet.size();
    std::vector<std::vector<int>> sendDisplsRed(redSize);
    std::vector<std::vector<int>> recvDisplsRed(redSize);

    for (std::size_t r = 0; r < redSize; ++r) {
        sendDisplsRed[r] = std::move(sendDispls[redSet[r]]);
        recvDisplsRed[r] = std::move(recvDispls[redSet[r]]);
    }

    sendDispls = std::move(sendDisplsRed);
    recvDispls = std::move(recvDisplsRed);
}


}
}




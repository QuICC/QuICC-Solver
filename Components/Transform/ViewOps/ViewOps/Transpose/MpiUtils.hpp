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

std::vector<std::vector<int>> getDispls(std::vector<point_t>& absCooNew, const std::vector<point_t>& absCooOld, const MPI_Comm comm = MPI_COMM_WORLD) {

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
                std::cout << "rank " << r << ": "
                    << lCoo[0] << ',' << lCoo[1] << " found " << r << '\n';
                sendDispls[r].push_back((*itLCoo).second);
                itLCoo = locOldIdx.erase(itLCoo);
                remNewIdx.erase(itRCoo);
            }
            else {
                ++itLCoo;
                std::cout << "not found\n";
            }
        }
    }
    return sendDispls;
}


}
}




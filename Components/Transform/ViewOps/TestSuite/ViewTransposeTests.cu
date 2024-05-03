#include <catch2/catch.hpp>

#include "ViewOps/Transpose/Transpose.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"

using namespace QuICC::Memory;
using namespace QuICC::View;

TEST_CASE("Serial DCCSC3D to DCCSC3DJIK", "SerialDCCSC3DtoDCCSC3DJIK")
{
    // FFT out -> AL in, gpu backend
    // Full data and type
    constexpr size_t M = 4;
    constexpr size_t N = 2;
    constexpr size_t K = 2;

    constexpr size_t S = M * N * K;
    std::array<double, S> dataIn = {/*k0*/ 1, 2, 3, 4,
        /*k0*/ 5, 6, 7, 8,
        /*k1*/ 9, 10, 11, 12,
        /*k1*/ 13, 14, 15, 16};

    std::array<double, S> dataOut;

    // view
    constexpr std::uint32_t rank = 3;
    std::array<std::uint32_t, rank> dimensionsIn{M, N, K};
    std::array<std::uint32_t, rank> dimensionsOut{N, M, K};
    // skip setting up pointers and indices
    // they are not used in the serial aka dense tranpose
    std::array<std::vector<std::uint32_t>, rank> pointers = {
        {{}, {}, {}}};
    std::array<std::vector<std::uint32_t>, rank> indices = {
        {{}, {}, {}}};
    using Tin = DCCSC3D;
    using Tout = DCCSC3DJIK;
    View<double, Tin> viewIn(dataIn, dimensionsIn, pointers, indices);
    View<double, Tout> viewOut(dataOut, dimensionsOut, pointers, indices);

    // device mem
    auto memDev = std::make_shared<QuICC::Memory::Cuda::Malloc>();
    QuICC::Memory::MemBlock<double> memBlockIn(S, memDev.get());
    QuICC::Memory::MemBlock<double> memBlockOut(S, memDev.get());

    // set device pointers and indices
    ViewBase<std::uint32_t> pointersDev[rank];
    ViewBase<std::uint32_t> indicesDev[rank];

    // set device views
    View<double, Tin> viewInDev(memBlockIn.data(), memBlockIn.size(), dimensionsIn.data(), pointersDev, indicesDev);
    View<double, Tout> viewOutDev(memBlockOut.data(), memBlockOut.size(), dimensionsOut.data(), pointersDev, indicesDev);

    // cpu -> gpu
    cudaErrChk(cudaMemcpy(viewInDev.data(), viewIn.data(),
        S * sizeof(double), cudaMemcpyHostToDevice));

    // Transpose op
    using namespace QuICC::Transpose::Cuda;
    using namespace QuICC::Transpose;
    auto transposeOp =
      std::make_unique<Op<View<double, Tout>, View<double, Tin>, p201_t>>();

    transposeOp->apply(viewOutDev, viewInDev);

    // gpu -> cpu
    cudaErrChk(cudaMemcpy(viewOut.data(), viewOutDev.data(),
      viewOut.size() * sizeof(double),
      cudaMemcpyDeviceToHost));

    // check
    for (std::uint64_t k = 0; k < K; ++k)
    {
        for (std::uint64_t n = 0; n < N; ++n)
        {
            for (std::uint64_t m = 0; m < M; ++m)
            {
                auto mnk = m + n*M + k*M*N;
                // plane is row major
                auto nkm = n*K + k + m*K*N;
                CHECK(viewIn[mnk] == viewOut[nkm]);
            }
        }
    }
}

TEST_CASE("Serial S1CLCSC3DJIK to DCCSC3DJIK", "SerialS1CLCSC3DJIKtoDCCSC3DJIK")
{
    // AL out -> JW in, gpu backend
    // Full data and type
    constexpr size_t M = 4;
    constexpr size_t N = 2;
    constexpr size_t K = 3;

    constexpr size_t S = (M + (M-1) + (M-2)) * N ;
    std::array<double, S> dataIn = {
        /*k0*/ 1, 5,
        /*k0*/ 2, 6,
        /*k0*/ 3, 7,
        /*k0*/ 4, 8,
        /*k1*/ 9, 12,
        /*k1*/ 10, 13,
        /*k1*/ 11, 14,
        /*k2*/ 15, 17,
        /*k2*/ 16, 18};

    // N K M
    std::array<double, S> dataRef = {
        /*m0*/ 1, 9, 15,
        /*m0*/ 5, 12, 17,
        /*m1*/ 2, 10, 16,
        /*m1*/ 6, 13, 18,
        /*m2*/ 3, 11,
        /*m2*/ 7, 14,
        /*m3*/ 4,
        /*m3*/ 8};

    std::array<double, S> dataOut;

    // view
    constexpr std::uint32_t rank = 3;
    std::array<std::uint32_t, rank> dimensionsIn{M, N, K};
    std::array<std::uint32_t, rank> dimensionsOut{N, M, K};
    // skip setting up pointers and indices
    // they are not used in the serial aka dense tranpose
    std::array<std::vector<std::uint32_t>, rank> pointers = {
        {{}, {}, {}}};
    std::array<std::vector<std::uint32_t>, rank> indices = {
        {{}, {}, {}}};

    using Tin = S1CLCSC3DJIK;
    using Tout = DCCSC3DJIK;

    View<double, Tin> viewIn(dataIn, dimensionsIn, pointers, indices);
    View<double, Tout> viewOut(dataOut, dimensionsOut, pointers, indices);

    // device mem
    auto memDev = std::make_shared<QuICC::Memory::Cuda::Malloc>();
    QuICC::Memory::MemBlock<double> memBlockIn(S, memDev.get());
    QuICC::Memory::MemBlock<double> memBlockOut(S, memDev.get());

    // set device pointers and indices
    ViewBase<std::uint32_t> pointersDev[rank];
    ViewBase<std::uint32_t> indicesDev[rank];

    // set device views
    View<double, Tin> viewInDev(memBlockIn.data(), memBlockIn.size(), dimensionsIn.data(), pointersDev, indicesDev);
    View<double, Tout> viewOutDev(memBlockOut.data(), memBlockOut.size(), dimensionsOut.data(), pointersDev, indicesDev);

    // cpu -> gpu
    cudaErrChk(cudaMemcpy(viewInDev.data(), viewIn.data(),
        S * sizeof(double), cudaMemcpyHostToDevice));

    // Transpose op
    using namespace QuICC::Transpose::Cuda;
    using namespace QuICC::Transpose;
    auto transposeOp =
      std::make_unique<Op<View<double, Tout>, View<double, Tin>, p201_t>>();

    transposeOp->apply(viewOutDev, viewInDev);

    // gpu -> cpu
    cudaErrChk(cudaMemcpy(viewOut.data(), viewOutDev.data(),
      viewOut.size() * sizeof(double),
      cudaMemcpyDeviceToHost));

    // check
    for (std::uint64_t s = 0; s < S; ++s)
    {
        CHECK(dataRef[s] == viewOut[s]);
    }
}

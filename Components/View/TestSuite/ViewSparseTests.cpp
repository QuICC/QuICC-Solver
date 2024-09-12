#include <catch2/catch.hpp>

#include "View/ViewSparse.hpp"

using namespace QuICC::View;

TEST_CASE("ViewOneDimSparse", "[ViewOneDimSparse]")
{
    constexpr size_t SF = 10;
    std::array<double, SF> fullData = {0,0,3,4,5,0,7,0,0,10};
    constexpr size_t S = 5;
    std::array<double, S> data = {3,4,5,7,10};
    std::array<std::uint32_t, 1> dimensions {SF};
    std::array<std::vector<std::uint32_t>, 1> pointers = {{{0,5}}};
    std::array<std::vector<std::uint32_t>, 1> indices = {{{2,3,4,6,9}}};

    using sparse1D = DimLevelType<compressed_t>;
    View<double, Attributes<sparse1D>> someView (data, dimensions, pointers, indices);

    CHECK(someView.rank() == 1);
    CHECK(someView.dims()[0] == SF);
    CHECK(someView.lds() == 0);

    // check data
    for (std::size_t i = 0; i < S; ++i)
    {
        CHECK(data[i] == someView.data()[i]);
    }

    for (std::size_t i = 0; i < S; ++i)
    {
        CHECK(fullData[indices[0][i]] == someView(indices[0][i]));
    }

    CHECK_THROWS(fullData[0] == someView(0));
}



TEST_CASE("ViewCopy", "[ViewCopy]")
{
    constexpr size_t SF = 10;
    std::array<double, SF> fullData = {0,0,3,4,5,0,7,0,0,10};
    constexpr size_t S = 5;
    std::array<double, S> data = {3,4,5,7,10};
    std::array<std::uint32_t, 1> dimensions {SF};
    std::array<std::vector<std::uint32_t>, 1> pointers = {{{0,5}}};
    std::array<std::vector<std::uint32_t>, 1> indices = {{{2,3,4,6,9}}};

    using sparse1D = DimLevelType<compressed_t>;
    // empty view
    View<double, Attributes<sparse1D>> someView;

    // move ctor
    someView = View<double, Attributes<sparse1D>> (data, dimensions, pointers, indices);

    CHECK(someView.rank() == 1);
    CHECK(someView.dims()[0] == SF);
    CHECK(someView.lds() == 0);

    // check data
    for (std::size_t i = 0; i < S; ++i)
    {
        CHECK(data[i] == someView.data()[i]);
    }

    for (std::size_t i = 0; i < S; ++i)
    {
        CHECK(fullData[indices[0][i]] == someView(indices[0][i]));
    }

    CHECK_THROWS(fullData[0] == someView(0));

}

TEST_CASE("View Two Diminsional CSR", "[ViewTwoDimCSR]")
{
    constexpr size_t M = 3;
    constexpr size_t N = 4;
    std::array<double, M*N> fullData = {1,0,0,3,
                                        0,0,0,0,
                                        2,0,0,0};
    constexpr size_t S = 3;
    std::array<double, S> data = {1,3,2};

    using CSR = DimLevelType<sparse_t, compressed_t>;
    std::array<std::uint32_t, 2> dimensions {M, N};
    std::array<std::vector<std::uint32_t>, 2> pointers = {{ {}, {0,2,2,3}}};
    std::array<std::vector<std::uint32_t>, 2> indices = {{ {}, {0,3,0}}};
    View<double, Attributes<CSR>> someView (data, dimensions, pointers, indices);

    CHECK(someView.rank() == 2);
    CHECK(someView.dims()[0] == M);
    CHECK(someView.dims()[1] == N);
    CHECK(someView.lds() == 0);

    for (std::size_t l = 0; l < S; ++l)
    {
        CHECK(data[l] == someView.data()[l]);
    }

    CHECK(fullData[0] == someView(0, 0));
    CHECK(fullData[3] == someView(0, 3));
    CHECK(fullData[8] == someView(2, 0));

    CHECK_THROWS(fullData[1] == someView(1, 0));
}

TEST_CASE("View Two Dimensional CSC", "[ViewTwoDimCSC]")
{
    constexpr size_t M = 3;
    constexpr size_t N = 4;
    std::array<double, M*N> fullData = {1,0,0,3,
                                        0,0,0,0,
                                        2,0,0,0};
    constexpr size_t S = 3;
    std::array<double, S> data = {1,2,3};

    std::array<std::uint32_t, 2> dimensions {M, N};
    std::array<std::vector<std::uint32_t>, 2> pointers = {{ {0,2,2,2,3}, {}}};
    std::array<std::vector<std::uint32_t>, 2> indices = {{ {0,2,0}, {}}};
    View<double, CSC> someView (data, dimensions, pointers, indices);

    CHECK(someView.rank() == 2);
    CHECK(someView.dims()[0] == M);
    CHECK(someView.dims()[1] == N);
    CHECK(someView.lds() == 0);

    for (std::size_t l = 0; l < S; ++l)
    {
        CHECK(data[l] == someView.data()[l]);
    }

    CHECK(fullData[0] == someView(0, 0));
    CHECK(fullData[3] == someView(0, 3));
    CHECK(fullData[8] == someView(2, 0));

    CHECK_THROWS(fullData[1] == someView(1, 0));
}


TEST_CASE("ViewThreeDimDtrClColMaj", "[ViewThreeDimDtrClColMaj]")
{
    constexpr size_t M = 3;
    constexpr size_t N = 2;
    constexpr size_t K = 3;
    constexpr size_t SF = M*N*K;
    std::array<double, SF> fullData = {1,2,3,
                                       4,5,6,
                                            0,7,8,
                                            0,10,11,
                                                   0,0,13,
                                                   0,0,16};

    constexpr size_t S = M*2+2*2+2;
    std::array<double, S> data = {1,2,3,
                                  4,5,6,
                                        7,8,
                                        10,11,
                                              13,
                                              16};


    // Akin full AL op
    // Dense triangular column/layer
    using DTRCL3D = Attributes<DimLevelType<step1K_t, dense_t, dense_t>>;
    std::array<std::uint32_t, 3> dimensions {M, N, K};
    std::array<std::vector<std::uint32_t>, 3> pointers = {{{}, {}, {}}};
    std::array<std::vector<std::uint32_t>, 3> indices = {{{}, {}, {}}};
    View<double, DTRCL3D> someView (data, dimensions, pointers, indices);

    CHECK(someView.rank() == 3);
    CHECK(someView.dims()[0] == M);
    CHECK(someView.dims()[1] == N);
    CHECK(someView.dims()[2] == K);

    for (std::size_t l = 0; l < S; ++l)
    {
        CHECK(data[l] == someView.data()[l]);
    }

    CHECK(fullData[0] == someView(0, 0, 0));
    CHECK(fullData[2] == someView(2, 0, 0));
    CHECK(fullData[1+M*1] == someView(1, 1, 0));
    CHECK(fullData[1+M*1+M*N*1] == someView(1, 1, 1));
    CHECK(fullData[2+M*1+M*N*2] == someView(2, 1, 2));

    CHECK_THROWS(fullData[2+M*1+M*N*1] == someView(0, 1, 1));
    CHECK_THROWS(fullData[0+M*1+M*N*2] == someView(0, 1, 2));
    CHECK_THROWS(fullData[0+M*1+M*N*2] == someView(1, 1, 2));
}

TEST_CASE("ViewThreeDim Compressed step 1 Column Layer ColMaj", "[ViewThreeDimCtrClColMaj]")
{
    constexpr size_t M = 4;
    constexpr size_t N = 2;
    constexpr size_t K = 4;
    constexpr size_t SF = M*N*K;
    std::array<double, SF> fullData = {1,2,3,4,
                                       5,6,7,8,
                                         0,0,0,0,
                                         0,0,0,0,
                                           0,0,0,0,
                                           0,0,0,0,
                                             0,0,0,9,
                                             0,0,0,10};

    constexpr size_t S = M*N+N;
    std::array<double, S> data = {1,2,3,4,
                                  5,6,7,8,
                                        9,
                                        10};


    // Akin distributed AL op
    // Dense triangular column, dense row, compressed triangular layer
    using CTRCL3D = Attributes<DimLevelType<step1K_t, dense_t, compressed_t>>;
    std::array<std::uint32_t, 3> dimensions {M, N, K};
    std::array<std::vector<std::uint32_t>, 3> pointers = {{{}, {}, {0,2}}};
    std::array<std::vector<std::uint32_t>, 3> indices = {{{}, {}, {0,3}}};
    View<double, CTRCL3D> someView (data, dimensions, pointers, indices);

    CHECK(someView.rank() == 3);
    CHECK(someView.dims()[0] == M);
    CHECK(someView.dims()[1] == N);
    CHECK(someView.dims()[2] == K);

    for (std::size_t l = 0; l < S; ++l)
    {
        CHECK(data[l] == someView.data()[l]);
    }

    CHECK(fullData[0+M*0+M*N*0] == someView(0, 0, 0));
    CHECK(fullData[2+M*0+M*N*0] == someView(2, 0, 0));
    CHECK(fullData[1+M*1+M*N*0] == someView(1, 1, 0));
    CHECK(fullData[3+M*1+M*N*3] == someView(3, 1, 3));

    CHECK_THROWS(fullData[2+M*1+M*N*1] == someView(2, 1, 1));
    CHECK_THROWS(fullData[1+M*1+M*N*1] == someView(1, 1, 1));
    CHECK_THROWS(fullData[0+M*1+M*N*3] == someView(0, 1, 3));
}

TEST_CASE("ViewThreeDim Compressed step 1 Column Layer, Row Column Layer layout",
    "[ViewThreeDimCtrClJIK]")
{
    constexpr size_t M = 4;
    constexpr size_t N = 2;
    constexpr size_t K = 4;
    constexpr size_t SF = M*N*K;
    std::array<double, SF> fullData = {1,5,
                                       2,6,
                                       3,7,
                                       4,8,
                                         0,0,
                                         0,0,
                                         0,0,
                                         0,0,
                                           0,0,
                                           0,0,
                                           0,0,
                                           0,0,
                                             9,10,
                                             0,0,
                                             0,0,
                                             0,0};

    constexpr size_t S = M*N+M;
    std::array<double, S> data = {1,5,
                                  2,6,
                                  3,7,
                                  4,8,
                                      9,10};


    // Akin distributed AL op
    // Dense triangular column, dense row, compressed triangular layer
    // row, column, layer layout in memory
    using CTRCL3D = DimLevelType<step1K_t, dense_t, compressed_t>;
    using JIK = LoopOrderType<j_t, i_t, k_t>;
    std::array<std::uint32_t, 3> dimensions {M, N, K};
    std::array<std::vector<std::uint32_t>, 3> pointers = {{{}, {}, {0,2}}};
    std::array<std::vector<std::uint32_t>, 3> indices = {{{}, {}, {0,3}}};
    View<double, Attributes<CTRCL3D, JIK>> someView (data, dimensions, pointers, indices);

    CHECK(someView.rank() == 3);
    CHECK(someView.dims()[0] == M);
    CHECK(someView.dims()[1] == N);
    CHECK(someView.dims()[2] == K);

    for (std::size_t l = 0; l < S; ++l)
    {
        CHECK(data[l] == someView.data()[l]);
    }

    CHECK(fullData[0] == someView(0, 0, 0));
    CHECK(fullData[2*N] == someView(2, 0, 0));
    CHECK(fullData[1*N+1] == someView(1, 1, 0));
    CHECK(fullData[0*N+1+M*N*3] == someView(0, 1, 3));

    CHECK_THROWS(fullData[2*N+1+M*N*1] == someView(2, 1, 1));
}

TEST_CASE("ViewThreeDim Dense Column, Compressed step 1 row/layer ColMaj ", "[ViewThreeDimCtrRlColMaj]")
{
    constexpr size_t M = 2;
    constexpr size_t N = 4;
    constexpr size_t K = 4;
    constexpr size_t SF = M*N*K;
    std::array<double, SF> fullData = {1,2,
                                       3,4,
                                       5,6,
                                       7,8,
                                         0,0,
                                         0,0,
                                         0,0,
                                         0,0,
                                           0,0,
                                           0,0,
                                           0,0,
                                           0,0,
                                             9,10,
                                             0,0,
                                             0,0,
                                             0,0};

    constexpr size_t S = M*N+M;
    std::array<double, S> data = {1,2,
                                  3,4,
                                  5,6,
                                  7,8,
                                    9,10};

    // Akin distributed AL op
    // Dense column, dense triangular row, compressed triangular layer
    std::array<std::uint32_t, 3> dimensions {M, N, K};
    std::array<std::vector<std::uint32_t>, 3> pointers = {{{}, {}, {0,2}}};
    std::array<std::vector<std::uint32_t>, 3> indices = {{{}, {}, {0,3}}};
    View<double, CS1RL3D> someView (data, dimensions, pointers, indices);

    CHECK(someView.rank() == 3);
    CHECK(someView.dims()[0] == M);
    CHECK(someView.dims()[1] == N);
    CHECK(someView.dims()[2] == K);

    for (std::size_t l = 0; l < S; ++l)
    {
        CHECK(data[l] == someView.data()[l]);
    }

    CHECK(fullData[0] == someView(0, 0, 0));
    CHECK(fullData[0+M*1+0] == someView(0, 1, 0));
    CHECK(fullData[1+M*1+0] == someView(1, 1, 0));
    CHECK(fullData[1+M*N*3] == someView(1, 0, 3));

    CHECK_THROWS(fullData[0+M*1+M*N*3] == someView(0, 1, 3));
    CHECK_THROWS(fullData[2+M*1+M*N*1] == someView(1, 1, 1));
}

TEST_CASE("ViewThreeDim Dense Column, Compressed step 1 row/layer, Row Column Layer layout ",
    "[ViewThreeDimCtrRlJIK]")
{
    constexpr size_t M = 2;
    constexpr size_t N = 4;
    constexpr size_t K = 4;
    constexpr size_t SF = M*N*K;
    std::array<double, SF> fullData = {1,3,5,7,
                                       2,4,6,8,
                                         0,0,0,0,
                                         0,0,0,0,
                                           0,0,0,0,
                                           0,0,0,0,
                                             9,0,0,0,
                                             10,0,0,0};

    constexpr size_t S = M*N+M;
    std::array<double, S> data = {1,3,5,7,
                                  2,4,6,8,
                                      9,
                                      10};


    // Akin distributed AL op
    // Dense column, dense triangular row, compressed triangular layer
    std::array<std::uint32_t, 3> dimensions {M, N, K};
    std::array<std::vector<std::uint32_t>, 3> pointers = {{{}, {}, {0,2}}};
    std::array<std::vector<std::uint32_t>, 3> indices = {{{}, {}, {0,3}}};
    View<double, CS1RL3DJIK> someView (data, dimensions, pointers, indices);

    CHECK(someView.rank() == 3);
    CHECK(someView.dims()[0] == M);
    CHECK(someView.dims()[1] == N);
    CHECK(someView.dims()[2] == K);
    CHECK(someView.lds() == 0);

    for (std::size_t l = 0; l < S; ++l)
    {
        CHECK(data[l] == someView.data()[l]);
    }

    CHECK(fullData[0] == someView(0, 0, 0));
    CHECK(fullData[1+0+0] == someView(0, 1, 0));
    CHECK(fullData[1*N+1+0] == someView(1, 1, 0));
    CHECK(fullData[1*N+M*N*3] == someView(1, 0, 3));

    CHECK_THROWS(fullData[0+1+M*N*3] == someView(0, 1, 3));
    CHECK_THROWS(fullData[1*N+1+M*N*1] == someView(1, 1, 1));
}

TEST_CASE("ViewThreeDim Dense Column, Compressed step 1 row/layer, Row Column Layer layout #2",
    "[ViewThreeDimCtrRlJIK2]")
{
    constexpr size_t M = 2;
    constexpr size_t N = 4;
    constexpr size_t K = 4;
    constexpr size_t SF = M*N*K;
    std::array<double, SF> fullData = {0,0,0,0,
                                       0,0,0,0,
                                         1,3,5,0,
                                         2,4,6,0,
                                           0,0,0,0,
                                           0,0,0,0,
                                             7,0,0,0,
                                             8,0,0,0};

    constexpr size_t S = M*(N-1)+M;
    std::array<double, S> data = {1,3,5,
                                  2,4,6,
                                      7,
                                      8};


    // Akin distributed AL op
    // Dense column, dense triangular row, compressed triangular layer
    std::array<std::uint32_t, 3> dimensions {M, N, K};
    std::array<std::vector<std::uint32_t>, 3> pointers = {{{}, {}, {0,2}}};
    std::array<std::vector<std::uint32_t>, 3> indices = {{{}, {}, {1,3}}};
    View<double, CS1RL3DJIK> someView (data, dimensions, pointers, indices);

    CHECK(someView.rank() == 3);
    CHECK(someView.dims()[0] == M);
    CHECK(someView.dims()[1] == N);
    CHECK(someView.dims()[2] == K);
    CHECK(someView.lds() == 0);

    for (std::size_t l = 0; l < S; ++l)
    {
        CHECK(data[l] == someView.data()[l]);
    }

    CHECK(fullData[0+M*N] == someView(0, 0, 1));
    CHECK(fullData[1+0+M*N] == someView(0, 1, 1));
    CHECK(fullData[1*N+1+M*N] == someView(1, 1, 1));
    CHECK(fullData[1*N+M*N*3] == someView(1, 0, 3));

    CHECK_THROWS(fullData[0+1+M*N*3] == someView(0, 1, 3));
    CHECK_THROWS(fullData[1*N+1+M*N*1] == someView(1, 1, 0));
}


TEST_CASE("ViewThreeDim Dense Column CSC Row/Layer ColMaj", "[ViewThreeDimDCCSCColMaj]")
{
    constexpr size_t M = 4;
    constexpr size_t N = 2;
    constexpr size_t K = 4;
    constexpr size_t SF = M*N*K;
    std::array<double, SF> fullData = {1,2,3,4,
                                       0,0,0,0,
                                         5,6,7,8,
                                         0,0,0,0,
                                           0,0,0,0,
                                           0,0,0,0,
                                             9,10,11,12,
                                             13,14,15,16};

    constexpr size_t S = M*N*2;
    std::array<double, S> data = {1,2,3,4,
                                  5,6,7,8,
                                        9,10,11,12,
                                        13,14,15,16};


    // Akin distributed Fourier op
    // 2D CSC column major matrix (N,K) which elements are 1D dense vectors,
    // i.e a 3D tensor (M,N,K) with fully populated columns
    std::array<std::uint32_t, 3> dimensions {M, N, K};
    std::array<std::vector<std::uint32_t>, 3> pointers = {{{}, {0,1,2,2,4}, {}}};
    std::array<std::vector<std::uint32_t>, 3> indices = {{{}, {0,0,0,1}, {}}};
    View<double, DCCSC3D> someView (data, dimensions, pointers, indices);

    CHECK(someView.rank() == 3);
    CHECK(someView.dims()[0] == M);
    CHECK(someView.dims()[1] == N);
    CHECK(someView.dims()[2] == K);
    CHECK(someView.lds() == M);

    for (std::size_t l = 0; l < S; ++l)
    {
        CHECK(data[l] == someView.data()[l]);
    }

    CHECK(fullData[0] == someView(0, 0, 0));
    CHECK(fullData[2] == someView(2, 0, 0));
    CHECK(fullData[1+M*0+M*N*1] == someView(1, 0, 1));
    CHECK(fullData[0+M*1+M*N*3] == someView(0, 1, 3));

    CHECK_THROWS(fullData[0+M*1] == someView(0, 1, 0));
}

TEST_CASE("ViewThreeDim Dense Column CSC Row/Layer ColMaj Padded", "[ViewThreeDimDCCSCColMajPadded]")
{
    constexpr size_t M = 4;
    constexpr size_t N = 2;
    constexpr size_t K = 4;
    constexpr size_t lds = 6;
    constexpr size_t SF = M*N*K;
    std::array<double, SF> fullData = {1,2,3,4,
                                       0,0,0,0,
                                         5,6,7,8,
                                         0,0,0,0,
                                           0,0,0,0,
                                           0,0,0,0,
                                             9,10,11,12,
                                             13,14,15,16};

    constexpr size_t S = lds*N*2;
    std::array<double, S> data = {1,2,3,4,0,0,
                                  5,6,7,8,0,0,
                                        9,10,11,12,0,0,
                                        13,14,15,16,0,0};


    // Akin distributed Fourier op
    // 2D CSC column major matrix (N,K) which elements are 1D dense vectors,
    // i.e a 3D tensor (M,N,K) with fully populated columns
    std::array<std::uint32_t, 3> dimensions {M, N, K};
    std::array<std::vector<std::uint32_t>, 3> pointers = {{{}, {0,1,2,2,4}, {}}};
    std::array<std::vector<std::uint32_t>, 3> indices = {{{}, {0,0,0,1}, {}}};
    View<double, DCCSC3D> someView (data, dimensions, pointers, indices, lds);

    CHECK(someView.rank() == 3);
    CHECK(someView.dims()[0] == M);
    CHECK(someView.dims()[1] == N);
    CHECK(someView.dims()[2] == K);
    CHECK(someView.lds() == lds);

    for (std::size_t l = 0; l < S; ++l)
    {
        CHECK(data[l] == someView.data()[l]);
    }

    CHECK(fullData[0] == someView(0, 0, 0));
    CHECK(fullData[2] == someView(2, 0, 0));
    CHECK(fullData[1+M*0+M*N*1] == someView(1, 0, 1));
    CHECK(fullData[0+M*1+M*N*3] == someView(0, 1, 3));

    CHECK_THROWS(fullData[0+M*1] == someView(0, 1, 0));
}

TEST_CASE("ViewThreeDim Dense Column CSC Row/Layer JIK", "[ViewThreeDimDCCSCJIK]")
{
    constexpr size_t M = 4;
    constexpr size_t N = 2;
    constexpr size_t K = 4;
    constexpr size_t SF = M*N*K;
    std::array<double, SF> fullData = {1,0,
                                       2,0,
                                       3,0,
                                       4,0,
                                         5,0,
                                         6,0,
                                         7,0,
                                         8,0,
                                           0,0,
                                           0,0,
                                           0,0,
                                           0,0,
                                             9,13,
                                             10,14,
                                             11,15,
                                             12,16};

    constexpr size_t S = M*N*2;
    std::array<double, S> data = {1,
                                  2,
                                  3,
                                  4,
                                    5,
                                    6,
                                    7,
                                    8,
                                      9,13,
                                      10,14,
                                      11,15,
                                      12,16};


    // Akin distributed Fourier op
    // 2D CSC column major matrix (N,K) which elements are 1D dense vectors,
    // i.e a 3D tensor (M,N,K) with fully populated columns
    std::array<std::uint32_t, 3> dimensions {M, N, K};
    std::array<std::vector<std::uint32_t>, 3> pointers = {{{}, {0,1,2,2,4}, {}}};
    std::array<std::vector<std::uint32_t>, 3> indices = {{{}, {0,0,0,1}, {}}};
    View<double, DCCSC3DJIK> someView (data, dimensions, pointers, indices);

    CHECK(someView.rank() == 3);
    CHECK(someView.dims()[0] == M);
    CHECK(someView.dims()[1] == N);
    CHECK(someView.dims()[2] == K);

    for (std::size_t l = 0; l < S; ++l)
    {
        CHECK(data[l] == someView.data()[l]);
    }

    CHECK(fullData[0] == someView(0, 0, 0));
    CHECK(fullData[N*2] == someView(2, 0, 0));
    CHECK(fullData[N*1+0+M*N*1] == someView(1, 0, 1));
    CHECK(fullData[0+1+M*N*3] == someView(0, 1, 3));

    CHECK_THROWS(fullData[0+M*1] == someView(0, 1, 0));
}

TEST_CASE("ViewThreeDim Dense Column/Row compressed layer ColMaj", "[ViewThreeDimCsl3DColMaj]")
{
    constexpr size_t M = 4;
    constexpr size_t N = 2;
    constexpr size_t K = 4;
    constexpr size_t SF = M*N*K;
    std::array<double, SF> fullData = {1,2,3,4,
                                       5,6,7,8,
                                         0,0,0,0,
                                         0,0,0,0,
                                           0,0,0,0,
                                           0,0,0,0,
                                             9,10,11,12,
                                             13,14,15,16};

    constexpr size_t S = M*N*2;
    std::array<double, S> data = {1,2,3,4,
                                  5,6,7,8,
                                        9,10,11,12,
                                        13,14,15,16};


    // Akin distributed JW projector operator
    // 1D compressed array with 2D dense elements
    // i.e a 3D tensor (M,N,K) with fully populated layers
    std::array<std::uint32_t, 3> dimensions {M, N, K};
    std::array<std::vector<std::uint32_t>, 3> pointers = {{{}, {}, {0,2}}};
    std::array<std::vector<std::uint32_t>, 3> indices = {{{}, {}, {0,3}}};
    View<double, CSL3D> someView (data, dimensions, pointers, indices);

    CHECK(someView.rank() == 3);
    CHECK(someView.dims()[0] == M);
    CHECK(someView.dims()[1] == N);
    CHECK(someView.dims()[2] == K);

    for (std::size_t l = 0; l < S; ++l)
    {
        CHECK(data[l] == someView.data()[l]);
    }

    CHECK(fullData[0] == someView(0, 0, 0));
    CHECK(fullData[2] == someView(2, 0, 0));
    CHECK(fullData[1+M*1+0] == someView(1, 1, 0));
    CHECK(fullData[2+M*1+M*N*3] == someView(2, 1, 3));

    CHECK_THROWS(fullData[0+M*1+M*N*1] == someView(0, 1, 1));
}

TEST_CASE("ViewThreeDim Dense Column/Row compressed layer Row Major", "[ViewThreeDimCsl3DRowMaj]")
{
    constexpr size_t M = 4;
    constexpr size_t N = 2;
    constexpr size_t K = 4;
    constexpr size_t SF = M*N*K;
    std::array<double, SF> fullData = {1,2,3,4,
                                       5,6,7,8,
                                         0,0,0,0,
                                         0,0,0,0,
                                           0,0,0,0,
                                           0,0,0,0,
                                             9,10,11,12,
                                             13,14,15,16};


    constexpr size_t S = M*N*2;
    std::array<double, S> data = {1,5,
                                  2,6,
                                  3,7,
                                  4,8,
                                        9,13,
                                        10,14,
                                        11,15,
                                        12,16};


    // Akin distributed JW projector operator
    // 1D compressed array with 2D dense elements
    // i.e a 3D tensor (M,N,K) with fully populated layers
    std::array<std::uint32_t, 3> dimensions {M, N, K};
    std::array<std::vector<std::uint32_t>, 3> pointers = {{{}, {}, {0,2}}};
    std::array<std::vector<std::uint32_t>, 3> indices = {{{}, {}, {0,3}}};
    View<double, CSL3DJIK> someView (data, dimensions, pointers, indices);

    CHECK(someView.rank() == 3);
    CHECK(someView.dims()[0] == M);
    CHECK(someView.dims()[1] == N);
    CHECK(someView.dims()[2] == K);

    for (std::size_t l = 0; l < S; ++l)
    {
        CHECK(data[l] == someView.data()[l]);
    }

    CHECK(fullData[0] == someView(0, 0, 0));
    CHECK(fullData[2] == someView(2, 0, 0));
    CHECK(fullData[1+M*1+0] == someView(1, 1, 0));
    CHECK(fullData[2+M*1+M*N*3] == someView(2, 1, 3));

    CHECK_THROWS(fullData[0+M*1+M*N*1] == someView(0, 1, 1));
}

TEST_CASE("ViewThreeDim Triangular Column/Layer CSC in NK plane ColMaj", "[ViewThreeDimTrClCSCS3DColMaj]")
{
    constexpr size_t M = 4;
    constexpr size_t N = 3;
    constexpr size_t K = 4;
    constexpr size_t SF = M*N*K;
    std::array<double, SF> fullData = {1,2,3,4,
                                       0,0,0,0,
                                       5,6,7,8,
                                         0,0,0,0,
                                         0,0,0,0,
                                         0,0,0,0,
                                           0,0,0,0,
                                           0,0,0,0,
                                           0,0,0,0,
                                             9,0,0,0,
                                             0,0,0,0,
                                             0,0,0,0};

    constexpr size_t S = M*2+1;
    std::array<double, S> data = {1,2,3,4,
                                  5,6,7,8,
                                        9};


    // Akin distributed AL projector input on cpu
    // Step 1 column layer, with (N,K) plane a 2D CSC column major matrix
    std::array<std::uint32_t, 3> dimensions {M, N, K};
    std::array<std::vector<std::uint32_t>, 3> pointers = {{{}, {0,2,2,2,3}, {}}};
    std::array<std::vector<std::uint32_t>, 3> indices = {{{}, {0,2,0}, {}}};
    View<double, S1CLCSC3D> someView (data, dimensions, pointers, indices);

    CHECK(someView.rank() == 3);
    CHECK(someView.dims()[0] == M);
    CHECK(someView.dims()[1] == N);
    CHECK(someView.dims()[2] == K);

    for (std::size_t l = 0; l < S; ++l)
    {
        CHECK(data[l] == someView.data()[l]);
    }

    CHECK(fullData[2] == someView(2, 0, 0));
    CHECK(fullData[3+M*2] == someView(3, 2, 0));
    CHECK(fullData[0+M*N*3] == someView(0, 0, 3));

    CHECK_THROWS(fullData[0+M*1+M*N*1] == someView(0, 1, 1));
    CHECK_THROWS(fullData[0+M*1+M*N*1] == someView(1, 0, 3));
}


TEST_CASE("ViewThreeDim Triangular Column/Layer CSC in NK plane JIK", "[ViewThreeDimTrClCSCS3DJIK]")
{
    constexpr size_t M = 4;
    constexpr size_t N = 3;
    constexpr size_t K = 4;
    constexpr size_t SF = M*N*K;
    std::array<double, SF> fullData = {1,0,5,
                                       2,0,6,
                                       3,0,7,
                                       4,0,8,
                                         0,0,0,
                                         0,0,0,
                                         0,0,0,
                                         0,0,0,
                                           0,0,0,
                                           0,0,0,
                                           0,0,0,
                                           0,0,0,
                                             9,0,0,
                                             0,0,0,
                                             0,0,0,
                                             0,0,0};

    constexpr size_t S = M*2+1;
    std::array<double, S> data = {1,5,
                                  2,6,
                                  3,7,
                                  4,8,
                                    9};


    // Akin distributed AL projector input on cpu
    // Step 1 column layer, with (N,K) plane a 2D CSC column major matrix
    std::array<std::uint32_t, 3> dimensions {M, N, K};
    std::array<std::vector<std::uint32_t>, 3> pointers = {{{}, {0,2,2,2,3}, {}}};
    std::array<std::vector<std::uint32_t>, 3> indices = {{{}, {0,2,0}, {}}};
    View<double, S1CLCSC3DJIK> someView (data, dimensions, pointers, indices);

    CHECK(someView.rank() == 3);
    CHECK(someView.dims()[0] == M);
    CHECK(someView.dims()[1] == N);
    CHECK(someView.dims()[2] == K);

    for (std::size_t l = 0; l < S; ++l)
    {
        CHECK(data[l] == someView.data()[l]);
    }

    CHECK(fullData[2*N] == someView(2, 0, 0));
    CHECK(fullData[3*N+2] == someView(3, 2, 0));
    CHECK(fullData[M*N*3] == someView(0, 0, 3));

    CHECK_THROWS(fullData[0+M*1+M*N*1] == someView(0, 1, 1));
    CHECK_THROWS(fullData[0+M*1+M*N*1] == someView(1, 0, 3));
}
#include <catch2/catch.hpp>

#include "View/ViewDense.hpp"

TEST_CASE("ViewOneDimDense", "[ViewOneDimDense]")
{
    using namespace QuICC::Memory;

    constexpr size_t S = 5;
    std::array<double, S> data = {1,2,3,4,5};

    using dense1D = DimLevelType<dense_t>;
    std::array<std::uint32_t, 1> dimensions {S};
    View<double, Attributes<dense1D>> someView (data , dimensions);

    CHECK(someView.rank() == 1);
    CHECK(someView.dims()[0] == S);

    for (std::size_t i = 0; i < S; ++i)
    {
        CHECK(data[i] == someView(i));
        CHECK(data[i] == someView.data()[i]);
        CHECK(data[i] == (&someView(0))[i]);
    }

}

TEST_CASE("ViewTwoDimDenseColMaj", "[ViewTwoDimDenseColMaj]")
{
    using namespace QuICC::Memory;

    constexpr size_t S = 3*2;
    std::array<double, S> data = {1,2,3,4,5,6};

    using dense2D = DimLevelType<dense_t, dense_t>;
    std::array<std::uint32_t, 2> dimensions {3, 2};
    View<double, Attributes<dense2D>> someView (data, dimensions);

    CHECK(someView.rank() == 2);
    CHECK(someView.dims()[0] == 3);
    CHECK(someView.dims()[1] == 2);


    for (std::size_t l = 0; l < S; ++l)
    {
        CHECK(data[l] == someView.data()[l]);
    }

    CHECK(data[0] == someView(0, 0));
    CHECK(data[2] == someView(2, 0));
    CHECK(data[4] == someView(1, 1));
}

TEST_CASE("ViewTwoDimDenseRowMaj", "[ViewTwoDimDenseRowMaj]")
{
    using namespace QuICC::Memory;

    constexpr size_t S = 3*2;
    std::array<double, S> data = {1,2,3,4,5,6};

    using dense2D = DimLevelType<dense_t, dense_t>;
    using rowMaj = LoopOrderType<j_t, i_t>;
    std::array<std::uint32_t, 2> dimensions {3, 2};
    View<double, Attributes<dense2D, rowMaj>> someView (data, dimensions);

    CHECK(someView.rank() == 2);
    CHECK(someView.dims()[0] == 3);
    CHECK(someView.dims()[1] == 2);

    for (std::size_t l = 0; l < S; ++l)
    {
        CHECK(data[l] == someView.data()[l]);
    }

    CHECK(data[0] == someView(0, 0));
    CHECK(data[2] == someView(1, 0));
    CHECK(data[3] == someView(1, 1));
}

TEST_CASE("ViewThreeDimDenseLayMaj", "[ViewThreeDimDenseLayMaj]")
{
    using namespace QuICC::Memory;

    constexpr size_t M = 3;
    constexpr size_t N = 2;
    constexpr size_t K = 2;
    constexpr size_t S = M*N*K;
    std::array<double, S> data = {1,2,
                                  3,4,
                                      5,6,
                                      7,8,
                                          9,10,
                                          11,12};

    using dense3D = DimLevelType<dense_t, dense_t, dense_t>;
    using layMaj = LoopOrderType<k_t, j_t, i_t>;
    std::array<std::uint32_t, 3> dimensions {M, N, K};

    View<double, Attributes<dense3D, layMaj>> someView (data, dimensions);

    CHECK(someView.rank() == 3);
    CHECK(someView.dims()[0] == M);
    CHECK(someView.dims()[1] == N);
    CHECK(someView.dims()[2] == K);

    for (std::size_t l = 0; l < S; ++l)
    {
        CHECK(data[l] == someView.data()[l]);
    }

    CHECK(data[0] == someView(0, 0, 0));
    CHECK(data[1] == someView(0, 0, 1));
    CHECK(data[3] == someView(0, 1, 1));
    CHECK(data[4] == someView(1, 0, 0));
    CHECK(data[9] == someView(2, 0, 1));

}

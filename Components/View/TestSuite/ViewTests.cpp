#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>


#include "View/View.hpp"

using namespace QuICC::Memory;

TEST_CASE("Level", "[Level]")
{
    CHECK(isLevelType_v<dense_t> == true);
    CHECK(isLevelType_v<compressed_t> == true);
    CHECK(isLevelType_v<triK_t> == true);
    CHECK(isLevelType_v<CSC_t> == true);
    CHECK(isLevelType_v<double> == false);

    CHECK(areLevelType_v<dense_t, dense_t> == true);
    CHECK(areLevelType_v<dense_t, double> == false);

    CHECK(isLevelTypeDense_v<dense_t> == true);
    CHECK(isLevelTypeDense_v<compressed_t> == false);
    CHECK(isLevelTypeDense_v<triK_t> == false);
    CHECK(isLevelTypeDense_v<CSC_t> == false);
    // this should trigger a static assert
    // CHECK(isLevelTypeDense_v<double> == false);

    CHECK(areLevelTypeDense_v<dense_t, dense_t> == true);
    CHECK(areLevelTypeDense_v<dense_t, compressed_t> == false);
    // this should trigger a static assert
    // CHECK(areLevelTypeDense_v<dense_t, double> == false);
}

TEST_CASE("DimLevelType", "[DimLevelType]")
{
    CHECK(isDimLevelType_v<DimLevelType<dense_t>> == true);
    CHECK(isDimLevelType_v<double> == false);
    // this should trigger a static assert
    // CHECK(isDimLevelType_v<DimLevelType<double>> == true);

    CHECK(isDimLevelTypeDense_v<DimLevelType<dense_t>> == true);
    CHECK(isDimLevelTypeDense_v<DimLevelType<dense_t, dense_t>> == true);
    // this should trigger a static assert
    // CHECK(isDimLevelTypeDense_v<dense_t> == false);

    CHECK(isDimLevelTypeDense_v<DimLevelType<compressed_t>> == false);
    CHECK(isDimLevelTypeDense_v<DimLevelType<dense_t, compressed_t>> == false);

    CHECK(howManyLevelTypeDense_v<dense_t> == 1);
    CHECK(howManyLevelTypeDense_v<dense_t, dense_t> == 2);
    CHECK(howManyLevelTypeDense_v<dense_t, compressed_t> == 1);
    // this should trigger a static assert
    // CHECK(howManyLevelTypeDense_v<double> == 1);

    CHECK(howManyDimLevelTypeDense_v<DimLevelType<dense_t>> == 1);
    CHECK(howManyDimLevelTypeDense_v<DimLevelType<dense_t, dense_t>> == 2);
    CHECK(howManyDimLevelTypeDense_v<DimLevelType<dense_t, compressed_t>> == 1);
    // this should trigger a static assert
    // CHECK(howManyDimLevelTypeDense_v<DimLevelType<double>> == 1);
}


TEST_CASE("Attributes", "[Attributes]")
{
    using levelTypeDcbc_t = DimLevelType<dense_t, dense_t, compressed_t>;
    using orderType_t = LoopOrderType<i_t, j_t, k_t>;

    CHECK(hasLevel<Attributes<orderType_t, levelTypeDcbc_t>>::value == true);
    CHECK(hasOrder<Attributes<orderType_t, levelTypeDcbc_t>>::value == true);
}

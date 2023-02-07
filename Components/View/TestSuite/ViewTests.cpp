#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>


#include "View/View.hpp"


TEST_CASE("Attributes", "[Attributes]")
{
    using namespace QuICC::View;

    using levelTypeC_t = DimLevelType<compressed_t>;
    levelTypeC_t c;
    using levelTypeDcbc_t = DimLevelType<dense_t, dense_t, compressed_t>;
    levelTypeDcbc_t ddc;
    using orderType_t = LoopOrderType<i_t, j_t, k_t>;
    orderType_t ijk;

    CHECK(hasLevel<Attributes<orderType_t, levelTypeDcbc_t>>::value == true);
    CHECK(hasOrder<Attributes<orderType_t, levelTypeDcbc_t>>::value == true);
}

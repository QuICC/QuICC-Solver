#include <catch2/catch.hpp>

#include "View/ViewUtils.hpp"

using namespace QuICC::Memory;

TEST_CASE("Index Mirror Map Odd", "[IndexMirrorMapOdd]")
{
    constexpr std::uint32_t M = 5;

    CHECK(mapInOrderIndex(0u, M) == 0);
    CHECK(mapInOrderIndex(1u, M) == 1);
    CHECK(mapInOrderIndex(2u, M) == 2);
    CHECK(mapInOrderIndex(3u, M) == 3);

    CHECK(mapInOrderIndex(M-1, M) == M-1);

}

TEST_CASE("Index Mirror Map Even", "[IndexMirrorMapEven]")
{
    constexpr std::uint32_t M = 20;

    CHECK(mapInOrderIndex(0u, M) == 0);
    CHECK(mapInOrderIndex(1u, M) == 1);
    CHECK(mapInOrderIndex(2u, M) == 2);

    CHECK(mapInOrderIndex(M-1, M) == M-1);

}

TEST_CASE("Index Mirror Map Even Padded", "[IndexMirrorMapEvenPadded]")
{
    constexpr std::uint32_t M = 20;
    constexpr std::uint32_t lds = 30;

    CHECK(mapInOrderIndex(0u, M, lds) == 0);
    CHECK(mapInOrderIndex(1u, M, lds) == 1);
    CHECK(mapInOrderIndex(2u, M, lds) == 2);

    CHECK(mapInOrderIndex(M-1, M, lds) == lds-1);

}


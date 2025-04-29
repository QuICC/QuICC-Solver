#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "Types/BasicTypes.hpp"
#include "Types/Math.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Math.hpp"
#include "Types/Internal/Literals.hpp"

using namespace QuICC;

TEST_CASE("Basic Types", "[BasicTypes]")
{
    CHECK(std::is_same_v<MHDFloat, double>);
    CHECK(sizeof(MHDComplex) == 2*sizeof(MHDFloat));
}

#ifdef QUICC_MULTPRECISION

TEST_CASE("Multiprecision Basic Types", "[MPBasicTypes]")
{
    CHECK(!std::is_same_v<MP::MHDFloat, double>);
    CHECK(sizeof(MP::MHDFloat) > sizeof(double));
    CHECK(sizeof(MP::MHDComplex) == 2*sizeof(MP::MHDFloat));
}

TEST_CASE("Internal Basic Types", "[InternalBasicTypes]")
{
    CHECK(std::is_same_v<Internal::MHDFloat, MP::MHDFloat>);
    CHECK(sizeof(Internal::MHDComplex) == 2*sizeof(Internal::MHDFloat));
}

TEST_CASE("Internal Literals", "[InternalLiterals]")
{
    using namespace Internal::Literals;
    auto one = 1.0_mp;
    CHECK(one == Internal::MHDFloat("1.0"));
    auto oneInt = 1_mp;
    CHECK(one - oneInt == 0.0_mp);
    auto oneDot = 1._mp;
    CHECK(one - oneDot == 0.0_mp);
}

TEST_CASE("PI", "[PI]")
{
    auto iPI = Internal::Math::PI;
    CHECK(sizeof(iPI) == sizeof(Internal::MHDFloat));
    auto eps = Internal::MHDFloat(std::numeric_limits<MHDFloat>::epsilon());
    CHECK(Math::PI - iPI <= eps);
}

#else

TEST_CASE("Internal Basic Types", "[InternalBasicTypes]")
{
    CHECK(std::is_same_v<Internal::MHDFloat, MHDFloat>);
    CHECK(sizeof(Internal::MHDComplex) == 2*sizeof(Internal::MHDFloat));
}

TEST_CASE("Internal Literals", "[InternalLiterals]")
{
    using namespace Internal::Literals;
    auto one = 1.0_mp;
    CHECK(one == Internal::MHDFloat(1.0));
    auto oneInt = 1_mp;
    CHECK(one - oneInt == 0.0_mp);
    auto oneDot = 1._mp;
    CHECK(one - oneDot == 0.0_mp);
}

TEST_CASE("PI", "[PI]")
{
    auto iPI = Internal::Math::PI;
    CHECK(sizeof(iPI) == sizeof(Internal::MHDFloat));
    auto eps = std::numeric_limits<MHDFloat>::epsilon();
    CHECK(Math::PI - iPI <= eps);
}

#endif

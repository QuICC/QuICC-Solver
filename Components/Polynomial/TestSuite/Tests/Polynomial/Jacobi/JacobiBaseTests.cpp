#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>


#include "QuICC/Polynomial/Jacobi/JacobiBase.hpp"

namespace QuICC {

template <class T>
std::pair<bool,T> checkNormal(const T data, const T ref, T refMod = 1.0)
{
    #ifdef QUICC_MULTPRECISION
    constexpr auto maxUlp = 305;
    #else
    constexpr auto maxUlp = 10;
    #endif


    auto epsilon = std::numeric_limits<T>::epsilon();
    auto tol = maxUlp * epsilon;

    bool isEqual = false;
    auto diff = Internal::Math::abs(data-ref);
    if(ref != 0.0)
    {
        refMod = Internal::Math::abs(ref);
    }

    if(diff < tol)
    {
        isEqual = diff < (tol * refMod);
    }
    else
    {
        isEqual = (diff / refMod ) < tol;
    }

    auto ulp = diff / (refMod * epsilon);

    return std::make_pair(isEqual, ulp);
}

#ifdef QUICC_MULTPRECISION
constexpr uint digits = 32;
#else
constexpr uint digits = 16;
#endif


template<int N>
struct refHolder
{
    using type = QuICC::Internal::MHDFloat;
    type a, b;
    std::array<type, N> res;
};

using namespace QuICC::Polynomial::Jacobi;
using namespace QuICC::Polynomial::Quadrature;

TEST_CASE("Natural P0ab", "[NaturalP0ab]")
{
    std::vector<refHolder<1>> ref{
        {MHD_MP(-0.5), MHD_MP(-0.5), {MHD_MP(1.0)}},
        {MHD_MP(-0.5), MHD_MP(0.5), {MHD_MP(1.0)}},
        {MHD_MP(0.5), MHD_MP(-0.5), {MHD_MP(1.0)}}};


    for(std::size_t i = 0; i < ref.size(); ++i)
    {
        auto cs = JacobiBase::P0ab<natural_t>(ref.at(i).a, ref.at(i).b);
        for(std::size_t n = 0; n < ref.at(i).res.size(); ++n)
        {
            auto err = checkNormal(cs(n), ref.at(i).res.at(n));
            {
              INFO( "checked normal value" );
              INFO( "a: " << std::scientific << std::setprecision(digits) << ref.at(i).a );
              INFO( "b: " << std::scientific << std::setprecision(digits) << ref.at(i).b );
              INFO( "ref: " << std::scientific << std::setprecision(digits) << ref.at(i).res.at(n) );
              INFO( "check: " << std::scientific << std::setprecision(digits) << cs(n) );
              INFO( "measured ulp: " << err.second);
              CHECK(err.first);
            }
        }
   }
}


TEST_CASE("Natural P1ab", "[NaturalP1ab]")
{
    std::vector<refHolder<3>> ref{
        {MHD_MP(-0.5), MHD_MP(-0.5), {MHD_MP(1.0), MHD_MP(0.0), MHD_MP(0.5)}},
        {MHD_MP(-0.5), MHD_MP(0.5), {MHD_MP(2.0), MHD_MP(-1.0), MHD_MP(0.5)}},
        {MHD_MP(0.5), MHD_MP(-0.5), {MHD_MP(2.0), MHD_MP(1.0), MHD_MP(0.5)}}};


    for(std::size_t i = 0; i < ref.size(); ++i)
    {
        auto cs = JacobiBase::P1ab<natural_t>(ref.at(i).a, ref.at(i).b);
        for(std::size_t n = 0; n < ref.at(i).res.size(); ++n)
        {
            auto err = checkNormal(cs(n), ref.at(i).res.at(n));
            {
              INFO( "checked normal value" );
              INFO( "a: " << std::scientific << std::setprecision(digits) << ref.at(i).a );
              INFO( "b: " << std::scientific << std::setprecision(digits) << ref.at(i).b );
              INFO( "ref: " << std::scientific << std::setprecision(digits) << ref.at(i).res.at(n) );
              INFO( "check: " << std::scientific << std::setprecision(digits) << cs(n) );
              INFO( "measured ulp: " << err.second);
              CHECK(err.first);
            }
        }
   }
}


TEST_CASE("Natural P2ab", "[NaturalP2ab]")
{
    std::vector<refHolder<4>> ref{
        {MHD_MP(-0.5), MHD_MP(-0.5),
            {MHD_MP(-0.375000000000000000000000000000000000000000000000000000000000000),
             MHD_MP(1.500000000000000000000000000000000000000000000000000000000000000),
             MHD_MP(0.0),
             MHD_MP(1.0)}},
        {MHD_MP(-0.5), MHD_MP(0.0),
            {MHD_MP(-0.388888888888888888888888888888888888888888888888888888888888889),
             MHD_MP(1.458333333333333333333333333333333333333333333333333333333333333),
             MHD_MP(0.0694444444444444444444444444444444444444444444444444444444444444),
             MHD_MP(1.0)}},
        {MHD_MP(-0.5), MHD_MP(0.5),
            {MHD_MP(-0.375000000000000000000000000000000000000000000000000000000000000),
             MHD_MP(1.500000000000000000000000000000000000000000000000000000000000000),
             MHD_MP(0.0),
             MHD_MP(1.0)}}
    };

    int n = 2;

    for(std::size_t i = 0; i < ref.size(); ++i)
    {
        auto cs = JacobiBase::Pnab<natural_t>(n, ref.at(i).a, ref.at(i).b);
        for(std::size_t c = 0; c < ref.at(i).res.size(); ++c)
        {
            auto err = checkNormal(cs(c), ref.at(i).res.at(c));
            {
              INFO( "checked normal value" );
              INFO( "a: " << std::scientific << std::setprecision(digits) << ref.at(i).a );
              INFO( "b: " << std::scientific << std::setprecision(digits) << ref.at(i).b );
              INFO( "ref: " << std::scientific << std::setprecision(digits) << ref.at(i).res.at(c) );
              INFO( "check: " << std::scientific << std::setprecision(digits) << cs(c) );
              INFO( "measured ulp: " << err.second);
              CHECK(err.first);
            }
        }
   }
}

TEST_CASE("Natural dP1ab", "[NaturaldP1ab]")
{
    std::vector<refHolder<1>> ref{
        {MHD_MP(-0.5), MHD_MP(-0.5), {MHD_MP(0.5)}},
        {MHD_MP(-0.5), MHD_MP(0.5), {MHD_MP(1.0)}},
        {MHD_MP(0.5), MHD_MP(-0.5), {MHD_MP(1.0)}}};


    for(std::size_t i = 0; i < ref.size(); ++i)
    {
        auto cs = JacobiBase::dP1ab<natural_t>(ref.at(i).a, ref.at(i).b);
        for(std::size_t n = 0; n < ref.at(i).res.size(); ++n)
        {
            auto err = checkNormal(cs(n), ref.at(i).res.at(n));
            {
              INFO( "checked normal value" );
              INFO( "a: " << std::scientific << std::setprecision(digits) << ref.at(i).a );
              INFO( "b: " << std::scientific << std::setprecision(digits) << ref.at(i).b );
              INFO( "ref: " << std::scientific << std::setprecision(digits) << ref.at(i).res.at(n) );
              INFO( "check: " << std::scientific << std::setprecision(digits) << cs(n) );
              INFO( "measured ulp: " << err.second);
              CHECK(err.first);
            }
        }
   }
}

TEST_CASE("Natural dP2ab", "[NaturaldP2ab]")
{
    std::vector<refHolder<3>> ref{
        {MHD_MP(-0.5), MHD_MP(-0.5), {MHD_MP(6.0), MHD_MP(0.0), MHD_MP(0.25)}},
        {MHD_MP(-0.5), MHD_MP(0.5), {MHD_MP(12.0), MHD_MP(-3.0), MHD_MP(0.25)}}
    };


    for(std::size_t i = 0; i < ref.size(); ++i)
    {
        auto cs = JacobiBase::dP2ab<natural_t>(ref.at(i).a, ref.at(i).b);
        for(std::size_t n = 0; n < ref.at(i).res.size(); ++n)
        {
            auto err = checkNormal(cs(n), ref.at(i).res.at(n));
            {
              INFO( "checked normal value" );
              INFO( "n: " << n );
              INFO( "a: " << std::scientific << std::setprecision(digits) << ref.at(i).a );
              INFO( "b: " << std::scientific << std::setprecision(digits) << ref.at(i).b );
              INFO( "ref: " << std::scientific << std::setprecision(digits) << ref.at(i).res.at(n) );
              INFO( "check: " << std::scientific << std::setprecision(digits) << cs(n) );
              INFO( "measured ulp: " << err.second);
              CHECK(err.first);
            }
        }
   }
}


TEST_CASE("Unity P0ab", "[UnityP0ab]")
{
    std::vector<refHolder<1>> ref{
        {MHD_MP(-0.5), MHD_MP(-0.5), {MHD_MP(0.5641895835477562869480794515607725858440506293289988568440857217)}},
        {MHD_MP(-0.5), MHD_MP(0.5), {MHD_MP(0.5641895835477562869480794515607725858440506293289988568440857217)}},
        {MHD_MP(0.5), MHD_MP(-0.5), {MHD_MP(0.5641895835477562869480794515607725858440506293289988568440857217)}},
        {MHD_MP(0.5), MHD_MP(0.0), {MHD_MP(0.72823765756098513042558094123662621856988478676737281554594382465)}},
        {MHD_MP(0.5), MHD_MP(1.0), {MHD_MP(0.814194453040788226251971229084579614171877058010907414680716)}},
        {MHD_MP(0.25), MHD_MP(0.0), {MHD_MP(0.724955350027552793005592872887516632480952114932488532219682)}}};

    for(std::size_t i = 0; i < ref.size(); ++i)
    {
        auto cs = JacobiBase::P0ab<unity_t>(ref.at(i).a, ref.at(i).b);
        for(std::size_t n = 0; n < ref.at(i).res.size(); ++n)
        {
            auto err = checkNormal(cs(n), ref.at(i).res.at(n));
            {
              INFO( "checked normal value" );
              INFO( "a: " << std::scientific << std::setprecision(digits) << ref.at(i).a );
              INFO( "b: " << std::scientific << std::setprecision(digits) << ref.at(i).b );
              INFO( "ref: " << std::scientific << std::setprecision(digits) << ref.at(i).res.at(n) );
              INFO( "check: " << std::scientific << std::setprecision(digits) << cs(n) );
              INFO( "measured ulp: " << err.second);
              CHECK(err.first);
            }
        }
   }
}


TEST_CASE("Unity P1ab", "[UnityP1ab]")
{
    std::vector<refHolder<3>> ref{
        {MHD_MP(-0.5), MHD_MP(-0.5), {MHD_MP(1.0), MHD_MP(0.0), MHD_MP(0.797884560802865355879892119868763736951717262329869315331851659)}},
        {MHD_MP(-0.5), MHD_MP(0.0), {MHD_MP(1.5), MHD_MP(-0.5), MHD_MP(0.664786987118123577031280357225814518825217305962988633471961216)}},
        {MHD_MP(-0.5), MHD_MP(0.5), {MHD_MP(2.0), MHD_MP(-1.0), MHD_MP(0.564189583547756286948079451560772585844050629328998856844085722)}}};


    for(std::size_t i = 0; i < ref.size(); ++i)
    {
        auto cs = JacobiBase::P1ab<unity_t>(ref.at(i).a, ref.at(i).b);
        for(std::size_t n = 0; n < ref.at(i).res.size(); ++n)
        {
            auto err = checkNormal(cs(n), ref.at(i).res.at(n));
            {
              INFO( "checked normal value" );
              INFO( "n: " << n );
              INFO( "a: " << std::scientific << std::setprecision(digits) << ref.at(i).a );
              INFO( "b: " << std::scientific << std::setprecision(digits) << ref.at(i).b );
              INFO( "ref: " << std::scientific << std::setprecision(digits) << ref.at(i).res.at(n) );
              INFO( "check: " << std::scientific << std::setprecision(digits) << cs(n) );
              INFO( "measured ulp: " << err.second);
              CHECK(err.first);
            }
        }
   }
}

TEST_CASE("Normalization Factor", "[NormFact]")
{
    struct refHolderN
    {
        using type = QuICC::Internal::MHDFloat;
        type n, a, b, res;
    };

    std::vector<refHolderN> ref{
        {MHD_MP(1), MHD_MP(-0.5), MHD_MP(-0.5), MHD_MP(0.785398163397448309615660845819875721049292349843776455243736148)},
        {MHD_MP(2), MHD_MP(-0.5), MHD_MP(-0.5), MHD_MP(0.883572933822129348317618451547360186180453893574248512149203167)},
        {MHD_MP(3), MHD_MP(-0.5), MHD_MP(-0.5), MHD_MP(0.92038847273138473783085255369516686060463947247317553348875330)},
        {MHD_MP(16), MHD_MP(-0.5), MHD_MP(-0.5), MHD_MP(0.9844989377871177684244116841956096822602866127568565908389632)},
        {MHD_MP(32), MHD_MP(-0.5), MHD_MP(-0.5), MHD_MP(0.9922182535855425116918532960127701416725371605079273446118359)},
        {MHD_MP(64), MHD_MP(-0.5), MHD_MP(-0.5), MHD_MP(0.996101409048731968484823368379419314990064678170488277338282)},
        {MHD_MP(128), MHD_MP(-0.5), MHD_MP(-0.5), MHD_MP(0.998048786064746639419685302726397264944637331237308355769909)},
        {MHD_MP(256), MHD_MP(-0.5), MHD_MP(-0.5), MHD_MP(0.999023914802248505635669383232439670509445515942105705918336)},
        {MHD_MP(512), MHD_MP(-0.5), MHD_MP(-0.5), MHD_MP(0.99951183801746160480578738080065059505696440146636599643552)},
        {MHD_MP(1024), MHD_MP(-0.5), MHD_MP(-0.5), MHD_MP(0.99975588918459612237048504360740628127637299927335371879644)},
        {MHD_MP(1), MHD_MP(-0.5), MHD_MP(0.5), MHD_MP(2.356194490192344928846982537459627163147877049531329365731208444)},
        {MHD_MP(2), MHD_MP(-0.5), MHD_MP(0.5), MHD_MP(2.20893233455532337079404612886840046545113473393562128037300792)},
        {MHD_MP(3), MHD_MP(-0.5), MHD_MP(0.5), MHD_MP(2.14757310303989772160532262528872267474415876910407624480709103)},
        {MHD_MP(64), MHD_MP(-0.5), MHD_MP(0.5), MHD_MP(2.007766902613850373977222101889767056776849116937390434009974)},
        {MHD_MP(128), MHD_MP(-0.5), MHD_MP(0.5), MHD_MP(2.003894828270624111959836896880344508521654641624908183069270)},
        {MHD_MP(192), MHD_MP(-0.5), MHD_MP(0.5), MHD_MP(2.002599087028509079893743254076522378451128147127350697049606)},
        {MHD_MP(256), MHD_MP(-0.5), MHD_MP(0.5), MHD_MP(2.00195026677169329449647809999313105848181855343086026225042)},
        {MHD_MP(512), MHD_MP(-0.5), MHD_MP(0.5), MHD_MP(2.00097584759355106430846106507942746080739943652934598895782)},
        {MHD_MP(1024), MHD_MP(-0.5), MHD_MP(0.5), MHD_MP(2.00048810247972407689172251401521042024930495655381032208389)},
        {MHD_MP(2048), MHD_MP(-0.5), MHD_MP(0.5), MHD_MP(2.0002440959269735529774897203976916698552378805558691998625)},
        {MHD_MP(16), MHD_MP(0.5), MHD_MP(0.0), MHD_MP(2.8284271247461900976033774484193961571393437507538961463533595)},
        {MHD_MP(32), MHD_MP(0.5), MHD_MP(0.0), MHD_MP(2.828427124746190097603377448419396157139343750753896146353359)},
        {MHD_MP(64), MHD_MP(0.5), MHD_MP(0.0), MHD_MP(2.828427124746190097603377448419396157139343750753896146353359)},
        {MHD_MP(128), MHD_MP(0.5), MHD_MP(0.0), MHD_MP(2.828427124746190097603377448419396157139343750753896146353359)},
        {MHD_MP(256), MHD_MP(0.5), MHD_MP(0.0), MHD_MP(2.828427124746190097603377448419396157139343750753896146353359)},
        {MHD_MP(512), MHD_MP(0.5), MHD_MP(0.0), MHD_MP(2.82842712474619009760337744841939615713934375075389614635336)},
        {MHD_MP(1024), MHD_MP(0.5), MHD_MP(0.0), MHD_MP(2.82842712474619009760337744841939615713934375075389614635336)},
        {MHD_MP(2048), MHD_MP(0.5), MHD_MP(0.0), MHD_MP(2.8284271247461900976033774484193961571393437507538961463534)},
        {MHD_MP(16), MHD_MP(0.5), MHD_MP(1.0), MHD_MP(5.495229842364026475343704756929112533870725001464712512915098)},
        {MHD_MP(32), MHD_MP(0.5), MHD_MP(1.0), MHD_MP(5.572423589052195416173818256587467055856617538798720765949902)},
        {MHD_MP(64), MHD_MP(0.5), MHD_MP(1.0), MHD_MP(5.613672155984804773869298752588114510352895993862694641617355)},
        {MHD_MP(128), MHD_MP(0.5), MHD_MP(1.0), MHD_MP(5.63501311339395401684688325631045720881815202852899772787001)},
        {MHD_MP(256), MHD_MP(0.5), MHD_MP(1.0), MHD_MP(5.64587006648365712686654760577696941658105898208738881252671)},
        {MHD_MP(512), MHD_MP(0.5), MHD_MP(1.0), MHD_MP(5.65134611487748985421823809557604762848094778631645072280145)},
        {MHD_MP(1024), MHD_MP(0.5), MHD_MP(1.0), MHD_MP(5.6540961538075960020350304917208796900396437728381151633587)},
        {MHD_MP(2048), MHD_MP(0.5), MHD_MP(1.0), MHD_MP(5.6554741923444191363643038710039938775101393952619987351823)},
    };

    for(std::size_t i = 0; i < ref.size(); ++i)
    {
        auto cs = QuICC::Polynomial::Jacobi::JacobiBase::normFact(ref.at(i).n, ref.at(i).a, ref.at(i).b);
        auto err = checkNormal(cs, ref.at(i).res);
        {
          INFO( "checked normal value" );
          INFO( "n: " << std::scientific << std::setprecision(digits) << ref.at(i).n );
          INFO( "a: " << std::scientific << std::setprecision(digits) << ref.at(i).a );
          INFO( "b: " << std::scientific << std::setprecision(digits) << ref.at(i).b );
          INFO( "ref: " << std::scientific << std::setprecision(digits) << ref.at(i).res );
          INFO( "check: " << std::scientific << std::setprecision(digits) << cs );
          INFO( "measured ulp: " << err.second);
          CHECK(err.first);
        }
   }
}


} // namespace QuICC

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>


#include "QuICC/Polynomial/Jacobi/JacobiAsymptotic.hpp"

namespace QuICC {

template <class T>
std::pair<bool,T> checkNormal(const T data, const T ref, T refMod = 1.0)
{
    constexpr auto maxUlp = 12;

    auto epsilon = std::numeric_limits<T>::epsilon();
    auto tol = maxUlp * epsilon;

    bool isEqual = false;
    auto diff = precision::abs(data-ref);
    if(ref != 0.0)
    {
        refMod = precision::abs(ref);
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

constexpr uint digits = 16;

TEST_CASE("Stirling's serie", "[Stirling]")
{
    struct refHolderN
    {
        using type = QuICC::internal::MHDFloat;
        type coeff;
    };

    std::vector<refHolderN> ref{
        {MHD_MP(1.000000000000000000000000000000000000000000000000000000000000000)},
        {MHD_MP(0.0833333333333333333333333333333333333333333333333333333333333333)},
        {MHD_MP(0.003472222222222222222222222222222222222222222222222222222222222222)},
        {MHD_MP(-0.002681327160493827160493827160493827160493827160493827160493827160)},
        {MHD_MP(-0.0002294720936213991769547325102880658436213991769547325102880658436)},
        {MHD_MP(0.000784039221720066627474034881442288849696257103664511071918479326)},
        {MHD_MP(0.0000697281375836585777429398828575783308293596359439980839157793890)},
        {MHD_MP(-0.000592166437353693882864836225604401187391585196797816825251667501)},
        {MHD_MP(-0.0000517179090826059219337057843002058822817853453427298877687537953)},
        {MHD_MP(0.000839498720672087279993357516764983445198182111593007628503572550)},
        // {MHD_MP(0.0000720489541602001055908571930225015052063451737975478952504489810)},
    };

    for(std::size_t i = 0; i < ref.size(); ++i)
    {
        auto val = QuICC::Polynomial::Jacobi::JacobiAsy::details::stirling();
        auto err = checkNormal(val.at(i), ref.at(i).coeff);
        {
          INFO( "checked normal value" );
          INFO( "ref: " << std::scientific << std::setprecision(digits) << ref.at(i).coeff );
          INFO( "check: " << std::scientific << std::setprecision(digits) << val.at(i) );
          INFO( "measured ulp: " << err.second);
          CHECK(err.first);
        }
   }
}

TEST_CASE("g(theta), g', g''", "[gTheta]")
{
    struct refHolderN
    {
        using type = QuICC::internal::MHDFloat;
        type t, a, b;
        std::array<type, 3> res;
    };

    std::vector<refHolderN> ref{
        {MHD_MP(0.1), MHD_MP(-0.5), MHD_MP(-0.5), {MHD_MP(0.0),MHD_MP(0.0),MHD_MP(0.0)}},
        {MHD_MP(0.1), MHD_MP(-0.5), MHD_MP(0.5), {MHD_MP(0.0),MHD_MP(0.0),MHD_MP(0.0)}},
        {Precision::PI/MHD_MP(2.0), MHD_MP(0.0), MHD_MP(0.0),
            {MHD_MP(-0.3183098861837906715377675267450287240689192914809128974953346881),
             MHD_MP(-0.2973576327153244571122410735805447221912824506555069023087819362),
             MHD_MP(-0.2580122754655959134753764215085095087445066891798246138766537777)}},
        {MHD_MP(0.3926990816987241548078304229099378605246461749218882276218680740),
            MHD_MP(0.0), MHD_MP(0.5),
            {MHD_MP(-0.016404671703700660022326339214096799841378602904779807912719635),
             MHD_MP(-0.04198991958121385825320631550394645795965236821486172341839776),
             MHD_MP(-0.0016564355663884440080409393332061391065818187009177776942674)}},
        {MHD_MP(0.04908738521234051935097880286374223256558077186523602845273350925),
            MHD_MP(0.0), MHD_MP(0.25),
            {MHD_MP(-0.00664825650589955979220529716323336282359967184302968099768144),
             MHD_MP(-0.135478184104041060780607956607584280704385308240952445406824),
             MHD_MP(-0.0025073901369084102346611694170402907154377338280083176805)}},
    };

    for(std::size_t i = 0; i < ref.size(); ++i)
    {
        auto val = QuICC::Polynomial::Jacobi::JacobiAsy::details::g(ref.at(i).a, ref.at(i).b, ref.at(i).t);
        for(std::size_t p = 0; p < ref.at(0).res.size(); ++p)
        {
            auto err = checkNormal(val.at(p), ref.at(i).res.at(p));
            {
              INFO( "checked normal value" );
              INFO( "t: " << std::scientific << std::setprecision(digits) << ref.at(i).t );
              INFO( "a: " << std::scientific << std::setprecision(digits) << ref.at(i).a );
              INFO( "b: " << std::scientific << std::setprecision(digits) << ref.at(i).b );
              INFO( "p: " << std::scientific << std::setprecision(digits) << p );
              INFO( "ref: " << std::scientific << std::setprecision(digits) << ref.at(i).res.at(p) );
              INFO( "check: " << std::scientific << std::setprecision(digits) << val.at(p) );
              INFO( "measured ulp: " << err.second);
              CHECK(err.first);
            }
        }
   }
}

TEST_CASE("ExpArgBnd", "[ExpArgBnd]")
{
    struct refHolderN
    {
        using type = QuICC::internal::MHDFloat;
        type n, a, res;
    };

    std::vector<refHolderN> ref{
        {10, MHD_MP(0.5), MHD_MP(0.01229672377903603218643124434322891538372345394565905064302464935)},
        {40, MHD_MP(0.25), MHD_MP(0.000779627463101987563721798136291069098611488688243252164597936297)},
        {60, MHD_MP(-0.125), MHD_MP(0.0001302988500949943164451558562232065302770877073177752516648790404)},

    };

    for(std::size_t i = 0; i < ref.size(); ++i)
    {
        auto val = QuICC::Polynomial::Jacobi::JacobiAsy::details::expArgBnd(ref.at(i).n, ref.at(i).a);
        auto err = checkNormal(val, ref.at(i).res);
        {
            INFO( "checked normal value" );
            INFO( "n: " << std::scientific << std::setprecision(digits) << ref.at(i).n );
            INFO( "a: " << std::scientific << std::setprecision(digits) << ref.at(i).a );
            INFO( "ref: " << std::scientific << std::setprecision(digits) << ref.at(i).res );
            INFO( "check: " << std::scientific << std::setprecision(digits) << val );
            INFO( "measured ulp: " << err.second);
            CHECK(err.first);
        }
   }
}

TEST_CASE("ExpArgInt", "[ExpArgInt]")
{
    struct refHolderN
    {
        using type = QuICC::internal::MHDFloat;
        type n, a, b, res;
    };

    std::vector<refHolderN> ref{
        {10, MHD_MP(0.25), MHD_MP(0.0), MHD_MP(0.00154324908052553084296781844409464546136742198813989277785319731)},
        {40, MHD_MP(0.5), MHD_MP(-0.5), MHD_MP(0.00625016277059004428605174525338263254019739521539231317180683707)},
        {60, MHD_MP(1.0), MHD_MP(0.5), MHD_MP(0.00102880948820345136328897842600290411137629781069498882438928415)},

    };





    for(std::size_t i = 0; i < ref.size(); ++i)
    {
        auto val = QuICC::Polynomial::Jacobi::JacobiAsy::details::expArgInt(ref.at(i).n, ref.at(i).a, ref.at(i).b);
        auto err = checkNormal(val, ref.at(i).res);
        {
            INFO( "checked normal value" );
            INFO( "n: " << std::scientific << std::setprecision(digits) << ref.at(i).n );
            INFO( "a: " << std::scientific << std::setprecision(digits) << ref.at(i).a );
            INFO( "b: " << std::scientific << std::setprecision(digits) << ref.at(i).b );
            INFO( "ref: " << std::scientific << std::setprecision(digits) << ref.at(i).res );
            INFO( "check: " << std::scientific << std::setprecision(digits) << val );
            INFO( "measured ulp: " << err.second);
            CHECK(err.first);
        }
   }
}


TEST_CASE("f(theta)", "[fTheta]")
{
    struct refHolderN
    {
        using type = double;
        type n, a, b, t, m, res;
    };

    std::vector<refHolderN> ref{
        { MHD_MP(10), MHD_MP(0.25), MHD_MP(0.0),
            MHD_MP(0.785398163397448309615660845819875721049292349843776455243736148),
            MHD_MP(2.0),
            MHD_MP(0.23690635892773916502897124025395384096053131295104928773679985)},
        { MHD_MP(20), MHD_MP(0.25), MHD_MP(0.0),
            MHD_MP(1.570796326794896619231321691639751442098584699687552910487472296),
            MHD_MP(5.0),
            MHD_MP(-52.169142768872091929830818974362263422912854916840817515620395)},
        { MHD_MP(10), MHD_MP(0.5), MHD_MP(0.5),
            MHD_MP(0.785398163397448309615660845819875721049292349843776455243736148),
            MHD_MP(2.0),
            MHD_MP(0.0)},
        { MHD_MP(10), MHD_MP(-0.5), MHD_MP(0.5),
            MHD_MP(0.785398163397448309615660845819875721049292349843776455243736148),
            MHD_MP(2.0),
            MHD_MP(0.0)},

    };





    for(std::size_t i = 0; i < ref.size(); ++i)
    {
        auto val = QuICC::Polynomial::Jacobi::JacobiAsy::details::
            f(ref.at(i).n, ref.at(i).a, ref.at(i).b, ref.at(i).t, ref.at(i).m);
        auto err = checkNormal(val, ref.at(i).res);
        {
            INFO( "checked normal value" );
            INFO( "n: " << std::scientific << std::setprecision(digits) << ref.at(i).n );
            INFO( "a: " << std::scientific << std::setprecision(digits) << ref.at(i).a );
            INFO( "b: " << std::scientific << std::setprecision(digits) << ref.at(i).b );
            INFO( "t: " << std::scientific << std::setprecision(digits) << ref.at(i).t );
            INFO( "m: " << std::scientific << std::setprecision(digits) << ref.at(i).m );
            INFO( "ref: " << std::scientific << std::setprecision(digits) << ref.at(i).res );
            INFO( "check: " << std::scientific << std::setprecision(digits) << val );
            INFO( "measured ulp: " << err.second);
            CHECK(err.first);
        }
   }
}

} // namespace QuICC

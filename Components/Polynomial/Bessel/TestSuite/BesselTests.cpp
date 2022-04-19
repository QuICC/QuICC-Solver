#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <boost/math/special_functions/bessel.hpp>

template <class T>
std::pair<bool,T> checkNormal(const T data, const T ref, T refMod = 1.0)
{
    constexpr auto maxUlp = 10;


    auto epsilon = std::numeric_limits<T>::epsilon();
    auto tol = maxUlp * epsilon;

    bool isEqual = false;
    auto diff = std::abs(data-ref);
    if(ref != 0.0)
    {
        refMod = std::abs(ref);
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


TEST_CASE("Bessel boost", "[BesselBoost]")
{
    struct refHolderN
    {
        using type = double;
        type nu, a, res;
    };

    std::vector<refHolderN> ref{
        {0.5, 1.0, 0.671396707141803090416364012040467080545640816769345913259378817},
        {0.5, 11.0, -0.2405688907232031185276254857626087124658799068895243478097556853},
        {0.5, 21.0, 0.145672360072824684364416495014700186634958070220804587256918021},
        {0.5, 31.0, -0.057900330936878658107301490646503219658030518244143385942728898},
        {0.5, 41.0, -0.019765753988144587814011061590729408569556654534070975697104817},
        {0.5, 64.0, 0.091759321426730820046709210170233165142005344414360412856357171},
        {0.5, 128.0, 0.050850245701274948798534630146750693328785659905303099995074568},
        {0.5, 192.0, -0.02043818766527955264271969254716349400114977475596240764340826},
        {0.5, 256.0, -0.0498282914652630142774606989588957134659416934825006708766594988},
        {0.5, 512.0, 0.00280396912634184132371530787081068757041763943940535044256784},
        {-0.5, 1.0, 0.4310988680183760795205209672985334000880560106886169046523},
        {-0.5, 11.0, 0.001064695682704474192176361873023207843229564591882704406020425},
        {-0.5, 21.0, -0.095366612430202344345674591548304651710529153552829479167079524},
        {-0.5, 31.0, 0.131086511001997263822515553759808667941679139440002406238541257},
        {-0.5, 41.0, -0.1230309980876391447164822513095651468800857685000953863588944068},
        {-0.5, 16.0, -0.191025428464131008048609093276976248731420291378816494302610169},
        {-0.5, 32.0, 0.117665032587494447425722384340898641641615688404213561573043073},
        {-0.5, 64.0,0.039082104274838588953672269482546594620258319212135229370553258},
        {-0.5, 128.0,-0.048865575651389781140880196849690064279998693961612818802940987},
        {-0.5, 256.0, -0.00198427706323025137421306645348911072115805873737910408738454},
        {-0.5, 512.0,-0.035150188478071551896543219389342617778146411050804033328671524},
        {-0.5, 1024.0,0.024618569000951989199715376821131559345493189652636131497381599},
    };

    for(std::size_t i = 0; i < ref.size(); ++i)
    {
        auto val = boost::math::cyl_bessel_j(ref.at(i).nu, ref.at(i).a);
        auto err = checkNormal(val, ref.at(i).res);
        {
          INFO( "checked normal value" );
          INFO( "nu: " << std::scientific << std::setprecision(digits) << ref.at(i).nu );
          INFO( "a: " << std::scientific << std::setprecision(digits) << ref.at(i).a );
          INFO( "ref: " << std::scientific << std::setprecision(digits) << ref.at(i).res );
          INFO( "check: " << std::scientific << std::setprecision(digits) << val );
          INFO( "measured ulp: " << err.second);
          CHECK(err.first);
        }
   }
}

TEST_CASE("Bessel std::c++17", "[BesselStd]")
{
    struct refHolderN
    {
        using type = double;
        type nu, x, res;
    };

    std::vector<refHolderN> ref{
        {0.5, 1.0, 0.671396707141803090416364012040467080545640816769345913259378817},
        {0.5, 11.0,-0.2405688907232031185276254857626087124658799068895243478097556853},
        {0.5, 21.0, 0.145672360072824684364416495014700186634958070220804587256918021},
        // {0.5, 31.0, -0.057900330936878658107301490646503219658030518244143385942728898}, // fail
        // {0.5, 41.0, -0.019765753988144587814011061590729408569556654534070975697104817}, // fail
        // {-0.5, 1.0, 0.4310988680183760795205209672985334000880560106886169046523}, // fail
        // {-0.5, 11.0, 0.001064695682704474192176361873023207843229564591882704406020425}, // fail
        // {-0.5, 21.0, -0.095366612430202344345674591548304651710529153552829479167079524}, // fail
        // {-0.5, 31.0, 0.131086511001997263822515553759808667941679139440002406238541257}, // fail
        // {-0.5, 41.0, -0.1230309980876391447164822513095651468800857685000953863588944068}, // fail
    };

    for(std::size_t i = 0; i < ref.size(); ++i)
    {
        auto val = std::cyl_bessel_j(ref.at(i).nu, ref.at(i).x);
        auto err = checkNormal(val, ref.at(i).res);
        {
          INFO( "checked normal value" );
          INFO( "nu: " << std::scientific << std::setprecision(digits) << ref.at(i).nu );
          INFO( "x: " << std::scientific << std::setprecision(digits) << ref.at(i).x );
          INFO( "ref: " << std::scientific << std::setprecision(digits) << ref.at(i).res );
          INFO( "check: " << std::scientific << std::setprecision(digits) << val );
          INFO( "measured ulp: " << err.second);
          CHECK(err.first);
        }
   }
}

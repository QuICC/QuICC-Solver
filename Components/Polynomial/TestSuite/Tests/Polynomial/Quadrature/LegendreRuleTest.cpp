/**
 * @file LegendreRuleTest.cpp
 * @brief Tests for the Legendre quadrature rule
 */


// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "Types/Precision.hpp"
#include "QuICC/TestSuite/Polynomial/Quadrature/Tester.hpp"
#include "QuICC/TestSuite/Polynomial/Quadrature/TestArgs.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"

namespace currentts = ::QuICC::TestSuite::Polynomial::Quadrature;

TEST_CASE( "Error for Legendre quadrature rule", "[LegendreRule]" ){
   // Set default arguments if required
   if(currentts::args().useDefault)
   {
      currentts::args().gridN.clear();
      currentts::args().gridN.push_back(47);
      currentts::args().gridN.push_back(48);
      currentts::args().gridN.push_back(256);
      currentts::args().gridN.push_back(255);
      currentts::args().gridN.push_back(511);
      currentts::args().gridN.push_back(512);
      currentts::args().gridN.push_back(1024);
      if(QuICC::Precision::hasMP)
      {
         currentts::args().gridN.push_back(1022);
      }
      //currentts::args().gridN.push_back(1025);
      //currentts::args().gridN.push_back(4098);
      //currentts::args().gridN.push_back(4093);
   }

   Catch::StringMaker<double>::precision = 15;

   REQUIRE( currentts::args().gridN.size() > 0 );
   REQUIRE( currentts::args().params.size() == 0.0 );

   typedef ::QuICC::Polynomial::Quadrature::LegendreRule QuadratureType;
   typedef typename currentts::Tester<QuadratureType> Tester;

   Tester tester("legendre.dat", currentts::args().dumpData);
   tester.setUlp(currentts::args().ulp);

   Tester::ParameterType params;

   for(int gN: currentts::args().gridN)
   {
      tester.validate(currentts::args().specN, gN, params, currentts::args().type);
   }
}

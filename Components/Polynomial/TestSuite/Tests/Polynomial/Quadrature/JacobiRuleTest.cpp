/**
 * @file JacobiRuleTest.cpp
 * @brief Tests for the Jacobi quadrature rule
 */


// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "QuICC/TestSuite/Polynomial/Quadrature/Tester.hpp"
#include "QuICC/TestSuite/Polynomial/Quadrature/TestArgs.hpp"
#include "QuICC/Polynomial/Quadrature/JacobiRule.hpp"

namespace currentts = ::QuICC::TestSuite::Polynomial::Quadrature;

TEST_CASE( "Error for Jacobi quadrature rule", "[JacobiRule]" ){
   // Set default arguments if required
   if(currentts::args().useDefault)
   {
      currentts::args().gridN.clear();
      currentts::args().gridN.push_back(2);
      currentts::args().gridN.push_back(3);
      currentts::args().gridN.push_back(4);
      currentts::args().gridN.push_back(12);
      currentts::args().gridN.push_back(24);
      currentts::args().gridN.push_back(47);
      currentts::args().gridN.push_back(48);
      currentts::args().gridN.push_back(64);
      currentts::args().gridN.push_back(128);
      currentts::args().gridN.push_back(256);
      currentts::args().gridN.push_back(255);
      // // alpha/beta
      currentts::args().ab.push_back({-1./2., -1./2.});
      currentts::args().ab.push_back({-1./2., 1./2.});
      currentts::args().ab.push_back({1./2., 0.0});
      currentts::args().ab.push_back({1./2., 1.0});
      currentts::args().ab.push_back({1./4., 1./4.});
      currentts::args().ab.push_back({1./4., 0.0});
   }

   Catch::StringMaker<double>::precision = 15;

   REQUIRE( currentts::args().gridN.size() > 0 );
   REQUIRE( currentts::args().ab.size() > 0 );

   typedef ::QuICC::Polynomial::Quadrature::JacobiRule QuadratureType;
   typedef typename currentts::Tester<QuadratureType,2> Tester;

   Tester tester("jacobi.dat", currentts::args().dumpData);
   tester.setUlp(currentts::args().ulp);

   Tester::ParameterType params = {-1, -1};

   for(auto&& gN: currentts::args().gridN)
   {
      for(auto&& ab: currentts::args().ab)
      {
         params.at(0) = ab.at(0);
         params.at(1) = ab.at(1);
         tester.validate(currentts::args().specN, gN, params, currentts::args().type);
      }
   }
}

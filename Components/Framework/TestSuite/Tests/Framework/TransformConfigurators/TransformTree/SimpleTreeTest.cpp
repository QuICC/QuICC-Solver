/**
 * @file SimpleTreeTest.cpp
 * @brief Tests for a simple transform tree
 */


// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>
#include <string>

// Project includes
//
#include "QuICC/TestSuite/Framework/TransformConfigurators/TransformTree/TestArgs.hpp"
#include "QuICC/Arithmetics/Add.hpp"
#include "QuICC/Arithmetics/Sub.hpp"
#include "QuICC/TransformConfigurators/TransformPath.hpp"
#include "QuICC/TransformConfigurators/TransformTreeTools.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"

TEST_CASE( "SimpleTree", "[Simple]" ){

   // TestSuite namespace
   namespace ns_ts = ::QuICC::TestSuite::Framework;
   // Test namespace
   namespace ns_test = ns_ts::TransformConfigurators::TransformTree;
   // Typedef for Test arguments
   typedef ns_test::TestArgs Args;
   // Typedef 
   namespace ns_QT = ::QuICC::Transform;

   // Set default arguments if required
   if(Args::useDefault)
   {
      // TODO
   }

   Catch::StringMaker<double>::precision = 15;

   // Typedef for used transform
   std::vector<ns_QT::TransformPath> t;

   t.push_back(ns_QT::TransformPath(::QuICC::FieldComponents::Physical::R, ::QuICC::FieldType::VECTOR));
   t.back().addEdge(1, 1, ::QuICC::Arithmetics::Add::id());
   t.back().addEdge(1);
   t.back().addEdge(1, ::QuICC::FieldComponents::Spectral::POL, ::QuICC::Arithmetics::Sub::id());

   t.push_back(ns_QT::TransformPath(::QuICC::FieldComponents::Physical::THETA, ::QuICC::FieldType::VECTOR));
   t.back().addEdge(1, 1, ::QuICC::Arithmetics::Add::id());
   t.back().addEdge(1);
   t.back().addEdge(1, ::QuICC::FieldComponents::Spectral::POL, ::QuICC::Arithmetics::Add::id());

   t.push_back(ns_QT::TransformPath(::QuICC::FieldComponents::Physical::PHI, ::QuICC::FieldType::VECTOR));
   t.back().addEdge(1, 1, ::QuICC::Arithmetics::Add::id());
   t.back().addEdge(1);
   t.back().addEdge(1, ::QuICC::FieldComponents::Spectral::POL, ::QuICC::Arithmetics::Add::id());

   std::map<size_t, std::vector<ns_QT::TransformPath> > mt;
   mt.insert(std::make_pair(::QuICC::PhysicalNames::Velocity::id(), t));

   std::vector<ns_QT::TransformTree> tree;
   ns_QT::TransformTreeTools::generateTrees(tree, mt, ::QuICC::TransformDirection::FORWARD, "simple");
   
}

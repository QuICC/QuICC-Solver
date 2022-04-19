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
#include "QuICC/TransformConfigurators/TransformPath.hpp"
#include "QuICC/TransformConfigurators/TransformTreeTools.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/Transform/Forward/P.hpp"

void createScalarTree(std::vector<::QuICC::Transform::TransformPath>& t);
void createScalarNLTree(std::vector<::QuICC::Transform::TransformPath>& t);

TEST_CASE( "Scalar Tree", "[Scalar]" ){

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

   std::vector<ns_QT::TransformPath> t;
   std::map<size_t, std::vector<ns_QT::TransformPath> > mt;
   std::vector<ns_QT::TransformTree> tree;

   createScalarTree(t);
   mt.insert(std::make_pair(::QuICC::PhysicalNames::Velocity::id(), t));
   ns_QT::TransformTreeTools::generateTrees(tree, mt, ::QuICC::TransformDirection::FORWARD, "scalar");
   t.clear();
   mt.clear();
   tree.clear();

   createScalarTree(t);
   mt.insert(std::make_pair(::QuICC::PhysicalNames::Velocity::id(), t));
   ns_QT::TransformTreeTools::generateTrees(tree, mt, ::QuICC::TransformDirection::FORWARD, "scalar_nl");
   t.clear();
   mt.clear();
   tree.clear();
   
}

void createScalarTree(std::vector<::QuICC::Transform::TransformPath>& t)
{
   namespace ns_QT = ::QuICC::Transform;
   typedef ::QuICC::Transform::TransformPath TP;

   t.push_back(TP(::QuICC::FieldComponents::Physical::SCALAR, ::QuICC::FieldType::SCALAR));
   t.back().addEdge(ns_QT::Forward::P::id());
   t.back().addEdge(ns_QT::Forward::P::id());
   t.back().addEdge(ns_QT::Forward::P::id(), ::QuICC::FieldComponents::Spectral::SCALAR, ::QuICC::Arithmetics::Add::id());
}

void createScalarNLTree(std::vector<::QuICC::Transform::TransformPath>& t)
{
   namespace ns_QT = ::QuICC::Transform;
   typedef ::QuICC::Transform::TransformPath TP;

   t.push_back(TP(::QuICC::FieldComponents::Physical::SCALAR, ::QuICC::FieldType::SCALAR));
   t.back().addEdge(ns_QT::Forward::P::id());
   t.back().addEdge(ns_QT::Forward::P::id());
   t.back().addEdge(ns_QT::Forward::P::id(), ::QuICC::FieldComponents::Spectral::SCALAR, ::QuICC::Arithmetics::Add::id());
}

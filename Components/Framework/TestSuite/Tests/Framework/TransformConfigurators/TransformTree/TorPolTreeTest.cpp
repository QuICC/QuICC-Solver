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
#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/T.hpp"
#include "QuICC/Transform/Forward/Overlaplh.hpp"
#include "QuICC/Transform/Forward/R1.hpp"
#include "QuICC/Transform/Forward/OverlaplhD1.hpp"
#include "QuICC/Transform/Forward/OverlaplhOversinDphi.hpp"
#include "QuICC/Transform/Forward/I4Q.hpp"
#include "QuICC/Transform/Forward/I4S.hpp"
#include "QuICC/Transform/Forward/I2Q.hpp"
#include "QuICC/Transform/Forward/I2S.hpp"
#include "QuICC/Transform/Forward/I2rQ.hpp"
#include "QuICC/Transform/Forward/I2rS.hpp"

void createTorPolTree(std::vector<::QuICC::Transform::TransformPath>& t);
void createTorPolNLTree(std::vector<::QuICC::Transform::TransformPath>& t);

TEST_CASE( "Forward Toroidal/Poloidal Tree", "[FwdTorPol]" ){

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

   createTorPolTree(t);
   mt.insert(std::make_pair(::QuICC::PhysicalNames::Velocity::id(), t));
   ns_QT::TransformTreeTools::generateTrees(tree, mt, ::QuICC::TransformDirection::FORWARD, "torpol");
   t.clear();
   mt.clear();
   tree.clear();

   createTorPolNLTree(t);
   mt.insert(std::make_pair(::QuICC::PhysicalNames::Velocity::id(), t));
   ns_QT::TransformTreeTools::generateTrees(tree, mt, ::QuICC::TransformDirection::FORWARD, "torpol_nl");
   t.clear();
   mt.clear();
   tree.clear();
}

void createTorPolTree(std::vector<::QuICC::Transform::TransformPath>& t)
{
   namespace ns_QT = ::QuICC::Transform;
   typedef ::QuICC::Transform::TransformPath TP;

   // Toroidal

   t.push_back(TP(::QuICC::FieldComponents::Physical::THETA, ::QuICC::FieldType::VECTOR));
   t.back().addEdge(ns_QT::Forward::P::id());
   t.back().addEdge(ns_QT::Forward::OverlaplhOversinDphi::id());
   t.back().addEdge(ns_QT::Forward::P::id(), ::QuICC::FieldComponents::Spectral::TOR, ::QuICC::Arithmetics::Add::id());

   t.push_back(TP(::QuICC::FieldComponents::Physical::PHI, ::QuICC::FieldType::VECTOR));
   t.back().addEdge(ns_QT::Forward::P::id());
   t.back().addEdge(ns_QT::Forward::OverlaplhD1::id());
   t.back().addEdge(ns_QT::Forward::P::id(), ::QuICC::FieldComponents::Spectral::TOR, ::QuICC::Arithmetics::Sub::id());

   // Poloidal

   t.push_back(TP(::QuICC::FieldComponents::Physical::R, ::QuICC::FieldType::VECTOR));
   t.back().addEdge(ns_QT::Forward::P::id());
   t.back().addEdge(ns_QT::Forward::Overlaplh::id());
   t.back().addEdge(ns_QT::Forward::R1::id(), ::QuICC::FieldComponents::Spectral::POL, ::QuICC::Arithmetics::Add::id());
}

void createTorPolNLTree(std::vector<::QuICC::Transform::TransformPath>& t)
{
   namespace ns_QT = ::QuICC::Transform;
   typedef ::QuICC::Transform::TransformPath TP;

   // Toroidal

   t.push_back(TP(::QuICC::FieldComponents::Physical::THETA, ::QuICC::FieldType::VECTOR));
   t.back().addEdge(ns_QT::Forward::P::id());
   t.back().addEdge(ns_QT::Forward::OverlaplhOversinDphi::id());
   t.back().addEdge(ns_QT::Forward::T::id(), ::QuICC::FieldComponents::Spectral::TOR, ::QuICC::Arithmetics::Add::id());

   t.push_back(TP(::QuICC::FieldComponents::Physical::PHI, ::QuICC::FieldType::VECTOR));
   t.back().addEdge(ns_QT::Forward::P::id());
   t.back().addEdge(ns_QT::Forward::OverlaplhD1::id());
   t.back().addEdge(ns_QT::Forward::T::id(), ::QuICC::FieldComponents::Spectral::TOR, ::QuICC::Arithmetics::Sub::id());

   // Poloidal

   t.push_back(TP(::QuICC::FieldComponents::Physical::R, ::QuICC::FieldType::VECTOR));
   t.back().addEdge(ns_QT::Forward::P::id());
   t.back().addEdge(ns_QT::Forward::P::id());
   t.back().addEdge(ns_QT::Forward::I4Q::id(), ::QuICC::FieldComponents::Spectral::POL, ::QuICC::Arithmetics::Sub::id());

   t.push_back(TP(::QuICC::FieldComponents::Physical::THETA, ::QuICC::FieldType::VECTOR));
   t.back().addEdge(ns_QT::Forward::P::id());
   t.back().addEdge(ns_QT::Forward::OverlaplhD1::id());
   t.back().addEdge(ns_QT::Forward::I4S::id(), ::QuICC::FieldComponents::Spectral::POL, ::QuICC::Arithmetics::Add::id());

   t.push_back(TP(::QuICC::FieldComponents::Physical::PHI, ::QuICC::FieldType::VECTOR));
   t.back().addEdge(ns_QT::Forward::P::id());
   t.back().addEdge(ns_QT::Forward::OverlaplhOversinDphi::id());
   t.back().addEdge(ns_QT::Forward::I4S::id(), ::QuICC::FieldComponents::Spectral::POL, ::QuICC::Arithmetics::Add::id());
}

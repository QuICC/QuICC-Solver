/**
 * @file CylinderTransformSteps.cpp
 * @brief Source of the implementation of the physical <-> spectral transform steps in a whole cylinder
 */

// Configuration includes
//

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/TransformConfigurators/CylinderTransformSteps.hpp"

// Project includes
//
#include "QuICC/Arithmetics/Add.hpp"
#include "QuICC/Arithmetics/Sub.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/D1.hpp"
#include "QuICC/Transform/Forward/I2.hpp"
#include "QuICC/Transform/Forward/I4.hpp"
#include "QuICC/Transform/Forward/I4D1.hpp"
#include "QuICC/Transform/Forward/Laplh_1.hpp"
#include "QuICC/Transform/Forward/I4R_1Pm.hpp"
#include "QuICC/Transform/Forward/I4R_1D1R1ZI2.hpp"
#include "QuICC/Transform/Forward/I6R_1Pm.hpp"
#include "QuICC/Transform/Forward/I6R_1D1R1ZI4.hpp"
#include "QuICC/Transform/Forward/I6LaplhZI4D1R1.hpp"
#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/D2.hpp"
#include "QuICC/Transform/Backward/R_1D1R1.hpp"
#include "QuICC/Transform/Backward/R_1Pm.hpp"
#include "QuICC/Transform/Backward/D1ZP.hpp"
#include "QuICC/Transform/Backward/R_1LaplhPm.hpp"
#include "QuICC/Transform/Backward/LaplhZR_1D1R1.hpp"
#include "QuICC/Transform/Backward/D1LaplhZD1R_1D1R1.hpp"

namespace QuICC {

namespace Transform {

   CylinderTransformSteps::CylinderTransformSteps(std::shared_ptr<const SpatialScheme::ISpatialScheme> spScheme)
      : ITransformSteps(spScheme)
   {
   }

   CylinderTransformSteps::~CylinderTransformSteps()
   {
   }

   bool CylinderTransformSteps::applicable(std::shared_ptr<const SpatialScheme::ISpatialScheme> spScheme)
   {
      return spScheme->has(SpatialScheme::Feature::CylinderGeometry);
   }

   std::vector<TransformPath>  CylinderTransformSteps::forwardScalar(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components) const
   {
      assert(components.size() == 1);
      std::vector<TransformPath> transform;

      FieldComponents::Spectral::Id scalId = components.at(0).first;

      transform.push_back(TransformPath(FieldComponents::Physical::SCALAR, FieldType::SCALAR));
      transform.back().addEdge(Forward::P::id());
      transform.back().addEdge(Forward::P::id());
      transform.back().addEdge(Forward::P::id(), scalId, Arithmetics::Add::id());

      return transform;
   }

   std::vector<TransformPath>  CylinderTransformSteps::forwardNLScalar(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components) const
   {
      assert(components.size() == 1);
      std::vector<TransformPath> transform;

      FieldComponents::Spectral::Id scalId = components.at(0).first;

      transform.push_back(TransformPath(FieldComponents::Physical::SCALAR, FieldType::SCALAR));
      transform.back().addEdge(Forward::P::id());
      transform.back().addEdge(Forward::P::id());
      transform.back().addEdge(Forward::P::id(), scalId, Arithmetics::Add::id());

      return transform;
   }

   std::vector<TransformPath>  CylinderTransformSteps::forwardVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         assert(components.size() == 2);
         FieldComponents::Spectral::Id curlId = components.at(0).first;
         int curlFlag = components.at(0).second;
         FieldComponents::Spectral::Id curlcurlId = components.at(1).first;
         int curlcurlFlag = components.at(1).second;

         if(curlFlag == 0 && curlcurlFlag == 0)
         {
            // Compute Toroidal component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::D1::id());
            transform.back().addEdge(Forward::Laplh_1::id(), curlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Laplh_1::id(), curlId, Arithmetics::Sub::id());

            // Compute Poloidal component
            transform.push_back(TransformPath(FieldComponents::Physical::Z, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Laplh_1::id(), curlcurlId, Arithmetics::Sub::id());
         } else
         {
            throw std::logic_error("Requested an unknown vector forward transform");
         }
      } else
      {
         transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), FieldComponents::Spectral::R, Arithmetics::Add::id());

         transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), FieldComponents::Spectral::THETA, Arithmetics::Add::id());

         transform.push_back(TransformPath(FieldComponents::Physical::Z, FieldType::VECTOR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), FieldComponents::Spectral::Z, Arithmetics::Add::id());
      }

      return transform;
   }

   std::vector<TransformPath>  CylinderTransformSteps::forwardNLVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         assert(components.size() == 2);
         FieldComponents::Spectral::Id curlId = components.at(0).first;
         int curlFlag = components.at(0).second;
         FieldComponents::Spectral::Id curlcurlId = components.at(1).first;
         int curlcurlFlag = components.at(1).second;

         // Integrate for curl equation, second order
         if(curlFlag == 0)
         {
            // Compute curl component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::I2::id());
            transform.back().addEdge(Forward::D1::id());
            transform.back().addEdge(Forward::I4R_1Pm::id(), curlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::I2::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::I4R_1D1R1ZI2::id(), curlId, Arithmetics::Sub::id());
         } else
         {
            throw std::logic_error("Requested an unknown vector forward transform");
         }

         // Integrate for double curl equation, fourth order
         if(curlcurlFlag == 0)
         {
            // Compute curlcurl Q component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::I4D1::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::I6R_1D1R1ZI4::id(), curlcurlId, Arithmetics::Add::id());

            // Compute curlcurl S component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::I4D1::id());
            transform.back().addEdge(Forward::D1::id());
            transform.back().addEdge(Forward::I6R_1Pm::id(), curlcurlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::Z, FieldType::VECTOR));
            transform.back().addEdge(Forward::I4::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::I6LaplhZI4D1R1::id(), curlcurlId, Arithmetics::Sub::id());

            // Integrate for double curl equation, second order (typically induction equation)
         } else if(curlcurlFlag == 1)
         {
            throw std::logic_error("Nonlinear transform for induction equation not implemented yet");
            // Compute curlcurl Q component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::I4D1::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::I6R_1D1R1ZI4::id(), curlcurlId, Arithmetics::Add::id());

            // Compute curlcurl S component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::I4D1::id());
            transform.back().addEdge(Forward::D1::id());
            transform.back().addEdge(Forward::I6R_1Pm::id(), curlcurlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::Z, FieldType::VECTOR));
            transform.back().addEdge(Forward::I4::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::I6LaplhZI4D1R1::id(), curlcurlId, Arithmetics::Sub::id());
         } else
         {
            throw std::logic_error("Requested an unknown vector forward transform");
         }
      } else
      {
         transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), FieldComponents::Spectral::R, Arithmetics::Add::id());

         transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), FieldComponents::Spectral::THETA, Arithmetics::Add::id());

         transform.push_back(TransformPath(FieldComponents::Physical::Z, FieldType::VECTOR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), FieldComponents::Spectral::Z, Arithmetics::Add::id());
      }

      return transform;
   }

   std::vector<TransformPath>  CylinderTransformSteps::backwardScalar(const std::map<FieldComponents::Physical::Id,bool>& req) const
   {
      std::vector<TransformPath> transform;

      transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::SCALAR));
      transform.back().addEdge(Backward::P::id());
      transform.back().addEdge(Backward::P::id());
      transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());

      return transform;
   }

   std::vector<TransformPath>  CylinderTransformSteps::backwardGradient(const std::map<FieldComponents::Physical::Id,bool>& req) const
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Backward::D1::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Backward::R_1Pm::id());
         transform.back().addEdge(Backward::D1::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
      }

      if(req.find(FieldComponents::Physical::Z)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::D1::id(), FieldComponents::Physical::Z, Arithmetics::Add::id());
      }

      return transform;
   }

   std::vector<TransformPath>  CylinderTransformSteps::backwardGradient2(const std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>& req) const
   {
      std::vector<TransformPath> transform;
      std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>  pairId;

      throw std::logic_error("Second derivative is not implementated yet!");

      pairId = std::make_pair(this->ss().physical().ONE(),this->ss().physical().ONE());
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
      }

      pairId = std::make_pair(this->ss().physical().ONE(),this->ss().physical().TWO());
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
      }

      pairId = std::make_pair(this->ss().physical().ONE(),this->ss().physical().THREE());
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
      }

      pairId = std::make_pair(this->ss().physical().TWO(),this->ss().physical().TWO());
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
      }

      pairId = std::make_pair(this->ss().physical().TWO(),this->ss().physical().THREE());
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
      }

      pairId = std::make_pair(this->ss().physical().THREE(),this->ss().physical().THREE());
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
      }

      return transform;
   }

   std::vector<TransformPath>  CylinderTransformSteps::backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         if(req.find(FieldComponents::Physical::R)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
            transform.back().addEdge(Backward::R_1Pm::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
            transform.back().addEdge(Backward::D1ZP::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::THETA)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
            transform.back().addEdge(Backward::D1ZP::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Sub::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
            transform.back().addEdge(Backward::R_1Pm::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::D1::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::Z)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
            transform.back().addEdge(Backward::LaplhZR_1D1R1::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::Z, Arithmetics::Sub::id());
         }
      } else
      {
         if(req.find(FieldComponents::Physical::R)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::R, FieldType::VECTOR));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::THETA)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::THETA, FieldType::VECTOR));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::Z)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::Z, FieldType::VECTOR));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::Z, Arithmetics::Add::id());
         }
      }

      return transform;
   }

   std::vector<TransformPath>  CylinderTransformSteps::backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         if(req.find(FieldComponents::Physical::R)->second)
         {
            transform.push_back(TransformPath(id, FieldType::GRADIENT));
            transform.back().addEdge(Backward::D1ZP::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::THETA)->second)
         {
            transform.push_back(TransformPath(id, FieldType::GRADIENT));
            transform.back().addEdge(Backward::R_1Pm::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::Z)->second)
         {
            transform.push_back(TransformPath(id, FieldType::GRADIENT));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id(), FieldComponents::Physical::Z, Arithmetics::Add::id());
         }
      } else
      {
         if(req.find(FieldComponents::Physical::R)->second)
         {
            transform.push_back(TransformPath(id, FieldType::GRADIENT));
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::THETA)->second)
         {
            transform.push_back(TransformPath(id, FieldType::GRADIENT));
            transform.back().addEdge(Backward::R_1Pm::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::Z)->second)
         {
            transform.push_back(TransformPath(id, FieldType::GRADIENT));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id(), FieldComponents::Physical::Z, Arithmetics::Add::id());
         }
      }

      return transform;
   }

   std::vector<TransformPath>  CylinderTransformSteps::backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         if(req.find(FieldComponents::Physical::R)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
            transform.back().addEdge(Backward::D1ZP::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id(), FieldComponents::Physical::R, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
            transform.back().addEdge(Backward::R_1LaplhPm::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Sub::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
            transform.back().addEdge(Backward::R_1Pm::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::D2::id(), FieldComponents::Physical::R, Arithmetics::Sub::id());
         }

         if(req.find(FieldComponents::Physical::THETA)->second)
         {
            // Toroidal part
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
            transform.back().addEdge(Backward::R_1Pm::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::D1::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());

            // Poloidal part
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
            transform.back().addEdge(Backward::D1LaplhZD1R_1D1R1::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());

            // Poloidal part
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
            transform.back().addEdge(Backward::D1ZP::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D2::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::Z)->second)
         {
            // Toroidal part
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
            transform.back().addEdge(Backward::LaplhZR_1D1R1::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::Z, Arithmetics::Sub::id());
         }
      } else
      {
         if(req.find(FieldComponents::Physical::R)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::THETA, FieldType::CURL));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id(), FieldComponents::Physical::R, Arithmetics::Sub::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::Z, FieldType::CURL));
            transform.back().addEdge(Backward::R_1Pm::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::THETA)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::R, FieldType::CURL));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::Z, FieldType::CURL));
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Sub::id());
         }

         if(req.find(FieldComponents::Physical::Z)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::R, FieldType::CURL));
            transform.back().addEdge(Backward::R_1Pm::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::Z, Arithmetics::Sub::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::THETA, FieldType::CURL));
            transform.back().addEdge(Backward::R_1D1R1::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::Z, Arithmetics::Add::id());
         }
      }

      return transform;
   }

   std::vector<TransformPath>  CylinderTransformSteps::backwardDivergence() const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         // The divergence is zero be construction in this case!
         throw std::logic_error("Divergence should not be used in Toroidal/Poloidal expansion");
      } else
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::R, FieldType::DIVERGENCE));
         transform.back().addEdge(Backward::R_1D1R1::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());

         transform.push_back(TransformPath(FieldComponents::Spectral::THETA, FieldType::DIVERGENCE));
         transform.back().addEdge(Backward::R_1Pm::id());
         transform.back().addEdge(Backward::D1::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());

         transform.push_back(TransformPath(FieldComponents::Spectral::Z, FieldType::DIVERGENCE));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::D1::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());
      }

      return transform;
   }

}
}

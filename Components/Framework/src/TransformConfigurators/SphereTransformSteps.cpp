/**
 * @file SphereTransformSteps.cpp
 * @brief Source of the implementation of the physical <-> spectral transform steps in a whole shell
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
#include "QuICC/TransformConfigurators/SphereTransformSteps.hpp"

// Project includes
//
#include "QuICC/Arithmetics/Add.hpp"
#include "QuICC/Arithmetics/Sub.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/T.hpp"
#include "QuICC/Transform/Forward/Laplh_1.hpp"
#include "QuICC/Transform/Forward/R1.hpp"
#include "QuICC/Transform/Forward/Laplh_1D1.hpp"
#include "QuICC/Transform/Forward/Laplh_1Sin_1Dphi.hpp"
#include "QuICC/Transform/Forward/Q4.hpp"
#include "QuICC/Transform/Forward/S4.hpp"
#include "QuICC/Transform/Forward/Q2.hpp"
#include "QuICC/Transform/Forward/S2.hpp"
#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/R_1.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/R_1D1R1.hpp"
#include "QuICC/Transform/Backward/Sin_1Dphi.hpp"
#include "QuICC/Transform/Backward/SLapl.hpp"
#include "QuICC/Transform/Backward/Sin_1D1Sin.hpp"
#include "QuICC/Transform/Backward/Laplh.hpp"

namespace QuICC {

namespace Transform {

   SphereTransformSteps::SphereTransformSteps(std::shared_ptr<const SpatialScheme::ISpatialScheme> spScheme)
      : ITransformSteps(spScheme)
   {
   }

   SphereTransformSteps::~SphereTransformSteps()
   {
   }

   bool SphereTransformSteps::applicable(std::shared_ptr<const SpatialScheme::ISpatialScheme> spScheme)
   {
      return spScheme->has(SpatialScheme::Feature::SphereGeometry);
   }

   std::vector<TransformPath>  SphereTransformSteps::forwardScalar(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components) const
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

   std::vector<TransformPath>  SphereTransformSteps::forwardNLScalar(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components) const
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

   std::vector<TransformPath>  SphereTransformSteps::forwardVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components) const
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
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Laplh_1Sin_1Dphi::id());
            transform.back().addEdge(Forward::P::id(), curlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Laplh_1D1::id());
            transform.back().addEdge(Forward::P::id(), curlId, Arithmetics::Sub::id());

            // Compute Poloidal component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Laplh_1::id());
            transform.back().addEdge(Forward::R1::id(), curlcurlId, Arithmetics::Add::id());
         } else
         {
            throw std::logic_error("Requested an unknown vector forward transform");
         }
      } else
      {
         assert(components.size() == 3);

         transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), FieldComponents::Spectral::R, Arithmetics::Add::id());

         transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), FieldComponents::Spectral::THETA, Arithmetics::Add::id());

         transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), FieldComponents::Spectral::PHI, Arithmetics::Add::id());
      }

      return transform;
   }

   std::vector<TransformPath>  SphereTransformSteps::forwardNLVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         assert(components.size() == 2);
         FieldComponents::Spectral::Id curlId = components.at(0).first;
         int curlFlag = components.at(0).second;
         FieldComponents::Spectral::Id curlcurlId = components.at(1).first;
         int curlcurlFlag = components.at(1).second;

         // Integrate for standard second order equation
         if(curlFlag == 0)
         {
            // Compute curl component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Laplh_1Sin_1Dphi::id());
            transform.back().addEdge(Forward::T::id(), curlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Laplh_1D1::id());
            transform.back().addEdge(Forward::T::id(), curlId, Arithmetics::Sub::id());
         } else
         {
            throw std::logic_error("Requested an unknown vector forward transform");
         }

         // Integrate for standard fourth order spherical equation
         if(curlcurlFlag == 0)
         {
            // Compute curlcurl Q component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Q4::id(), curlcurlId, Arithmetics::Sub::id());

            // Compute curlcurl S component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Laplh_1D1::id());
            transform.back().addEdge(Forward::S4::id(), curlcurlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Laplh_1Sin_1Dphi::id());
            transform.back().addEdge(Forward::S4::id(), curlcurlId, Arithmetics::Add::id());

         // Integrate for second order spherical equation
         } else if(curlcurlFlag == 1)
         {
            // Compute curlcurl Q component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Q2::id(), curlcurlId, Arithmetics::Add::id());

            // Compute curlcurl S component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Laplh_1D1::id());
            transform.back().addEdge(Forward::S2::id(), curlcurlId, Arithmetics::Sub::id());

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Laplh_1Sin_1Dphi::id());
            transform.back().addEdge(Forward::S2::id(), curlcurlId, Arithmetics::Sub::id());
         } else
         {
            throw std::logic_error("Requested an unknown vector forward transform");
         }

      // The following assumes the physical values are obtained from a primitive formulation
      } else
      {
         assert(components.size() == 3);

         transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), FieldComponents::Spectral::R, Arithmetics::Add::id());

         transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), FieldComponents::Spectral::THETA, Arithmetics::Add::id());

         transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), FieldComponents::Spectral::PHI, Arithmetics::Add::id());
      }

      return transform;
   }

   std::vector<TransformPath>  SphereTransformSteps::backwardScalar(const std::map<FieldComponents::Physical::Id,bool>& req) const
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::SCALAR)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::SCALAR));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());
      }

      return transform;
   }

   std::vector<TransformPath>  SphereTransformSteps::backwardGradient(const std::map<FieldComponents::Physical::Id,bool>& req) const
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
         transform.back().addEdge(Backward::R_1::id());
         transform.back().addEdge(Backward::D1::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Backward::R_1::id());
         transform.back().addEdge(Backward::Sin_1Dphi::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
      }

      return transform;
   }

   std::vector<TransformPath>  SphereTransformSteps::backwardGradient2(const std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>& req) const
   {
      std::vector<TransformPath> transform;
      std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>  pairId;

      throw std::logic_error("Second derivative is not implementated yet!");

      pairId = std::make_pair(FieldComponents::Physical::R,FieldComponents::Physical::R);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
      }

      pairId = std::make_pair(FieldComponents::Physical::R,FieldComponents::Physical::THETA);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
      }

      pairId = std::make_pair(FieldComponents::Physical::R,FieldComponents::Physical::PHI);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
      }

      pairId = std::make_pair(FieldComponents::Physical::THETA,FieldComponents::Physical::THETA);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
      }

      pairId = std::make_pair(FieldComponents::Physical::THETA,FieldComponents::Physical::PHI);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
      }

      pairId = std::make_pair(FieldComponents::Physical::PHI,FieldComponents::Physical::PHI);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
      }

      return transform;
   }

   std::vector<TransformPath>  SphereTransformSteps::backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         if(req.find(FieldComponents::Physical::R)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
            transform.back().addEdge(Backward::R_1::id());
            transform.back().addEdge(Backward::Laplh::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::THETA)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::Sin_1Dphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
            transform.back().addEdge(Backward::R_1D1R1::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::PHI)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Sub::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
            transform.back().addEdge(Backward::R_1D1R1::id());
            transform.back().addEdge(Backward::Sin_1Dphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
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

         if(req.find(FieldComponents::Physical::PHI)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::PHI, FieldType::VECTOR));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
         }
      }

      return transform;
   }

   std::vector<TransformPath>  SphereTransformSteps::backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
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
            transform.back().addEdge(Backward::R_1::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::PHI)->second)
         {
            transform.push_back(TransformPath(id, FieldType::GRADIENT));
            transform.back().addEdge(Backward::R_1::id());
            transform.back().addEdge(Backward::Sin_1Dphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
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
            transform.back().addEdge(Backward::R_1::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::PHI)->second)
         {
            transform.push_back(TransformPath(id, FieldType::GRADIENT));
            transform.back().addEdge(Backward::R_1::id());
            transform.back().addEdge(Backward::Sin_1Dphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
         }
      }

      return transform;
   }

   std::vector<TransformPath>  SphereTransformSteps::backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         if(req.find(FieldComponents::Physical::R)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
            transform.back().addEdge(Backward::R_1::id());
            transform.back().addEdge(Backward::Laplh::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::THETA)->second)
         {
            // Toroidal part
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
            transform.back().addEdge(Backward::R_1D1R1::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());

            // Poloidal part
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
            transform.back().addEdge(Backward::SLapl::id());
            transform.back().addEdge(Backward::Sin_1Dphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Sub::id());
         }

         if(req.find(FieldComponents::Physical::PHI)->second)
         {
            // Toroidal part
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
            transform.back().addEdge(Backward::R_1D1R1::id());
            transform.back().addEdge(Backward::Sin_1Dphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());

            // Poloidal part
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
            transform.back().addEdge(Backward::SLapl::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
         }
      } else
      {
         if(req.find(FieldComponents::Physical::R)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::THETA, FieldType::CURL));
            transform.back().addEdge(Backward::R_1::id());
            transform.back().addEdge(Backward::Sin_1Dphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Sub::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::PHI, FieldType::CURL));
            transform.back().addEdge(Backward::R_1::id());
            transform.back().addEdge(Backward::Sin_1D1Sin::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::R, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::THETA)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::R, FieldType::CURL));
            transform.back().addEdge(Backward::R_1::id());
            transform.back().addEdge(Backward::Sin_1Dphi::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::PHI, FieldType::CURL));
            transform.back().addEdge(Backward::R_1D1R1::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Sub::id());
         }

         if(req.find(FieldComponents::Physical::PHI)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::R, FieldType::CURL));
            transform.back().addEdge(Backward::R_1::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Sub::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::THETA, FieldType::CURL));
            transform.back().addEdge(Backward::R_1D1R1::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::PHI, Arithmetics::Add::id());
         }
      }

      return transform;
   }

   std::vector<TransformPath>  SphereTransformSteps::backwardDivergence() const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         // The divergence is zero be construction in this case!
         throw std::logic_error("Divergence should not be used in Toroidal/Poloidal expansion");
      } else
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::R, FieldType::DIVERGENCE));
         transform.back().addEdge(Backward::D1::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());

         transform.push_back(TransformPath(FieldComponents::Spectral::THETA, FieldType::DIVERGENCE));
         transform.back().addEdge(Backward::R_1::id());
         transform.back().addEdge(Backward::Sin_1D1Sin::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());

         transform.push_back(TransformPath(FieldComponents::Spectral::PHI, FieldType::DIVERGENCE));
         transform.back().addEdge(Backward::R_1::id());
         transform.back().addEdge(Backward::Sin_1Dphi::id());
         transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());
      }

      return transform;
   }

}
}

/**
 * @file CartesianTransformSteps.cpp
 * @brief Source of the implementation of the cartesian transform steps
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
#include "QuICC/TransformConfigurators/CartesianTransformSteps.hpp"

// Project includes
//
#include "QuICC/Arithmetics/Add.hpp"
#include "QuICC/Arithmetics/Sub.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/P0.hpp"
#include "QuICC/Transform/Forward/Pm.hpp"
#include "QuICC/Transform/Forward/D1.hpp"
#include "QuICC/Transform/Forward/D1ZP0.hpp"
#include "QuICC/Transform/Forward/DfLaplh_1.hpp"
#include "QuICC/Transform/Forward/Laplh.hpp"
#include "QuICC/Transform/Forward/Laplh_1.hpp"
#include "QuICC/Transform/Forward/I2.hpp"
#include "QuICC/Transform/Forward/I2D1.hpp"
#include "QuICC/Transform/Forward/I2ZI2D1.hpp"
#include "QuICC/Transform/Forward/I4.hpp"
#include "QuICC/Transform/Forward/I4D1.hpp"
#include "QuICC/Transform/Forward/I4D1ZI2.hpp"
#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/P0.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/D2.hpp"
#include "QuICC/Transform/Backward/Laplh.hpp"
#include "QuICC/Transform/Backward/DfLaplh.hpp"
#include "QuICC/Transform/Backward/DsLaplh.hpp"

namespace QuICC {

namespace Transform {

   CartesianTransformSteps::CartesianTransformSteps(std::shared_ptr<const SpatialScheme::ISpatialScheme> spScheme)
      : ITransformSteps(spScheme)
   {
   }

   CartesianTransformSteps::~CartesianTransformSteps()
   {
   }

   bool CartesianTransformSteps::applicable(std::shared_ptr<const SpatialScheme::ISpatialScheme> spScheme)
   {
      return spScheme->has(SpatialScheme::Feature::CartesianGeometry);
   }

   std::vector<TransformPath>  CartesianTransformSteps::forwardScalar(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components) const
   {
      assert(components.size() == 1);
      std::vector<TransformPath> transform;

      FieldComponents::Spectral::Id scalId = components.at(0).first;

      transform.push_back(TransformPath(FieldComponents::Physical::SCALAR, FieldType::SCALAR));
      if(this->ss().dimension() > 2) transform.back().addEdge(Forward::P::id());
      if(this->ss().dimension() > 1) transform.back().addEdge(Forward::P::id());
      transform.back().addEdge(Forward::P::id(), scalId, Arithmetics::Add::id());

      return transform;
   }

   std::vector<TransformPath>  CartesianTransformSteps::forwardNLScalar(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components) const
   {
      assert(components.size() == 1);
      std::vector<TransformPath> transform;

      FieldComponents::Spectral::Id scalId = components.at(0).first;

      transform.push_back(TransformPath(FieldComponents::Physical::SCALAR, FieldType::SCALAR));
      if(this->ss().dimension() > 2) transform.back().addEdge(Forward::P::id());
      if(this->ss().dimension() > 1) transform.back().addEdge(Forward::P::id());
      transform.back().addEdge(Forward::P::id(), scalId, Arithmetics::Add::id());

      return transform;
   }

   std::vector<TransformPath>  CartesianTransformSteps::forwardVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components) const
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
            // Extract Toroidal component
            transform.push_back(TransformPath(FieldComponents::Physical::X, FieldType::VECTOR));
            transform.back().addEdge(Forward::D1::id());
            transform.back().addEdge(Forward::Laplh_1::id());
            transform.back().addEdge(Forward::P::id(), curlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::Y, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::DfLaplh_1::id());
            transform.back().addEdge(Forward::P::id(), curlId, Arithmetics::Sub::id());

            // Extract Poloidal component
            transform.push_back(TransformPath(FieldComponents::Physical::Z, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Laplh_1::id());
            transform.back().addEdge(Forward::P::id(), curlcurlId, Arithmetics::Sub::id());

            // Extract X mean into Toroidal component
            transform.push_back(TransformPath(FieldComponents::Physical::X, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P0::id());
            transform.back().addEdge(Forward::P::id(), curlId, Arithmetics::Add::id());

            // Extract Y mean into Poloidal component
            transform.push_back(TransformPath(FieldComponents::Physical::Y, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P0::id());
            transform.back().addEdge(Forward::P::id(), curlcurlId, Arithmetics::Add::id());
         } else
         {
            throw std::logic_error("Requested an unknown vector forward transform");
         }
      } else
      {
         assert(components.size() == 3);
         transform.push_back(TransformPath(this->ss().physical().ONE(), FieldType::VECTOR));
         if(this->ss().dimension() > 2) transform.back().addEdge(Forward::P::id());
         if(this->ss().dimension() > 1) transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), this->ss().spectral().ONE(), Arithmetics::Add::id());

         if(this->ss().dimension() > 1)
         {
            transform.push_back(TransformPath(this->ss().physical().TWO(), FieldType::VECTOR));
            if(this->ss().dimension() > 2) transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id(), this->ss().spectral().TWO(), Arithmetics::Add::id());
         }

         if(this->ss().dimension() == 3)
         {
            transform.push_back(TransformPath(this->ss().physical().THREE(), FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id(), this->ss().spectral().THREE(), Arithmetics::Add::id());
         }
      }

      return transform;
   }

   std::vector<TransformPath>  CartesianTransformSteps::forwardNLVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components) const
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
            // Compute curl component and mean in X
            transform.push_back(TransformPath(FieldComponents::Physical::X, FieldType::VECTOR));
            transform.back().addEdge(Forward::D1ZP0::id());
            transform.back().addEdge(Forward::Pm::id());
            transform.back().addEdge(Forward::I2::id(), curlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::Y, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::D1::id());
            transform.back().addEdge(Forward::I2::id(), curlId, Arithmetics::Sub::id());

            // Integrate for curl of nonlinear and mean Y (induction)
         } else if(curlFlag == 1)
         {
            transform.push_back(TransformPath(FieldComponents::Physical::X, FieldType::VECTOR));
            transform.back().addEdge(Forward::D1ZP0::id());
            transform.back().addEdge(Forward::Pm::id());
            transform.back().addEdge(Forward::I2ZI2D1::id(), curlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::Y, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::D1::id());
            transform.back().addEdge(Forward::I2::id(), curlId, Arithmetics::Sub::id());

         } else
         {
            throw std::logic_error("Requested an unknown vector forward transform");
         }

         // Integrate for standard fourth order double curl equation
         if(curlcurlFlag == 0)
         {
            // Compute curlcurl with Dz component and mean in Y
            transform.push_back(TransformPath(FieldComponents::Physical::X, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::D1::id());
            transform.back().addEdge(Forward::I4D1::id(), curlcurlId, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Physical::Y, FieldType::VECTOR));
            transform.back().addEdge(Forward::D1ZP0::id());
            transform.back().addEdge(Forward::Pm::id());
            transform.back().addEdge(Forward::I4D1ZI2::id(), curlcurlId, Arithmetics::Add::id());

            // Compute curlcurl without Dz component
            transform.push_back(TransformPath(FieldComponents::Physical::Z, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Laplh::id());
            transform.back().addEdge(Forward::I4::id(), curlcurlId, Arithmetics::Sub::id());

            // Integrate for double curl of nonlinear term and mean X (induction)
         } else if(curlcurlFlag == 1)
         {
            // Compute curlcurl with Dz component
            transform.push_back(TransformPath(FieldComponents::Physical::X, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::D1::id());
            transform.back().addEdge(Forward::I2D1::id(), curlcurlId, Arithmetics::Sub::id());

            transform.push_back(TransformPath(FieldComponents::Physical::Y, FieldType::VECTOR));
            transform.back().addEdge(Forward::D1ZP0::id());
            transform.back().addEdge(Forward::Pm::id());
            transform.back().addEdge(Forward::I2D1::id(), curlcurlId, Arithmetics::Sub::id());

            // Compute curlcurl without Dz component
            transform.push_back(TransformPath(FieldComponents::Physical::Z, FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::Laplh::id());
            transform.back().addEdge(Forward::I2::id(), curlcurlId, Arithmetics::Add::id());

         } else
         {
            throw std::logic_error("Requested an unknown vector forward transform");
         }

         // The following assumes the physical values are obtained from a Toroidal/Poloidal decomposition
      } else
      {
         assert(components.size() == 3);
         transform.push_back(TransformPath(this->ss().physical().ONE(), FieldType::VECTOR));
         if(this->ss().dimension() > 2) transform.back().addEdge(Forward::P::id());
         if(this->ss().dimension() > 1) transform.back().addEdge(Forward::P::id());
         transform.back().addEdge(Forward::P::id(), this->ss().spectral().ONE(), Arithmetics::Add::id());

         if(this->ss().dimension() > 1)
         {
            transform.push_back(TransformPath(this->ss().physical().TWO(), FieldType::VECTOR));
            if(this->ss().dimension() > 2) transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id(), this->ss().spectral().TWO(), Arithmetics::Add::id());
         }

         if(this->ss().dimension() == 3)
         {
            transform.push_back(TransformPath(this->ss().physical().THREE(), FieldType::VECTOR));
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id());
            transform.back().addEdge(Forward::P::id(), this->ss().spectral().THREE(), Arithmetics::Add::id());
         }
      }

      return transform;
   }

   std::vector<TransformPath>  CartesianTransformSteps::backwardScalar(const std::map<FieldComponents::Physical::Id,bool>& req) const
   {
      std::vector<TransformPath> transform;

      transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::SCALAR));
      if(this->ss().dimension() > 2) transform.back().addEdge(Backward::P::id());
      if(this->ss().dimension() > 1) transform.back().addEdge(Backward::P::id());
      transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());

      return transform;
   }

   std::vector<TransformPath>  CartesianTransformSteps::backwardGradient(const std::map<FieldComponents::Physical::Id,bool>& req) const
   {
      std::vector<TransformPath> transform;

      if(req.find(this->ss().physical().ONE())->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         if(this->ss().dimension() > 1) transform.back().addEdge(Backward::D1::id());
         if(this->ss().dimension() > 2) transform.back().addEdge(Backward::P::id());
         if(this->ss().dimension() > 1)
         {
            transform.back().addEdge(Backward::P::id(), this->ss().physical().ONE(), Arithmetics::Add::id());
         } else
         {
            transform.back().addEdge(Backward::D1::id(), this->ss().physical().ONE(), Arithmetics::Add::id());
         }
      }

      if(this->ss().dimension() > 1 && req.find(this->ss().physical().TWO())->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Backward::P::id());
         if(this->ss().dimension() > 2)
         {
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), this->ss().physical().TWO(), Arithmetics::Add::id());
         } else
         {
            transform.back().addEdge(Backward::D1::id(), this->ss().physical().TWO(), Arithmetics::Add::id());
         }
      }

      if(this->ss().dimension() > 2 && req.find(this->ss().physical().THREE())->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::P::id());
         transform.back().addEdge(Backward::D1::id(), this->ss().physical().THREE(), Arithmetics::Add::id());
      }

      return transform;
   }

   std::vector<TransformPath>  CartesianTransformSteps::backwardGradient2(const std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>& req) const
   {
      std::vector<TransformPath> transform;
      std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>  pairId;

      pairId = std::make_pair(this->ss().physical().ONE(),this->ss().physical().ONE());
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         if(this->ss().dimension() > 1) transform.back().addEdge(Backward::D2::id());
         if(this->ss().dimension() > 2) transform.back().addEdge(Backward::P::id());
         if(this->ss().dimension() > 1)
         {
            transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
         } else
         {
            transform.back().addEdge(Backward::D2::id(), pairId, Arithmetics::Add::id());
         }
      }

      if(this->ss().dimension() > 1)
      {
         pairId = std::make_pair(this->ss().physical().ONE(),this->ss().physical().TWO());
         if(req.find(pairId)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
            transform.back().addEdge(Backward::D1::id());
            if(this->ss().dimension() > 2)
            {
               transform.back().addEdge(Backward::D1::id());
               transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
            } else
            {
               transform.back().addEdge(Backward::D1::id(), pairId, Arithmetics::Add::id());
            }
         }
      }

      if(this->ss().dimension() > 2)
      {
         pairId = std::make_pair(this->ss().physical().ONE(),this->ss().physical().THREE());
         if(req.find(pairId)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id(), pairId, Arithmetics::Add::id());
         }
      }

      if(this->ss().dimension() > 1)
      {
         pairId = std::make_pair(this->ss().physical().TWO(),this->ss().physical().TWO());
         if(req.find(pairId)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
            transform.back().addEdge(Backward::P::id());
            if(this->ss().dimension() > 2)
            {
               transform.back().addEdge(Backward::D2::id());
               transform.back().addEdge(Backward::P::id(), pairId, Arithmetics::Add::id());
            } else
            {
               transform.back().addEdge(Backward::D2::id(), pairId, Arithmetics::Add::id());
            }
         }
      }

      if(this->ss().dimension() > 2)
      {
         pairId = std::make_pair(this->ss().physical().TWO(),this->ss().physical().THREE());
         if(req.find(pairId)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::D1::id(), pairId, Arithmetics::Add::id());
         }

         pairId = std::make_pair(this->ss().physical().THREE(),this->ss().physical().THREE());
         if(req.find(pairId)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D2::id(), pairId, Arithmetics::Add::id());
         }
      }

      return transform;
   }

   std::vector<TransformPath>  CartesianTransformSteps::backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         if(req.find(FieldComponents::Physical::X)->second)
         {
            // Toroidal contribution
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id(), FieldComponents::Physical::X, Arithmetics::Add::id());

            // Mean contribution
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P0::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::X, Arithmetics::Add::id());

            // Poloidal contribution
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::X, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::Y)->second)
         {
            // Toroidal contribution
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::Y, Arithmetics::Sub::id());

            // Poloidal contribution
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id(), FieldComponents::Physical::Y, Arithmetics::Add::id());

            // Mean contribution
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P0::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::Y, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::Z)->second)
         {
            // Poloidal contributions
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::Laplh::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::Z, Arithmetics::Sub::id());
         }
      } else
      {
         if(req.find(this->ss().physical().ONE())->second)
         {
            transform.push_back(TransformPath(this->ss().spectral().ONE(), FieldType::VECTOR));
            if(this->ss().dimension() > 1) transform.back().addEdge(Backward::P::id());
            if(this->ss().dimension() > 2) transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), this->ss().physical().ONE(), Arithmetics::Add::id());
         }

         if(this->ss().dimension() > 1)
         {
            if(req.find(this->ss().physical().TWO())->second)
            {
               transform.push_back(TransformPath(this->ss().spectral().TWO(), FieldType::VECTOR));
               transform.back().addEdge(Backward::P::id());
               if(this->ss().dimension() > 2) transform.back().addEdge(Backward::P::id());
               transform.back().addEdge(Backward::P::id(), this->ss().physical().TWO(), Arithmetics::Add::id());
            }
         }

         if(this->ss().dimension() > 2)
         {
            if(req.find(this->ss().physical().THREE())->second)
            {
               transform.push_back(TransformPath(this->ss().spectral().THREE(), FieldType::VECTOR));
               transform.back().addEdge(Backward::P::id());
               transform.back().addEdge(Backward::P::id());
               transform.back().addEdge(Backward::P::id(), this->ss().physical().THREE(), Arithmetics::Add::id());
            }
         }
      }

      return transform;
   }

   std::vector<TransformPath>  CartesianTransformSteps::backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         throw std::logic_error("Vector gradient for Toroidal/Poloidal field not implemented");
      } else
      {
         if(req.find(this->ss().physical().ONE())->second)
         {
            transform.push_back(TransformPath(id, FieldType::GRADIENT));
            if(this->ss().dimension() > 1) transform.back().addEdge(Backward::D1::id());
            if(this->ss().dimension() > 2) transform.back().addEdge(Backward::P::id());
            if(this->ss().dimension() > 1)
            {
               transform.back().addEdge(Backward::P::id(), this->ss().physical().ONE(), Arithmetics::Add::id());
            } else
            {
               transform.back().addEdge(Backward::D1::id(), this->ss().physical().ONE(), Arithmetics::Add::id());
            }
         }

         if(this->ss().dimension() > 1)
         {
            if(req.find(this->ss().physical().TWO())->second)
            {
               transform.push_back(TransformPath(id, FieldType::GRADIENT));
               transform.back().addEdge(Backward::P::id());
               if(this->ss().dimension() > 2)
               {
                  transform.back().addEdge(Backward::D1::id());
                  transform.back().addEdge(Backward::P::id(), this->ss().physical().TWO(), Arithmetics::Add::id());
               } else
               {
                  transform.back().addEdge(Backward::D1::id(), this->ss().physical().TWO(), Arithmetics::Add::id());
               }
            }
         }

         if(this->ss().dimension() > 2)
         {
            if(req.find(this->ss().physical().THREE())->second)
            {
               transform.push_back(TransformPath(id, FieldType::GRADIENT));
               transform.back().addEdge(Backward::P::id());
               transform.back().addEdge(Backward::P::id());
               transform.back().addEdge(Backward::D1::id(), this->ss().physical().THREE(), Arithmetics::Add::id());
            }
         }
      }

      return transform;
   }

   std::vector<TransformPath>  CartesianTransformSteps::backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req) const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         if(req.find(FieldComponents::Physical::X)->second)
         {
            // Toroidal contribution
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::X, Arithmetics::Add::id());

            // Poloidal contributions
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
            transform.back().addEdge(Backward::D2::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id(), FieldComponents::Physical::X, Arithmetics::Sub::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::DsLaplh::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::X, Arithmetics::Sub::id());

            // Mean contribution
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P0::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::X, Arithmetics::Sub::id());
         }

         if(req.find(FieldComponents::Physical::Y)->second)
         {
            // Toroidal contribution
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id(), FieldComponents::Physical::Y, Arithmetics::Add::id());

            // Mean contribution
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P0::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::Y, Arithmetics::Add::id());

            // Poloidal contributions
            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
            transform.back().addEdge(Backward::D2::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::Y, Arithmetics::Add::id());

            transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::DfLaplh::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::Y, Arithmetics::Add::id());
         }

         if(req.find(FieldComponents::Physical::Z)->second)
         {
            transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::Laplh::id());
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::Z, Arithmetics::Sub::id());
         }
      } else
      {
         if(this->ss().dimension() != 3)
         {
            throw std::logic_error("Curl cannot be computed!");
         }

         if(req.find(this->ss().physical().ONE())->second)
         {
            transform.push_back(TransformPath(this->ss().spectral().TWO(), FieldType::CURL));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id(), this->ss().physical().ONE(), Arithmetics::Sub::id());

            transform.push_back(TransformPath(this->ss().spectral().THREE(), FieldType::CURL));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), this->ss().physical().ONE(), Arithmetics::Add::id());
         }

         if(req.find(this->ss().physical().TWO())->second)
         {
            transform.push_back(TransformPath(this->ss().spectral().ONE(), FieldType::CURL));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id(), this->ss().physical().TWO(), Arithmetics::Add::id());

            transform.push_back(TransformPath(this->ss().spectral().THREE(), FieldType::CURL));
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), this->ss().physical().TWO(), Arithmetics::Sub::id());
         }

         if(req.find(this->ss().physical().THREE())->second)
         {
            transform.push_back(TransformPath(this->ss().spectral().ONE(), FieldType::CURL));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id(), this->ss().physical().THREE(), Arithmetics::Sub::id());

            transform.push_back(TransformPath(this->ss().spectral().TWO(), FieldType::CURL));
            transform.back().addEdge(Backward::D1::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id(), this->ss().physical().THREE(), Arithmetics::Add::id());
         }
      }

      return transform;
   }

   std::vector<TransformPath>  CartesianTransformSteps::backwardDivergence() const
   {
      std::vector<TransformPath> transform;

      if(this->ss().formulation() == VectorFormulation::TORPOL)
      {
         // The divergence is zero be construction in this case!
         throw std::logic_error("Divergence should not be used in Toroidal/Poloidal expansion");
      } else
      {
         transform.push_back(TransformPath(this->ss().spectral().ONE(), FieldType::DIVERGENCE));
         if(this->ss().dimension() > 1) transform.back().addEdge(Backward::D1::id());
         if(this->ss().dimension() > 2) transform.back().addEdge(Backward::P::id());
         if(this->ss().dimension() > 1)
         {
            transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());
         } else
         {
            transform.back().addEdge(Backward::D1::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());
         }

         if(this->ss().dimension() > 1)
         {
            transform.push_back(TransformPath(this->ss().spectral().TWO(), FieldType::DIVERGENCE));
            transform.back().addEdge(Backward::P::id());
            if(this->ss().dimension() > 2)
            {
               transform.back().addEdge(Backward::D1::id());
               transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());
            } else
            {
               transform.back().addEdge(Backward::D1::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());
            }
         }

         if(this->ss().dimension() > 2)
         {
            transform.push_back(TransformPath(this->ss().spectral().THREE(), FieldType::DIVERGENCE));
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::P::id());
            transform.back().addEdge(Backward::D1::id(), FieldComponents::Physical::SCALAR, Arithmetics::Add::id());
         }
      }

      return transform;
   }

}
}

/** 
 * @file IPhysicalModel.cpp
 * @brief source of the implementation of a physical model
 */

// System includes
//
#include <string>
#include <vector>
#include <set>
#include <memory>

// Project includes
//
#include "QuICC/Model/IPhysicalModel.hpp"
#include "QuICC/SpatialScheme/Feature.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "QuICC/Arithmetics/registerAll.hpp"
#include "QuICC/ModelOperator/registerAll.hpp"
#include "QuICC/ModelOperatorBoundary/registerAll.hpp"
#include "QuICC/NonDimensional/registerAll.hpp"
#include "QuICC/PhysicalNames/registerAll.hpp"
#include "QuICC/RuntimeStatus/registerAll.hpp"
#include "QuICC/SolveTiming/registerAll.hpp"
#include "QuICC/PseudospectralTag/registerAll.hpp"

namespace QuICC {

namespace Model {

   void IPhysicalModel::init()
   {
      this->registerNames();
   }

   void IPhysicalModel::registerNames()
   {
      // Arithmetics names
      Arithmetics::registerAll();
      // ModelOperator names
      ModelOperator::registerAll();
      // ModelOperatorBoundary names
      ModelOperatorBoundary::registerAll();
      // NonDimensional names
      NonDimensional::registerAll();
      // Physical names
      PhysicalNames::registerAll();
      // RuntimeStatus names
      RuntimeStatus::registerAll();
      // SolveTiming names
      SolveTiming::registerAll();
      // PseudospectralTag names
      PseudospectralTag::registerAll();
   }

   std::vector<std::size_t> IPhysicalModel::extraFieldIds() const
   {
      std::vector<std::size_t> extra;

      return extra;
   }

   std::map<std::string, std::map<std::string,int> > IPhysicalModel::configTags() const
   {
      std::map<std::string, std::map<std::string,int> > tags;

      return tags;
   }

   void IPhysicalModel::configure(const std::set<SpatialScheme::Feature>& f)
   {
      // Propagate Galerkin flag
      this->mpBackend->enableGalerkin(f.count(SpatialScheme::Feature::GalerkinBasis));

      // Propagate split 4th order equations flag
      this->mpBackend->enableSplitEquation(f.count(SpatialScheme::Feature::SplitFourthOrder));
   }

   const IModelBackend& IPhysicalModel::backend() const
   {
      return *this->mpBackend;
   }

   std::shared_ptr<IModelBackend> IPhysicalModel::spBackend() const
   {
      return this->mpBackend;
   }

}
}

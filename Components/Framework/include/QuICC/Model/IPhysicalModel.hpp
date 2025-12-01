/** 
 * @file IPhysicalModel.hpp
 * @brief Interface for implementation of a physical model
 */

#ifndef QUICC_MODEL_IPHYSICALMODEL_HPP
#define QUICC_MODEL_IPHYSICALMODEL_HPP

// System includes
//
#include <string>
#include <vector>
#include <set>
#include <memory>

// Project includes
//
#include "QuICC/Enums/VectorFormulation.hpp"
#include "QuICC/SpatialScheme/Feature.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"

namespace QuICC {

namespace Model {

   /**
    * @brief Interface for the implementation of a physical models
    */
   class IPhysicalModel
   {
      public:
         /**
          * @brief Constructor
          */
         IPhysicalModel() = default;

         /**
          * @brief Destructor
          */
         virtual ~IPhysicalModel() = default;

         /**
          * @brief Tune the spatial scheme (for example change mesher)
          *
          * @param spScheme   Spatial scheme
          */
         virtual void tuneScheme(std::shared_ptr<SpatialScheme::ISpatialScheme> spScheme) {};

         /**
          * @brief Initialize model
          */
         virtual void init();

         /**
          * @brief Formulation used for vector fields
          */
         virtual VectorFormulation::Id SchemeFormulation() = 0;

         /**
          * @brief Version string of model
          */
         virtual std::string version() const = 0;

         /**
          * @brief XML configuration tags for model
          */
         virtual std::map<std::string,std::map<std::string,int> > configTags() const;

         /**
          * @brief Configure additional features set at run time
          */
         virtual void configure(const std::set<SpatialScheme::Feature>& f);

         /**
          * @brief Add extra field IDs (example: imposed fields)
          */
         virtual std::vector<std::size_t> extraFieldIds() const;

         /**
          * @brief Get model generator
          */
         const IModelBackend& backend() const;

         /**
          * @brief Get model generator
          */
         std::shared_ptr<IModelBackend> spBackend() const;

         /**
          * @brief Interface to adding ASCII output file
          */
         template <typename T, typename TApp> std::shared_ptr<T> enableAsciiFile(const std::string tag, const std::string prefix, const std::size_t id, std::shared_ptr<TApp> spSim);

         /**
          * @brief Interface to adding ASCII output file (anelastic case)
          */
         template <typename T, typename TApp> std::shared_ptr<T> enableAsciiFile(const std::string tag, 
                                                                                 const std::string prefix, 
                                                                                 const std::size_t id, 
                                                                                 std::shared_ptr<TApp> spSim,
                                                                                 std::vector<std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction>> pF);
      protected:
         /**
          * @brief Register Named IDs needed for simulation
          */
         virtual void registerNames();

         /**
          * @brief Model generator
          */
         std::shared_ptr<IModelBackend> mpBackend;

      private:
   };

   template <typename T, typename TApp> std::shared_ptr<T> IPhysicalModel::enableAsciiFile(const std::string tag, const std::string prefix, const std::size_t id, std::shared_ptr<TApp> spSim)
   {
      if(spSim->config().model(tag).at("enable"))
      {
         auto spFile = std::make_shared<T>(prefix, spSim->ss().tag());
         spFile->expect(id);
         if((spSim->config().model(tag).count("numbered") > 0) && spSim->config().model(tag).at("numbered"))
         {
            spFile->numberOutput();
         }
         if(spSim->config().model(tag).count("only_every") > 0)
         {
            spFile->onlyEvery(spSim->config().model(tag).at("only_every"));
         }
         spSim->addAsciiOutputFile(spFile);

         return spFile;
      }
      else
      {
         return nullptr;
      }
   }

   // anelastic case overload:
   // this version accepts a vector of pointers to DenseSM profiles
   template <typename T, typename TApp> std::shared_ptr<T> IPhysicalModel::enableAsciiFile(const std::string tag, 
                                                                                           const std::string prefix, 
                                                                                           const std::size_t id, 
                                                                                           std::shared_ptr<TApp> spSim,
                                                                                           std::vector<std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction>> pF)
   {
      if(spSim->config().model(tag).at("enable"))
      {
         auto spFile = std::make_shared<T>(prefix, spSim->ss().tag(), pF);
         spFile->expect(id);
         if((spSim->config().model(tag).count("numbered") > 0) && spSim->config().model(tag).at("numbered"))
         {
            spFile->numberOutput();
         }
         if(spSim->config().model(tag).count("only_every") > 0)
         {
            spFile->onlyEvery(spSim->config().model(tag).at("only_every"));
         }
         spSim->addAsciiOutputFile(spFile);

         return spFile;
      }
      else
      {
         return nullptr;
      }
   }

}
}

#endif // QUICC_MODEL_IPHYSICALMODEL_HPP

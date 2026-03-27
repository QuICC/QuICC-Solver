/**
 * @file IFieldEquation._decl.hpp
 * @brief Base building block for the implementation of an equation
 */

#ifndef QUICC_EQUATIONS_IFIELDEQUATION_DECL_HPP
#define QUICC_EQUATIONS_IFIELDEQUATION_DECL_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Equations/IEquation.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "Arithmetics/Basic.hpp"
#include "QuICC/Equations/SolutionUpdater.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Base building block for the implementation of an equation
    */
   class IFieldEquation : public IEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Equation parameters
          * @param spScheme   Spatial scheme
          * @param spBackend  Model backend
          */
         explicit IFieldEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Equation parameters
          * @param spScheme   Spatial scheme
          * @param spBackend  Model backend
          * @param spOptions  Additional options
          */
         explicit IFieldEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend, std::shared_ptr<EquationOptions> spOptions);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IFieldEquation() = default;

         /**
          * @brief Apply generic constraint on spectral data
          *
          * @param compId  ID of the spectral component
          * @param timeId  Timing of the constraint
          * @return contraint changed solution?
          */
         virtual bool applyConstraint(FieldComponents::Spectral::Id compId, const std::size_t timeId);

         /**
          * @brief Get correction corresponding to generic constraint on spectral data
          *
          * @param compId  ID of the spectral component
          * @param timeId  Timing of the constraint
          */
         virtual std::vector<std::tuple<MHDVariant,int,int,int>> correctionConstraint(FieldComponents::Spectral::Id compId, const std::size_t timeId);

         /**
          * @brief Get solution updater kernel
          *
          * @param compId  ID of the spectral component
          */
         std::shared_ptr<SolutionUpdater> solutionUpdater(FieldComponents::Spectral::Id compId) const;

         /**
          * @brief Get source term kernel
          *
          * @param compId  ID of the spectral component
          */
         Spectral::Kernel::SharedISpectralKernel sourceKernel(FieldComponents::Spectral::Id compId) const;

         /**
          * @brief Generic source term implementation
          *
          * @param compId  ID of the spectral component
          * @param i       Fastest index
          * @param j       Second index
          * @param k       Slowest index
          */
         virtual MHDVariant sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

         /**
          * @brief Get boundary value kernel
          *
          * @param compId  ID of the spectral component
          */
         Spectral::Kernel::SharedISpectralKernel boundaryKernel(FieldComponents::Spectral::Id compId) const;

         /**
          * @brief Generic boundary value implementation
          *
          * @param compId  ID of the spectral component
          * @param i       Fastest index
          * @param j       Second index
          * @param k       Slowest index
          */
         virtual MHDVariant boundaryValue(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

      protected:
      private:

   };

   /// Typedef for a smart IFieldEquation
   typedef std::shared_ptr<IFieldEquation> SharedIFieldEquation;

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_IFIELDEQUATION_DECL_HPP

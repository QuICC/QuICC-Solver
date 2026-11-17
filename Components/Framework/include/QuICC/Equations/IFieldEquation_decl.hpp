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
#include "QuICC/Equations/details/StoreSolutionFunctor.hpp"

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
          * @brief Generic source term implementation
          *
          * @param compId  ID of the spectral component
          * @param i       Fastest index
          * @param j       Second index
          * @param k       Slowest index
          */
         virtual MHDVariant sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

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
         /**
          * @brief Templated passthrough update the stored value with the solver solution
          */
         virtual MHDVariant updateStoredSolution(const MHDVariant newData, FieldComponents::Spectral::Id compId, const int i, const int j, const int k);

         /**
          * @brief Transfer solver solution to equation unknown
          *
          * @param field   Scalar or vector field
          * @param compId  Component ID
          * @param storage Solver solution
          * @param matIdx  Index of the given data
          * @param start   Start index for the storage
          */
         template <typename TData, typename TField> void storeSolutionImpl(TField& field, FieldComponents::Spectral::Id compId, const TData& storage, const int matIdx, const int start);

         template <CouplingIndexType IndexType> friend class details::StoreSolutionFunctor;
      private:

   };

   /// Typedef for a smart IFieldEquation
   typedef std::shared_ptr<IFieldEquation> SharedIFieldEquation;

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_IFIELDEQUATION_DECL_HPP

#include "QuICC/ScalarFields/ScalarFieldSetup.hpp"
#include "QuICC/ScalarFields/FlatScalarField.hpp"

namespace details {

   enum class SetupType {
      UniformUp,
      UniformDown,
      TriangularUp,
      TriangularDown
   };

   std::shared_ptr<QuICC::Datatypes::ScalarFieldSetup> createSetup(const SetupType type, const int dim1D, const int dim3D);

   template <typename T> std::shared_ptr<QuICC::Datatypes::FlatScalarField<T>> createFlatScalarField(const SetupType type, const int dim1D, const int dim3D);

   template <typename T> std::shared_ptr<QuICC::Datatypes::FlatScalarField<T>> createFlatScalarField(const SetupType type, const int dim1D, const int dim3D)
   {
      auto spSetup = details::createSetup(details::SetupType::UniformUp, dim1D, dim3D);

      auto spField = std::make_shared<QuICC::Datatypes::FlatScalarField<T>>(spSetup);
      auto&& data = spField->rData();

      data.setConstant(-42.42);
      for(int k = 0; k < spSetup->nBlock(); k++)
      {
         for(int j = 0; j < spSetup->blockCols(k); j++)
         {
            for(int i = 0; i < spSetup->blockRows(k); i++)
            {
               data(i, spSetup->colIdx(j,k)) = i + 100*j + 100*100*k;
            }
         }
      }

      return spField;
   }
}

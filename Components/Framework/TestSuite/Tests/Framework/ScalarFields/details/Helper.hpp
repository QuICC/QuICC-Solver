#include <memory>
#include "Memory/Cpu/NewDelete.hpp"
#include "Memory/MemoryResource.hpp"
#include "QuICC/ScalarFields/ScalarFieldSetup.hpp"
#include "QuICC/ScalarFields/FlatScalarField.hpp"
#include "QuICC/ScalarFields/ViewScalarField.hpp"

namespace details {

   enum class SetupType {
      UniformUp,
      UniformDown,
      TriangularUp,
      TriangularDown
   };

   std::shared_ptr<QuICC::Datatypes::ScalarFieldSetup> createSetup(const SetupType type, const int dim1D, const int dim3D);

   template <typename T> std::shared_ptr<QuICC::Datatypes::FlatScalarField<T>> createFlatScalarField(const SetupType type, const int dim1D, const int dim3D);

   template <typename T> std::shared_ptr<QuICC::Datatypes::ViewScalarField<T>> createViewScalarField(const SetupType type, const int dim1D, const int dim3D);

   template <typename T> T fieldValueA(const int i, const int j, const int k);
   template <typename T> T fieldValueB(const int i, const int j, const int k);

   template <typename T> T fieldValueA(const int i, const int j, const int k)
   {
      return i + 100*j + 100*100*k;
   }

   template <typename T> T fieldValueB(const int i, const int j, const int k)
   {
      return -700 - k*100*100 - j - i*0.001;
   }
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
               data(i, spSetup->colIdx(j,k)) = fieldValueA<T>(i,j,k);
            }
         }
      }

      return spField;
   }

   template <typename T> std::shared_ptr<QuICC::Datatypes::ViewScalarField<T>> createViewScalarField(const SetupType type, const int dim1D, const int dim3D)
   {
      auto spSetup = details::createSetup(details::SetupType::UniformUp, dim1D, dim3D);

      auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();

      auto spField = std::make_shared<QuICC::Datatypes::ViewScalarField<T>>(spSetup, mem);
      auto&& data = spField->rDataView();

      for(int k = 0; k < spSetup->nBlock(); k++)
      {
         for(int j = 0; j < spSetup->blockCols(k); j++)
         {
            for(int i = 0; i < spSetup->blockRows(k); i++)
            {
               data(i, j, k) = fieldValueA<T>(i,j,k);                                                    ;
            }
         }
      }

      return spField;
   }
}

#include "Helper.hpp"

namespace details {

   std::shared_ptr<QuICC::Datatypes::ScalarFieldSetup> createSetup(const SetupType type, const int dim1D, const int dim3D)
   {
      auto spDim1D = std::make_shared<QuICC::ArrayI>(dim3D);
      auto spDim2D = std::make_shared<QuICC::ArrayI>(dim3D);

      if(type == SetupType::UniformUp)
      {
         spDim1D->setConstant(dim1D);
         for(int i = 0; i < dim3D; i++)
         {
            (*spDim2D)(i) = i + 2;
         }
      }
      else if(type == SetupType::UniformDown)
      {
         spDim1D->setConstant(dim1D);
         for(int i = 0; i < dim3D; i++)
         {
            (*spDim2D)(i) = dim1D - i;
         }
      }
      else if(type == SetupType::TriangularUp)
      {
         for(int i = 0; i < dim3D; i++)
         {
            (*spDim1D)(i) = i + dim1D;
            (*spDim2D)(i) = i + 1;
         }
      }
      else if(type == SetupType::TriangularDown)
      {
         for(int i = 0; i < dim3D; i++)
         {
            (*spDim1D)(i) = dim1D - i;
            (*spDim2D)(i) = dim1D - 2*i;
         }
      }
      else
      {
         throw std::logic_error("Unknown setup type");
      }

      auto spSetup = std::make_shared<QuICC::Datatypes::ScalarFieldSetup>(spDim1D, spDim2D, dim3D);

      return spSetup;
   }
}


/**
 * @file CorrectSolutionFunctorS.cpp
 * @brief Implementation of the CorrectSolution functor fpr SINGLE
 */

// System includes
//

// Project includes
//
#include "QuICC/Equations/details/CorrectSolutionFunctorS.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <>
void CorrectSolutionFunctor<CouplingIndexType::SINGLE>::init(const std::vector<std::tuple<MHDVariant, int, int, int>>& corrections)
{
   if(corrections.size() > 0)
   {
      assert(matIdx == 0);

      const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
      const auto& sRes = eq->res().sim();

      // Copy data
      int l, k_, j_, dimK, dimJ;

      switch(sRes.ss().dimension())
      {
         case 3:
            dimK = sRes.dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL)*sRes.dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
            dimJ = sRes.dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
            break;
         case 2:
            dimK = 1;
            dimJ = sRes.dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
            break;
         case 1:
            dimK = 1;
            dimJ = 1;
            break;
         default:
            dimK = -1;
            dimJ = -1;
            throw  std::logic_error("Spatial scheme has unknown dimension!");
      }

      for (auto&& c: corrections)
      {
         int k = std::get<3>(c);
         int j = std::get<2>(c);
         int i = std::get<1>(c);
         k_ = tRes.template idx<Dimensions::Data::DAT3D>(k)*dimK;
         j_ = tRes.template idx<Dimensions::Data::DAT2D>(j,k)*dimJ;
         l = k_ + j_ + i;
         corr.emplace_back(std::get<MHDComplex>(std::get<0>(c)), l, 0);
      }
   }
}

} // namespace details
} // namespace Equations
} // namespace QuICC

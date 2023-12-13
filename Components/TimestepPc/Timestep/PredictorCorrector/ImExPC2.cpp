/**
 * @file ImExPC2.cpp
 * @brief Implementation of an implicit/explicit predictor-corrector scheme of
 * order 2
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "Timestep/PredictorCorrector/ImExPC2.hpp"

namespace QuICC {

namespace Timestep {

namespace PredictorCorrector {

// Scheme requires 3 substeps
int ImExPC2::steps() const
{
   return 2;
}

// Scheme order
int ImExPC2::order() const
{
   return 2;
}

// Scheme has embedded lower order scheme?
bool ImExPC2::hasEmbedded() const
{
   return true;
}

// Name of the scheme
std::string ImExPC2::name() const
{
   return "ImExPC2";
}

ImExPC2::ImExPC2() : IImExPCScheme(), mAIm({}), mCEx({})
{
   this->init();
}

MHDFloat ImExPC2::aIm(const int i) const
{
   return this->mAIm[i];
}

MHDFloat ImExPC2::cEx(const int i) const
{
   return this->mCEx[i];
}

void ImExPC2::init()
{
   // Initialize implicit a factors
   std::fill(this->mAIm.begin(), this->mAIm.end(), 1. / 2.);

   // Initialize step fractions
   std::fill(this->mCEx.begin(), this->mCEx.end(), 1);
}

} // namespace PredictorCorrector
} // namespace Timestep
} // namespace QuICC

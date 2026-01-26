/**
 * @file ImExEuler.cpp
 * @brief Implementation of an implicit/explicit predictor-corrector scheme of
 * order 2
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "Timestep/PredictorCorrector/ImExEuler.hpp"

namespace QuICC {

namespace Timestep {

namespace PredictorCorrector {

// Scheme requires 1 substep
int ImExEuler::steps() const
{
   return 1;
}

// Scheme order
int ImExEuler::order() const
{
   return 1;
}

// Scheme has embedded lower order scheme?
bool ImExEuler::hasEmbedded() const
{
   return true;
}

// Name of the scheme
std::string ImExEuler::name() const
{
   return "ImExEuler";
}

ImExEuler::ImExEuler() : IImExPCScheme(), mAIm({}), mCEx({})
{
   this->init();
}

MHDFloat ImExEuler::aIm(const int i) const
{
   return this->mAIm[i];
}

MHDFloat ImExEuler::cEx(const int i) const
{
   return this->mCEx[i];
}

void ImExEuler::init()
{
   // Initialize implicit a factors
   std::fill(this->mAIm.begin(), this->mAIm.end(), 1.);

   // Initialize step fractions
   std::fill(this->mCEx.begin(), this->mCEx.end(), 1);
}

} // namespace PredictorCorrector
} // namespace Timestep
} // namespace QuICC

/**
 * @file ImExRKCB3e.cpp
 * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order 3e
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "Timestep/RungeKuttaCB/ImExRKCB3e.hpp"

namespace QuICC {

namespace Timestep {

namespace RungeKuttaCB {

// Scheme requires 3 substeps
int ImExRKCB3e::steps() const
{
   return 4;
}

// Scheme order
int ImExRKCB3e::order() const
{
   return 3;
}

// Scheme has embedded lower order scheme?
bool ImExRKCB3e::hasEmbedded() const
{
   return false;
}

// Name of the scheme
std::string ImExRKCB3e::name() const
{
   return "ImExRKCB3e";
}

ImExRKCB3e::ImExRKCB3e() :
    IImExRKCBScheme(),
    mAIm(Eigen::Array<MHDFloat, 4, 4>::Zero()),
    mBIm(Eigen::Array<MHDFloat, 4, 1>::Zero()),
    mAEx(Eigen::Array<MHDFloat, 4, 4>::Zero()),
    mBEx(Eigen::Array<MHDFloat, 4, 1>::Zero()),
    mCEx(Eigen::Array<MHDFloat, 4, 1>::Zero()),
    mBImErr(Eigen::Array<MHDFloat, 4, 1>::Zero()),
    mBExErr(Eigen::Array<MHDFloat, 4, 1>::Zero())
{
   this->init();
}

MHDFloat ImExRKCB3e::aIm(const int i, const int j) const
{
   return ImExRKCB3e::mAIm(i, j);
}

MHDFloat ImExRKCB3e::bIm(const int i) const
{
   return ImExRKCB3e::mBIm(i);
}

MHDFloat ImExRKCB3e::aEx(const int i, const int j) const
{
   return ImExRKCB3e::mAEx(i, j);
}

MHDFloat ImExRKCB3e::bEx(const int i) const
{
   return ImExRKCB3e::mBEx(i);
}

MHDFloat ImExRKCB3e::cEx(const int i) const
{
   return ImExRKCB3e::mCEx(i);
}

MHDFloat ImExRKCB3e::bImErr(const int i) const
{
   return ImExRKCB3e::mBImErr(i);
}

MHDFloat ImExRKCB3e::bExErr(const int i) const
{
   return ImExRKCB3e::mBExErr(i);
}

void ImExRKCB3e::init()
{
   // Initialize implicit a factors
   ImExRKCB3e::mAIm(1, 1) = 1. / 3.;
   ImExRKCB3e::mAIm(2, 1) = 1. / 2.;
   ImExRKCB3e::mAIm(2, 2) = 1. / 2.;
   ImExRKCB3e::mAIm(3, 1) = 3. / 4.;
   ImExRKCB3e::mAIm(3, 2) = -1. / 4.;
   ImExRKCB3e::mAIm(3, 3) = 1. / 2.;

   // Initialize implicit b factors
   ImExRKCB3e::mBIm(1) = 3. / 4.;
   ImExRKCB3e::mBIm(2) = -1. / 4.;
   ImExRKCB3e::mBIm(3) = 1. / 2.;

   // Initialize explicit a factors
   ImExRKCB3e::mAEx(1, 0) = 1. / 3.;
   ImExRKCB3e::mAEx(2, 1) = 1.;
   ImExRKCB3e::mAEx(3, 1) = 3. / 4.;
   ImExRKCB3e::mAEx(3, 2) = 1. / 4.;

   // Initialize explicit b factors
   ImExRKCB3e::mBEx(1) = 3. / 4.;
   ImExRKCB3e::mBEx(2) = -1. / 4.;
   ImExRKCB3e::mBEx(3) = 1. / 2.;

   // Initialize explicit c factors
   ImExRKCB3e::mCEx(1) = 1. / 3.;
   ImExRKCB3e::mCEx(2) = 1.;
   ImExRKCB3e::mCEx(3) = 1.;
}

} // namespace RungeKuttaCB
} // namespace Timestep
} // namespace QuICC

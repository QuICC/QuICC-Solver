/**
 * @file ImExRKCB4.cpp
 * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order 4
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "Timestep/RungeKuttaCB/ImExRKCB4.hpp"

namespace QuICC {

namespace Timestep {

namespace RungeKuttaCB {

// Scheme requires 3 substeps
int ImExRKCB4::steps() const
{
   return 6;
}

// Scheme order
int ImExRKCB4::order() const
{
   return 4;
}

// Scheme has embedded lower order scheme?
bool ImExRKCB4::hasEmbedded() const
{
   return true;
}

// Name of the scheme
std::string ImExRKCB4::name() const
{
   return "ImExRKCB4";
}

ImExRKCB4::ImExRKCB4() :
    IImExRKCBScheme(),
    mAIm(Eigen::Array<MHDFloat, 6, 6>::Zero()),
    mBIm(Eigen::Array<MHDFloat, 6, 1>::Zero()),
    mAEx(Eigen::Array<MHDFloat, 6, 6>::Zero()),
    mBEx(Eigen::Array<MHDFloat, 6, 1>::Zero()),
    mCEx(Eigen::Array<MHDFloat, 6, 1>::Zero()),
    mBImErr(Eigen::Array<MHDFloat, 6, 1>::Zero()),
    mBExErr(Eigen::Array<MHDFloat, 6, 1>::Zero())
{
   this->init();
}

MHDFloat ImExRKCB4::aIm(const int i, const int j) const
{
   return ImExRKCB4::mAIm(i, j);
}

MHDFloat ImExRKCB4::bIm(const int i) const
{
   return ImExRKCB4::mBIm(i);
}

MHDFloat ImExRKCB4::aEx(const int i, const int j) const
{
   return ImExRKCB4::mAEx(i, j);
}

MHDFloat ImExRKCB4::bEx(const int i) const
{
   return ImExRKCB4::mBEx(i);
}

MHDFloat ImExRKCB4::cEx(const int i) const
{
   return ImExRKCB4::mCEx(i);
}

MHDFloat ImExRKCB4::bImErr(const int i) const
{
   return ImExRKCB4::mBImErr(i);
}

MHDFloat ImExRKCB4::bExErr(const int i) const
{
   return ImExRKCB4::mBExErr(i);
}

void ImExRKCB4::init()
{
   MHDFloat b1 = 232049084587. / 1377130630063.;
   MHDFloat b2 = 322009889509. / 2243393849156.;
   MHDFloat b3 = -195109672787. / 1233165545817.;
   MHDFloat b4 = -340582416761. / 705418832319.;
   MHDFloat b5 = 463396075661. / 409972144477.;
   MHDFloat b6 = 323177943294. / 1626646580633.;

   MHDFloat c2 = 1. / 4.;
   MHDFloat c3 = 3. / 4.;
   MHDFloat c4 = 3. / 8.;
   MHDFloat c5 = 1. / 2.;

   // Initialize implicit a factors
   ImExRKCB4::mAIm(1, 0) = c2 / 2.0;
   ImExRKCB4::mAIm(1, 1) = c2 / 2.0;
   ImExRKCB4::mAIm(2, 0) = 216145252607. / 961230882893.;
   ImExRKCB4::mAIm(2, 1) = 257479850128. / 1143310606989.;
   ImExRKCB4::mAIm(2, 2) = 30481561667. / 101628412017.;
   ImExRKCB4::mAIm(3, 0) = b1;
   ImExRKCB4::mAIm(3, 1) = -381180097479. / 1276440792700.;
   ImExRKCB4::mAIm(3, 2) = -54660926949. / 461115766612.;
   ImExRKCB4::mAIm(3, 3) = 344309628413. / 552073727558.;
   ImExRKCB4::mAIm(4, 0) = b1;
   ImExRKCB4::mAIm(4, 1) = b2;
   ImExRKCB4::mAIm(4, 2) = -100836174740. / 861952129159.;
   ImExRKCB4::mAIm(4, 3) = -250423827953. / 1283875864443.;
   ImExRKCB4::mAIm(4, 4) = 1. / 2.;
   ImExRKCB4::mAIm(5, 0) = b1;
   ImExRKCB4::mAIm(5, 1) = b2;
   ImExRKCB4::mAIm(5, 2) = b3;
   ImExRKCB4::mAIm(5, 3) = b4;
   ImExRKCB4::mAIm(5, 4) = b5;
   ImExRKCB4::mAIm(5, 5) = b6;

   // Initialize implicit b factors
   ImExRKCB4::mBIm(0) = b1;
   ImExRKCB4::mBIm(1) = b2;
   ImExRKCB4::mBIm(2) = b3;
   ImExRKCB4::mBIm(3) = b4;
   ImExRKCB4::mBIm(4) = b5;
   ImExRKCB4::mBIm(5) = b6;

   // Initialize explicit a factors
   ImExRKCB4::mAEx(1, 0) = c2;
   ImExRKCB4::mAEx(2, 0) = 153985248130. / 1004999853329.;
   ImExRKCB4::mAEx(2, 1) = 902825336800. / 1512825644809.;
   ImExRKCB4::mAEx(3, 0) = b1;
   ImExRKCB4::mAEx(3, 1) = 99316866929. / 820744730663.;
   ImExRKCB4::mAEx(3, 2) = 82888780751. / 969573940619.;
   ImExRKCB4::mAEx(4, 0) = b1;
   ImExRKCB4::mAEx(4, 1) = b2;
   ImExRKCB4::mAEx(4, 2) = 57501241309. / 765040883867.;
   ImExRKCB4::mAEx(4, 3) = 76345938311. / 676824576433.;
   ImExRKCB4::mAEx(5, 0) = b1;
   ImExRKCB4::mAEx(5, 1) = b2;
   ImExRKCB4::mAEx(5, 2) = b3;
   ImExRKCB4::mAEx(5, 3) = -4099309936455. / 6310162971841.;
   ImExRKCB4::mAEx(5, 4) = 1395992540491. / 933264948679.;

   // Initialize explicit b factors
   ImExRKCB4::mBEx(0) = b1;
   ImExRKCB4::mBEx(1) = b2;
   ImExRKCB4::mBEx(2) = b3;
   ImExRKCB4::mBEx(3) = b4;
   ImExRKCB4::mBEx(4) = b5;
   ImExRKCB4::mBEx(5) = b6;

   // Initialize explicit c factors
   ImExRKCB4::mCEx(1) = c2;
   ImExRKCB4::mCEx(2) = c3;
   ImExRKCB4::mCEx(3) = c4;
   ImExRKCB4::mCEx(4) = c5;
   ImExRKCB4::mCEx(5) = 1.;

   MHDFloat b1e = 5590918588. / 49191225249.;
   MHDFloat b2e = 92380217342. / 122399335103.;
   MHDFloat b3e = -29257529014. / 55608238079.;
   MHDFloat b4e = -126677396901. / 66917692409.;
   MHDFloat b5e = 384446411890. / 169364936833.;
   MHDFloat b6e = 58325237543. / 207682037557.;

   // Initialize embedded implicit b factors
   ImExRKCB4::mBImErr(0) = b1e;
   ImExRKCB4::mBImErr(1) = b2e;
   ImExRKCB4::mBImErr(2) = b3e;
   ImExRKCB4::mBImErr(3) = b4e;
   ImExRKCB4::mBImErr(4) = b5e;
   ImExRKCB4::mBImErr(5) = b6e;

   // Initialize embedded explicit b factors
   ImExRKCB4::mBExErr(0) = b1e;
   ImExRKCB4::mBExErr(1) = b2e;
   ImExRKCB4::mBExErr(2) = b3e;
   ImExRKCB4::mBExErr(3) = b4e;
   ImExRKCB4::mBExErr(4) = b5e;
   ImExRKCB4::mBExErr(5) = b6e;
}

} // namespace RungeKuttaCB
} // namespace Timestep
} // namespace QuICC

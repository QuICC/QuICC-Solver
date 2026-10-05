/**
 * @file EnergyR2viewCpu_t.hpp.inc
 * @brief Wrapper of the Worland EnergyR2 Reductor
 */
#pragma once


// System includes
//
#include <cstdio>
#include <vector>
#include <complex>


// Project includes
//
#include "Types/Typedefs.hpp"
#include "View/View.hpp"
#include "Operator/Unary.hpp"
#include "Operator/Nary.hpp"
#include "Memory/MemoryResource.hpp"
#include "Memory/Memory.hpp"
#include "Common/include/QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Resolution.hpp"

namespace QuICC {

namespace Timestep  {

namespace Exponential {

namespace Functors {

/**
* @brief Implementation of the Worland based EnergyR2 Reductor
*/
class EnergyFunctor
{
public:
    /**
     * @brief Constructor
     */
    EnergyFunctor(std::shared_ptr<Resolution> spRes, const int gSize, std::shared_ptr<Memory::memory_resource> mem);

    /**
     * @brief Destructor
     */
    ~EnergyFunctor() = default;

    using mods_t = View::View<std::complex<double>, View::DCCSC3D>;

    /// @brief Compute inner product (u, u)
    /// @param in
    /// @para id  Field ID
    double operator()(mods_t in, const SpectralFieldId& id) const;

    /// @brief Compute inner product (u,v)
    /// @param u
    /// @para id  Field ID
    double operator()(mods_t u, mods_t v, const SpectralFieldId& id) const;

    /// @brief Convert to energy basis
    void convert(mods_t out, mods_t in, const SpectralFieldId& id) const;

    /// @brief Initialize operators
    /// @param gSize Grid size
    void init(mods_t in) const;

    /// @brief Energy basis View
    mutable mods_t mEModsView;

private:
    using power_t = View::View<double, View::CSC>;
    template <typename TScaleL, typename TScaleM>
       double computeEnergy(power_t in, const std::unique_ptr<TScaleL>& pScaleL, const std::unique_ptr<TScaleM>& pScaleM) const;

    std::shared_ptr<Resolution> spRes;

    /// @brief Grid size
    const int mcGridSize;

   /**
    * @brief
    */
   std::shared_ptr<Memory::memory_resource> _mem;

    using phys_t = View::View<std::complex<double>, View::DCCSC3D>;
    using modsAbs_t = View::View<double, View::DCCSC3D>;
    using op_t = View::View<double, View::CSL3DJIK>;

    /// @brief Projector backend class pointer
    mutable std::unique_ptr<QuICC::Operator::UnaryOp<phys_t, mods_t>> mPrjR2;

    /// @brief Projector backend class pointer
    mutable std::unique_ptr<QuICC::Operator::UnaryOp<phys_t, mods_t>> mPrjOverR1;

    /// @brief Projector backend class pointer
    mutable std::unique_ptr<QuICC::Operator::UnaryOp<phys_t, mods_t>> mPrjOverR1D1R1;

    /// @brief Integrator backend class pointer
    mutable std::unique_ptr<QuICC::Operator::UnaryOp<mods_t, phys_t>> mIntR2;

    /// @brief Integrator backend class pointer
    mutable std::unique_ptr<QuICC::Operator::UnaryOp<mods_t, phys_t>> mIntOverR1;

    /// @brief Integrator backend class pointer
    mutable std::unique_ptr<QuICC::Operator::UnaryOp<mods_t, phys_t>> mIntOverR1D1R1;

    /// @brief Pointwise backend class pointer
    mutable std::unique_ptr<QuICC::Operator::NaryOp<modsAbs_t, mods_t, mods_t>> mDot;

    /// @brief Reduction backend class pointer
    mutable std::unique_ptr<QuICC::Operator::UnaryOp<power_t, modsAbs_t>> mRed;

    /// @brief Input View
    mutable mods_t mModsView;

    /// @brief temp/output View
    mutable modsAbs_t mModsAbsView;

    /// @brief temp/output View
    mutable power_t mModsRedView;

    /// @brief temp/output View
    mutable phys_t mPhysView;

    /// @brief Temporary storage for flattened input
    mutable Memory::MemBlock<std::complex<double>> mModsData;

    /// @brief Temporary storage for flattened input
    mutable Memory::MemBlock<std::complex<double>> mEModsData;

    /// @brief Temporary storage for flattened temp/output
    mutable Memory::MemBlock<double> mModsAbsData;

    /// @brief Temporary storage for flattened temp/output
    mutable Memory::MemBlock<double> mModsRedData;

    /// @brief Temporary storage for flattened temp/output
    mutable Memory::MemBlock<std::complex<double>> mPhysData;
};

///@brief Scale for m = 0
struct ScaleM
{
   ScaleM() = default;
   ~ScaleM() = default;
   double operator()(const std::uint32_t l, const std::uint32_t m) const;
};

///@brief No scaling
struct NoScale
{
   NoScale() = default;
   ~NoScale() = default;
   template <typename ...Args>
   double operator()(const Args ...) const{return 1.0;};
};

///@brief Energy scaling for q component
struct ScaleLq
{
   ScaleLq() = default;
   ~ScaleLq() = default;
   double operator()(const std::uint32_t l) const;
};

///@brief Energy scaling for s and t components
struct ScaleLst
{
   ScaleLst() = default;
   ~ScaleLst() = default;
   double operator()(const std::uint32_t l) const;
};

template <typename TScaleL, typename TScaleM> double EnergyFunctor::computeEnergy(power_t in, const std::unique_ptr<TScaleL>& pScaleL, const std::unique_ptr<TScaleM>& pScaleM) const
{
   const auto& scaleL = *pScaleL;
   const auto& scaleM = *pScaleM;

   double energy = 0.0;
   double lfactor = 1.0;
   double factor = 1.0;
   std::uint32_t P = in.dims()[1];
   for(std::uint32_t l = 0; l < P; l++)
   {
      lfactor = scaleL(l);
      for(std::uint32_t j = 0; j < in.pointers()[0][l+1] - in.pointers()[0][l]; j++)
      {
         auto i = in.pointers()[0][l] + j;
         auto m = in.indices()[0][i];

         factor = scaleM(l, m) * lfactor;

         energy += factor*mModsRedView(m, l);
      }
   }

   return energy;
}

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

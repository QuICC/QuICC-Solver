/**
 * @file Utils.cpp
 * @brief Pseudospectral utils
 */

// System includes
//
#include <stdexcept>
#include <type_traits>
#include <boost/functional/hash.hpp>

// Project includes
//
#include "QuICC/Pseudospectral/Utils.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"

namespace QuICC {

namespace Pseudospectral {

std::size_t hash_combine(const std::size_t a, const std::size_t b)
{
   std::size_t seed = 0;
   boost::hash_combine(seed, a);
   boost::hash_combine(seed, b);
   return seed;
}

namespace details
{
   ptrAndIdxBlock getMeta(const TransformResolution& res, const std::uint32_t maxLayers, std::shared_ptr<Memory::memory_resource> mem)
   {
      std ::uint32_t nLayers = res.dim<QuICC::Dimensions::Data::DAT3D>();

      // ptr
      Memory::MemBlock<std::uint32_t> ptrBlock(maxLayers+1, mem.get());
      View::ViewBase<std::uint32_t> ptr(ptrBlock.data(), ptrBlock.size());

      using namespace QuICC::Memory;
      tempOnHostMemorySpace ConverterP(ptr, TransferMode::write | TransferMode::block);

      std::uint32_t cumLayerSize = 0;
      std::uint32_t layerCounter = 0;
      ptr[0] = 0;
      for(std::uint32_t l = 0; l < maxLayers; ++l)
      {
         std::uint32_t layerSize = 0;
         if (layerCounter < nLayers)
         {
            auto layerIndex = static_cast<std::uint32_t>(res.idx<QuICC::Dimensions::Data::DAT3D>(layerCounter));
            if (l == layerIndex)
            {
               layerSize = res.dim<QuICC::Dimensions::Data::DAT2D>(layerCounter);
               ++layerCounter;
            }
         }
         ptr[l+1] = ptr[l] + layerSize;
         cumLayerSize += layerSize;
      }

      // idx
      Memory::MemBlock<std::uint32_t> idxBlock(cumLayerSize, mem.get());
      View::ViewBase<std::uint32_t> idx(idxBlock.data(), idxBlock.size());

      tempOnHostMemorySpace ConverterI(idx, TransferMode::write);

      std::uint32_t l = 0;
      for(std::uint32_t k = 0; k < nLayers; ++k)
      {
         for(int j = 0; j < res.dim<QuICC::Dimensions::Data::DAT2D>(k); ++j)
         {
            // auto layerIndex = res.idx<QuICC::Dimensions::Data::DAT3D>(k);
            auto columnIndex = res.idx<QuICC::Dimensions::Data::DAT2D>(j, k);
            // auto columnHeight = res.dim<QuICC::Dimensions::Data::DATB1D>(j, k);
            idx[l] = columnIndex;
            ++l;
         }
      }

      ptrAndIdxBlock ret;
      ret.ptr = std::move(ptrBlock);
      ret.idx = std::move(idxBlock);
      return ret;
   }

   template <class SCALAROUT, class SCALARIN, class VIEWATT>
   void copyEig2View(QuICC::View::View<SCALAROUT, VIEWATT> view, const Eigen::Matrix<SCALARIN, -1, -1>& eig, const TransformResolution& res)
   {
      throw std::logic_error("trying to copy different types");
   }

   template <class SCALAR, class VIEWATT>
   void copyEig2View(QuICC::View::View<SCALAR, VIEWATT> view, const Eigen::Matrix<SCALAR, -1, -1>& eig, const TransformResolution& res)
   {
      std ::uint32_t nLayers = res.dim<QuICC::Dimensions::Data::DAT3D>();

      auto pointers = view.pointers()[1];

      using namespace QuICC::Memory;

      /// \todo permanent host buffer
      View::ViewBase viewb(view.data(), view.size());
      tempOnHostMemorySpace Converter(viewb, TransferMode::write | TransferMode::block);

      /// \todo permanent transfer / check width only in debug mode
      tempOnHostMemorySpace ConverterP(pointers, TransferMode::read | TransferMode::block);


      /// copy data to view
      int start = 0;
      std::int64_t offSet = 0;
      for(std::uint32_t p = 0; p < nLayers; ++p)
      {
         // layer width
         int cols = res.dim<QuICC::Dimensions::Data::DAT2D>(p);
         auto widthView = static_cast<std::uint32_t>(cols);
         // check only in debug mode
         // auto widthView = pointers[p+1] - pointers[p];
         // assert(widthView == static_cast<std::uint32_t>(cols));
         // layer height (DCCSC3D only)
         /// \todo extend to other view types
         int inRows;
         if constexpr (std::is_same_v<SCALAR, double>)
         {
            // phys
            inRows = res.dim<QuICC::Dimensions::Data::DATF1D>(0, p);
         }
         else
         {
            // mods
            inRows = res.dim<QuICC::Dimensions::Data::DATB1D>(0, p);
         }

         auto heightView = view.lds();
         if constexpr (std::is_same_v<VIEWATT, View::DCCSC3D>)
         {
            assert(heightView == static_cast<std::uint32_t>(inRows));
         }

         const Eigen::Ref<const Eigen::Matrix<SCALAR, -1, -1>> inB = eig.block(0, start, inRows, cols);

         for (std::int64_t j = 0; j < inB.cols(); ++j)
         {
            for (std::int64_t i = 0; i < inB.rows(); ++i)
            {
               if constexpr (std::is_same_v<VIEWATT, View::DCCSC3DJIK>)
               {
                  // copy padded to flattened and transpose
                  viewb[offSet + i*widthView+j] = inB.data()[i+j*eig.rows()];
               }
               else if constexpr (std::is_same_v<VIEWATT, View::DCCSC3D>)
               {
                  // copy padded to flattened column
                  viewb[offSet + i+j*heightView] = inB.data()[i+j*eig.rows()];
               }
               else
               {
                  throw std::logic_error("copy 2 view not implemented for this layout");
               }
            }
         }
         offSet += inB.size();
         start += cols;
      }
   }

   template <class SCALAROUT, class SCALARIN, class VIEWATT>
   void copyView2Eig(Eigen::Matrix<SCALAROUT, -1, -1>& eig, const QuICC::View::View<SCALARIN, VIEWATT> view, const TransformResolution& res)
   {
      throw std::logic_error("trying to copy different types");
   }

   template <class SCALAR, class VIEWATT>
   void copyView2Eig(Eigen::Matrix<SCALAR, -1, -1>& eig, const QuICC::View::View<SCALAR, VIEWATT> view, const TransformResolution& res)
   {
      std ::uint32_t nLayers = res.dim<QuICC::Dimensions::Data::DAT3D>();

      auto pointers = view.pointers()[1];

      using namespace QuICC::Memory;
      /// \todo permanent host buffer
      View::ViewBase viewb(view.data(), view.size());
      tempOnHostMemorySpace Converter(viewb, TransferMode::read);

      /// \todo permanent transfer / check width only in debug mode
      tempOnHostMemorySpace ConverterP(pointers, TransferMode::read | TransferMode::block);

      int start = 0;
      std::uint64_t offSet = 0;
      for(std::uint32_t p = 0; p < nLayers; p++)
      {
         // layer width
         int cols = res.dim<QuICC::Dimensions::Data::DAT2D>(p);
         auto widthView = static_cast<std::uint32_t>(cols);
         // check only in debug mode
         // auto widthView = pointers[p+1] - pointers[p];
         // assert(widthView == static_cast<std::uint32_t>(cols));
         // layer height (DCCSC3D only)
         /// \todo extend to other view types
         int outRows;
         if constexpr (std::is_same_v<SCALAR, double>)
         {
            // phys
            outRows = res.dim<QuICC::Dimensions::Data::DATF1D>(0, p);
         }
         else
         {
            // mods
            outRows = res.dim<QuICC::Dimensions::Data::DATB1D>(0, p);
         }

         auto heightView = view.lds();
         if constexpr (std::is_same_v<VIEWATT, View::DCCSC3D>)
         {
            assert(heightView == static_cast<std::uint32_t>(outRows));
         }

         Eigen::Ref<Eigen::Matrix<SCALAR, -1, -1>> outB = eig.block(0, start, outRows, cols);

         for (std::int64_t j = 0; j < outB.cols(); ++j)
         {
            for (std::int64_t i = 0; i < outB.rows(); ++i)
            {
               if constexpr (std::is_same_v<VIEWATT, View::DCCSC3DJIK>)
               {
                  // copy padded from flattened column and transpose
                  outB.data()[i+j*eig.rows()] = viewb[offSet + i*widthView+j];
               }
               else if constexpr (std::is_same_v<VIEWATT, View::DCCSC3D>)
               {
                  // copy padded from flattened column
                  outB.data()[i+j*eig.rows()] = viewb[offSet + i+j*heightView];
               }
               else
               {
                  throw std::logic_error("copy from view not implemented for this layout");
               }
           }
         }

         offSet += outB.size();
         start += cols;
      }
   }

   void copyScalar2View(Graph::varData_t vVar, Framework::Selector::VariantSharedScalarVariable sVar, const TransformResolution& res)
   {
      std::visit(
         [&](auto& Tv, auto&& p)
         {
            auto& ptrTemp = p->rDom(0).rPerturbation();
            details::copyEig2View(Tv, ptrTemp.data(), res);
         }, vVar, sVar);
   }

   void copyVector2View(Graph::varData_t vVar0, Graph::varData_t vVar1, Framework::Selector::VariantSharedVectorVariable vecVar, const TransformResolution& res)
   {
      std::visit(
         [&](auto& Torv, auto& Polv, auto&& p)
         {
            auto& ptrTor = p->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR);
            details::copyEig2View(Torv, ptrTor.data(), res);
            auto& ptrPol = p->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::POL);
            details::copyEig2View(Polv, ptrPol.data(), res);
         }, vVar0, vVar1, vecVar);
   }

   void copyView2Scalar(Framework::Selector::VariantSharedScalarVariable sVar, Graph::varData_t vVar, const TransformResolution& res)
   {
      std::visit(
         [&](auto&& p, auto& Tv)
         {
            auto& ptrTemp = p->rDom(0).rPerturbation();
            p->rDom(0).rPerturbation().setZeros();
            details::copyView2Eig(ptrTemp.rData(), Tv, res);
         }, sVar, vVar);
   }

   void copyView2Vector(Framework::Selector::VariantSharedVectorVariable vecVar, Graph::varData_t vVar0, Graph::varData_t vVar1, const TransformResolution& res)
   {
      std::visit(
         [&](auto&& p, auto& Torv, auto& Polv)
         {
            auto& ptrTor = p->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR);
            p->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::TOR).setZeros();
            details::copyView2Eig(ptrTor.rData(), Torv, res);
            auto& ptrPol = p->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::POL);
            p->rDom(0).rPerturbation().rComp(FieldComponents::Spectral::POL).setZeros();
            details::copyView2Eig(ptrPol.rData(), Polv, res);
         }, vecVar, vVar0, vVar1);

   }

   void copyView2Vector(Framework::Selector::VariantSharedVectorVariable vecVar, Graph::varData_t vVar0, Graph::varData_t vVar1, Graph::varData_t vVar2, const TransformResolution& res)
   {
      std::visit(
            [&](auto&& p, auto& Urv, auto& Uthetav, auto& Uphiv)
            {
               auto& ptrUr = p->rDom(0).rPhys().rComp(FieldComponents::Physical::R);
               p->rDom(0).rPhys().rComp(FieldComponents::Physical::R).setZeros();
               details::copyView2Eig(ptrUr.rData(), Urv, res);
               auto& ptrUtheta = p->rDom(0).rPhys().rComp(FieldComponents::Physical::THETA);
               p->rDom(0).rPhys().rComp(FieldComponents::Physical::THETA).setZeros();
               details::copyView2Eig(ptrUtheta.rData(), Uthetav, res);
               auto& ptrUphi = p->rDom(0).rPhys().rComp(FieldComponents::Physical::PHI);
               p->rDom(0).rPhys().rComp(FieldComponents::Physical::PHI).setZeros();
               details::copyView2Eig(ptrUphi.rData(), Uphiv, res);
            }, vecVar, vVar0, vVar1, vVar2);

   }


} // namespace details
} // Pseudospectral
} // QuICC

/**
 * @file Utils.hpp
 * @brief Pseudospectral Utils
 */

#ifndef QUICC_PSEUDOSPECTRAL_UTILS_HPP
#define QUICC_PSEUDOSPECTRAL_UTILS_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "View/View.hpp"

namespace QuICC {

namespace Pseudospectral {

namespace details
{
   std::size_t hash_combine(const std::size_t a, const std::size_t b);


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

} // namespace details
} // Pseudospectral
} // QuICC

#endif // QUICC_PSEUDOSPECTRAL_UTILS_HPP

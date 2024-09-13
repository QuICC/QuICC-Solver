/**
 * @file Attributes.hpp
 * @brief Attributes/traits for wrapper for the sparse solver
 */

#pragma once

// System includes
//
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <variant>

// Project includes
//

namespace QuICC {
namespace View {

   //
   // Level types tags
   //

   /// @brief tag type for dense level
   /// a dense level has implicit coordinates
   /// zeros are stored explicitly
   struct dense_t {};

   /// @brief tag type for compressed level
   /// a compressed level has explicit coordinates
   struct compressed_t {};

   /// @brief tag type for sparse level
   /// a sparse level has implicit coordinates
   /// zeros are not stored explicitly
   /// for instance <sparse_t, compressed_t> -> CSR
   /// this tag needs always be associated with a compressed_t tag
   struct sparse_t {};

   /// @brief tag type for dense level with implicit starts from k
   /// to the full dimension, so that the length is implicity equal
   /// to the full dimension minus K (k being the 3rd index)
   /// \todo template on tag index and coeff?
   struct step1K_t {};

   /// @brief tag type for dense level with implicit start starts from 2*k
   /// to the full dimension, so that the length is implicity equal
   /// to the full dimension minus 2*K (k being the 3rd index)
   /// \todo template on tag index and coeff?
   struct step2K_t {};

   /// @brief check if a type is a valid level type attribute
   /// @tparam T type to be checked
   template <class T>
   struct isLevelType : std::false_type {};

   /// @brief enable dense_t tag
   template <>
   struct isLevelType<dense_t> : std::true_type{};

   /// @brief enable compressed_t tag
   template <>
   struct isLevelType<compressed_t> : std::true_type{};

   /// @brief enable sparse_t tag
   template <>
   struct isLevelType<sparse_t> : std::true_type{};

   /// @brief enable step1K_t tag
   template <>
   struct isLevelType<step1K_t> : std::true_type{};

   /// @brief enable step2K_t tag
   template <>
   struct isLevelType<step2K_t> : std::true_type{};

   /// @brief helper
   /// @tparam T level tag
   template <class T>
   inline constexpr bool isLevelType_v = isLevelType<T>::value;

   /// @brief check all levels
   /// @tparam ...Ts
   template<class ... Ts>
   struct areLevelType {
      static constexpr bool value {(isLevelType_v<Ts> && ...)};
   };

   /// @brief helper
   /// @tparam Ts
   template <class ... Ts>
   inline constexpr bool areLevelType_v = areLevelType<Ts...>::value;


   /// @brief generic template to check if a level is dense
   /// @tparam T
   template <class T>
   struct isLevelTypeDense : std::false_type
   {
      static_assert(isLevelType_v<T>, "type is not a level type");
   };

   /// @brief specialiazed template to check if a level is dense (for dense_t)
   template <>
   struct isLevelTypeDense<dense_t> : std::true_type{};

   /// @brief helper
   /// @tparam T
   template <class T>
   inline constexpr bool isLevelTypeDense_v = isLevelTypeDense<T>::value;

   /// @brief check if of level are dense, one by one
   /// @tparam ...Ts
   template<class... Ts>
   struct areLevelTypeDense {
      static constexpr bool value {(isLevelTypeDense_v<Ts> && ...)};
   };

   /// @brief helper
   /// @tparam T
   template <class... Ts>
   inline constexpr bool areLevelTypeDense_v = areLevelTypeDense<Ts...>::value;

   /// @brief count how many dense levels
   /// @tparam ...Ts
   template<class... Ts>
   struct howManyLevelTypeDense {
      static constexpr std::uint16_t value {(static_cast<uint16_t>(isLevelTypeDense_v<Ts>) + ...)};
   };

   /// @brief helper
   /// @tparam T
   template <class... Ts>
   inline constexpr std::uint16_t howManyLevelTypeDense_v = howManyLevelTypeDense<Ts...>::value;

   namespace details
   {
      /// @brief check with human readable error message
      /// @tparam ...T
      /// @return
      template <class... T>
      constexpr bool assertLevel()
      {
         static_assert(areLevelType<T...>::value, "unknown level type, must be either dense_t, compressed_t, triDense_t, triCompressed_t or CSC_t");
         return true;
      }
   }

   /// @brief aggregate dimension levels type, i.e. contained in a std::variant
   /// @tparam ...Types
   template <class... Types>
   using DimLevelType = std::enable_if_t<details::assertLevel<Types...>(), std::variant<Types...>>;

   /// @brief generic template to check if aggregate level types is the correct type
   /// @tparam Tvar
   template<class T>
   struct isDimLevelType : std::false_type {};

   /// @brief specialized template to check if aggregate level types is the correct type
   /// @tparam ...Ts
   template<class... Ts>
   struct isDimLevelType<std::variant<Ts...>> : std::true_type
   {
      static_assert(details::assertLevel<Ts...>());
   };

   /// @brief helper
   /// @tparam T
   template <class T>
   inline constexpr bool isDimLevelType_v = isDimLevelType<T>::value;

   /// @brief generic template to check if aggregate level types (aka DimLevelType) are all dense
   /// @tparam Tvar
   template<class T>
   struct isDimLevelTypeDense : std::false_type
   {
      static_assert(isDimLevelType_v<T>, "Type not a DimLevel.");
   };

   /// @brief specialized template to check if aggregate level types (aka DimLevelType) are all dense_t
   /// @tparam ...Ts
   template<class... Ts>
   struct isDimLevelTypeDense<std::variant<Ts...>> : public areLevelTypeDense<Ts...> {};

   /// @brief helper
   /// @tparam T
   template <class T>
   inline constexpr bool isDimLevelTypeDense_v = isDimLevelTypeDense<T>::value;


   /// @brief generic template to check if aggregate level types (aka DimLevelType) are all dense
   /// @tparam Tvar
   template<class T>
   struct howManyDimLevelTypeDense
   {
      static_assert(isDimLevelType_v<T>, "Type not a DimLevel.");
      static constexpr std::uint16_t value = 0;
   };

   /// @brief specialized template to check if aggregate level types (aka DimLevelType) are all dense_t
   /// @tparam ...Ts
   template<class... Ts>
   struct howManyDimLevelTypeDense<std::variant<Ts...>> : public howManyLevelTypeDense<Ts...> {};

   /// @brief helper
   /// @tparam T
   template <class T>
   inline constexpr std::uint16_t howManyDimLevelTypeDense_v = howManyDimLevelTypeDense<T>::value;


   //
   // Loop attribute tags (for now only 3D)
   //

   /// @brief tag type for logical first index
   struct i_t {};
   /// @brief tag type for logical first index
   /// with in-order mapping (complex DFT) [-N...N] -> [0...N,-N...-1]
   /// which requires special mapping for padding
   struct iInOrder_t {};
   /// @brief tag type for logical second index
   struct j_t {};
   /// @brief tag type for logical third index
   struct k_t {};
   /// @brief tag type for logical fourth index
   struct l_t {};

   /// @brief check if a type is a valid loop type attribute tag
   /// @tparam T
   template <class T>
   struct isLoopType : std::false_type {};

   /// @brief enable i_t
   template <>
   struct isLoopType<i_t> : std::true_type{};

   /// @brief enable iInOrder_t
   template <>
   struct isLoopType<iInOrder_t> : std::true_type{};

   /// @brief enable j_t
   template <>
   struct isLoopType<j_t> : std::true_type{};

   /// @brief enable k_t
   template <>
   struct isLoopType<k_t> : std::true_type{};

   /// @brief enable l_t
   template <>
   struct isLoopType<l_t> : std::true_type{};

   /// @brief helper
   /// @tparam T
   template <class T>
   inline constexpr bool isLoopType_v = isLoopType<T>::value;

   /// @brief check that all tags are valid loop tags
   /// @tparam ...T
   template<class ... T>
   struct areLoopType {
      static constexpr bool value {(isLoopType_v<T> && ...)};
   };

   /// @brief helper
   /// @tparam T
   template <class T>
   inline constexpr bool areLoopType_v = areLoopType<T>::value;

   namespace details
   {
      /// @brief check with human readable error message
      /// @tparam ...T
      /// @return
      template <class... T>
      constexpr bool assertOrder()
      {
         // add unique tag check
         static_assert(areLoopType<T...>::value, "unknown level type, must be i_t, j_t, k_t or l_t");
         return true;
      }
   }

   /// @brief aggregate loop order type, i.e. contained in a std::variant
   /// @tparam ...Types
   template <class... Types>
   using LoopOrderType = std::enable_if_t<details::assertOrder<Types...>(), std::variant<Types...>>;


   /// @brief generic template for default loop order, i.e. column major (leftmost index has stride 1)
   /// @tparam I
   template <std::uint32_t I>
   struct DefaultLoopOrder;

   /// @brief specialized 1D default loop order
   template <>
   struct DefaultLoopOrder<1>
   {
      using type = LoopOrderType<i_t>;
   };

   /// @brief specialized 2D default loop order
   template <>
   struct DefaultLoopOrder<2>
   {
      using type = LoopOrderType<i_t, j_t>;
   };

   /// @brief specialized 3D default loop order
   template <>
   struct DefaultLoopOrder<3>
   {
      using type = LoopOrderType<i_t, j_t, k_t>;
   };

      /// @brief specialized 4D default loop order
   template <>
   struct DefaultLoopOrder<4>
   {
      using type = LoopOrderType<i_t, j_t, k_t, l_t>;
   };

   /// @brief helper
   template <std::uint32_t I>
   using DefaultLoopOrder_t = typename DefaultLoopOrder<I>::type;


   /// Attributes wrapper

   /** @brief attributes wrapper.
    *  At the moment only Level and Order attributes are implemented
    *  as modifiable attributes. It can be extend as needed.
    * @tparam ...
    */
   template<class...>
   struct Attributes {};

   /** @brief attributes wrapper to store Level and Order
    *  \todo SFINAE std::enable_if_t<areLevelType_v<Level> && areLoopType_v<Order>> = true
    *  @tparam Level
    *  @tparam Order
    *  @tparam ...Rest
    */
   template<class Level, class Order, class... Rest>
   struct Attributes<Level, Order, Rest...>
   {
      /// @brief store aggregate level tags
      using level = Level;
      /// @brief store aggregate order tags
      using order = Order;
      /// @brief default index type
      using index = std::uint32_t;
   };

   /** @brief attributes wrapper to store Level and Order
    *  \todo SFINAE std::enable_if_t<areLevelType_v<Level>> = true
    *  @tparam Level
    *  @tparam ...Rest
    */
   template<class Level, class... Rest>
   struct Attributes<Level, Rest...>
   {
      /// @brief store aggregate level tags
      using level = Level;
      /// @brief default aggregate order tags
      using order = DefaultLoopOrder_t<std::variant_size_v<Level>>;
      /// @brief default index type
      using index = std::uint32_t;
   };


   /** @brief check if type has level attribute.
    *  Generic template handles types that have no nested ::level type member
    *  @tparam type to check
    *  @tparam SFINAE param
    */
   template <class, class = std::void_t<>> struct hasLevel : std::false_type
   {
   };

   /** @brief check if type has level attribute.
    *  Specialization recognizes types that do have a nested ::level type member
    *  @tparam T type to check
    */
   template <class T>
   struct hasLevel<T, std::void_t<typename T::level>> : std::true_type
   {
   };

   /** @brief check if type has order attribute.
    *  Generic template handles types that have no nested ::order type member
    *  @tparam type to check
    *  @tparam SFINAE param
    */
   template <class, class = std::void_t<>> struct hasOrder : std::false_type
   {
   };

   /** @brief check if type has order attribute.
    *  Specialization recognizes types that do have a nested ::order type member
    *  @tparam T type to check
    */
   template <class T>
   struct hasOrder<T, std::void_t<typename T::order>> : std::true_type
   {
   };

   //
   // Common formats
   //

   /**
    * @brief 1D dense tensor
    */
   using dense1D = Attributes<DimLevelType<dense_t>>;

   /**
    * @brief 2D dense tensor with column major memory layout
    */
   using dense2D = Attributes<DimLevelType<dense_t, dense_t>>;

   /**
    * @brief 2D dense tensor with row major memory layout
    */
   using dense2DRM = Attributes<DimLevelType<dense_t, dense_t>, LoopOrderType<j_t, i_t>>;

   /**
    * @brief 2D CSC column major matrix (M,N)
    */
   using CSC = Attributes<DimLevelType<compressed_t, sparse_t>>;

   /**
    * @brief 2D CSC column major matrix (N,K) which elements are 1D dense vectors,
    * i.e a 3D tensor (M,N,K) with fully populated columns
    *
    * used as
    * AL projector output type on cpu
    * FFT in/out types on cpu / gpu
    *
    */
   using DCCSC3D = Attributes<DimLevelType<dense_t, compressed_t, sparse_t>>;

   /**
    * @brief 2D CSC column major matrix (N,K) which elements are 1D dense vectors,
    * i.e a 3D tensor (M,N,K) with fully populated columns
    *
    * used as
    * AL projector output type on cpu
    * FFT in/out types on cpu / gpu
    *
    */
   using DCCSC3DInOrder = Attributes<DimLevelType<dense_t, compressed_t, sparse_t>,
      LoopOrderType<iInOrder_t, j_t, k_t>>;


   /**
    * @brief 2D CSC column major matrix (N,K) which elements are 1D dense vectors,
    * i.e a 3D tensor (M,N,K) with fully populated columns
    *
    * used as
    * AL projector output type on gpu
    *
    */
   using DCCSC3DJIK = Attributes<DimLevelType<dense_t, compressed_t, sparse_t>,
      LoopOrderType<j_t, i_t, k_t>>;


   /**
    * @brief JW projector/integrator operator type on cpu
    *  with uniform truncation
    *  1D compressed array with 2D dense elements
    *  i.e a 3D tensor (M,N,K) with fully populated layers
    */
   using CSL3DJIK = Attributes<DimLevelType<dense_t, dense_t, compressed_t>,
      LoopOrderType<j_t, i_t, k_t>>;

   /**
    * @brief JW projector/integrator operator type on gpu
    *  with uniform truncation
    *  1D compressed array with 2D dense elements
    *  i.e a 3D tensor (M,N,K) with fully populated layers
    */
   using CSL3D = Attributes<DimLevelType<dense_t, dense_t, compressed_t>>;

   /**
    * @brief AL projector operator type on cpu
    *  Compressed step 1 row layer
    */
   using CS1RL3DJIK = Attributes<DimLevelType<dense_t, step1K_t, compressed_t>,
      LoopOrderType<j_t, i_t, k_t>>;

   /**
    * @brief AL projector operator type on gpu
    *  Compressed step 1 row layer
    */
   using CS1RL3D = Attributes<DimLevelType<dense_t, step1K_t, compressed_t>>;

   /**
    * @brief AL projector input type on cpu
    * AL integrator op
    * Step 1 column layer, with (N,K) plane a 2D CSC column major matrix
    */
   using S1CLCSC3D = Attributes<DimLevelType<step1K_t, compressed_t, sparse_t>>;

   /**
    * @brief AL projector input type on gpu
    * AL integrator op
    * Step 1 column layer, with (N,K) plane a 2D CSC column major matrix
    */
   using S1CLCSC3DJIK = Attributes<DimLevelType<step1K_t, compressed_t, sparse_t>,
      LoopOrderType<j_t, i_t, k_t>>;

   /**
    * @brief AL projector operator type on cpu
    *  Compressed step 2 row layer
    */
   using CS2RL3DJIK = Attributes<DimLevelType<dense_t, step2K_t, compressed_t>,
      LoopOrderType<j_t, i_t, k_t>>;

   /**
    * @brief AL projector operator type on gpu
    *  Compressed step 2 row layer
    */
   using CS2RL3D = Attributes<DimLevelType<dense_t, step2K_t, compressed_t>>;

   /**
    * @brief JW projector input type on cpu
    * Step 2 column layer, with (N,K) plane a 2D CSC column major matrix
    */
   using S2CLCSC3D = Attributes<DimLevelType<step2K_t, compressed_t, sparse_t>>;

   /**
    * @brief JW projector input type on gpu
    * Step 2 column layer, with (N,K) plane a 2D CSC column major matrix
    */
   using S2CLCSC3DJIK = Attributes<DimLevelType<step2K_t, compressed_t, sparse_t>,
      LoopOrderType<j_t, i_t, k_t>>;

   /**
    * @brief AL projector operator type on cpu
    *  Compressed step 2 row layer
    */
   using CS2RL3DJIK = Attributes<DimLevelType<dense_t, step2K_t, compressed_t>,
      LoopOrderType<j_t, i_t, k_t>>;

   /**
    * @brief AL projector operator type on gpu
    *  Compressed step 2 row layer
    */
   using CS2RL3D = Attributes<DimLevelType<dense_t, step2K_t, compressed_t>>;

   /**
    * @brief JW projector input type on cpu
    * Step 2 column layer, with (N,K) plane a 2D CSC column major matrix
    */
   using S2CLCSC3D = Attributes<DimLevelType<step2K_t, compressed_t, sparse_t>>;

   /**
    * @brief JW projector input type on gpu
    * Step 2 column layer, with (N,K) plane a 2D CSC column major matrix
    */
   using S2CLCSC3DJIK = Attributes<DimLevelType<step2K_t, compressed_t, sparse_t>,
      LoopOrderType<j_t, i_t, k_t>>;


} // namespace View
} // namespave QuICC

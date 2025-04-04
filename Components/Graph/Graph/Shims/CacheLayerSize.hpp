/**
 * @file CacheLayerSize.hpp
 * @brief Cache layer size
 */
#pragma once

// External includes
//

#include <map>

// Project includes
//


namespace QuICC
{
namespace Graph
{

namespace details
{

/// @brief CacheLayerSize class
/// Computing the size on the device and copying to the host in oreder
/// to allocate the temporary is expensive.
/// We cache the size of the layers in a map based on the metadata pointer.
/// @tparam T type of ptr to cache
template <class T>
class CacheLayerSize
{
public:
    /// @brief delete copy ctor
    /// @param
    CacheLayerSize(CacheLayerSize const&) = delete;
    /// @brief delete assignment op
    /// @param
    void operator=(CacheLayerSize const&) = delete;
    /// @brief create/return singleton
    /// @return singleton
    static CacheLayerSize& getInstance();
    /// @brief access the memory resource
    /// @return
    std::map<T*, std::size_t>& getMap();
    /// @brief invalidate the entire cache
    void invalidate();
    /// @brief invalidate a specific entry in the cache
    /// @param ptr pointer to the metadata to remove
    void invalidate(const T* ptr);
private:
    /// @brief map pointe to size
    std::map<T*, std::size_t> _ptr2size;
    /// @brief Default ctor
    CacheLayerSize() = default;
    /// @brief dtor
    ~CacheLayerSize() = default;
};

template <class T>
CacheLayerSize<T>& CacheLayerSize<T>::getInstance()
{
    static CacheLayerSize instance;
    return instance;
}

template <class T>
std::map<T*, std::size_t>& CacheLayerSize<T>::getMap()
{
    return _ptr2size;
}

template <class T>
void CacheLayerSize<T>::invalidate()
{
    _ptr2size.clear();
}

template <class T>
void CacheLayerSize<T>::invalidate(const T* ptr)
{
    _ptr2size.erase(ptr);
}

} // namespace details
} // namespace Graph
} // namespace QuICC
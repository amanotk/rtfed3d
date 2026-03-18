// -*- C++ -*-
#ifndef _ARRAY_TYPES_HPP_
#define _ARRAY_TYPES_HPP_

#include <algorithm>
#include <array>
#include <cstddef>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>

#include <experimental/mdspan>

template <class T, std::size_t Rank>
class Array
{
public:
  using T_numtype    = T;
  using extents_type = std::dextents<int, Rank>;
  using mdspan_type  = std::mdspan<T, extents_type, std::layout_right>;

  Array() : extents_{}, view_()
  {
  }

  Array(const Array& other) : extents_(other.extents_), storage_(other.storage_), view_()
  {
    update_view(std::make_index_sequence<Rank>{});
  }

  Array(Array&& other) noexcept
      : extents_(other.extents_), storage_(std::move(other.storage_)), view_()
  {
    update_view(std::make_index_sequence<Rank>{});
  }

  Array& operator=(const Array& other)
  {
    if (this != &other) {
      extents_ = other.extents_;
      storage_ = other.storage_;
      update_view(std::make_index_sequence<Rank>{});
    }
    return *this;
  }

  Array& operator=(Array&& other) noexcept
  {
    if (this != &other) {
      extents_ = other.extents_;
      storage_ = std::move(other.storage_);
      update_view(std::make_index_sequence<Rank>{});
    }
    return *this;
  }

  template <class... Dims, typename std::enable_if<sizeof...(Dims) == Rank, int>::type = 0>
  void resize(Dims... dims)
  {
    extents_ = {static_cast<int>(dims)...};

    std::size_t count = 1;
    for (std::size_t r = 0; r < Rank; r++) {
      count *= static_cast<std::size_t>(extents_[r]);
    }

    storage_.resize(count);
    update_view(std::make_index_sequence<Rank>{});
  }

  template <class... Indices, typename std::enable_if<sizeof...(Indices) == Rank, int>::type = 0>
  T& operator()(Indices... indices)
  {
    return view_(static_cast<int>(indices)...);
  }

  template <class... Indices, typename std::enable_if<sizeof...(Indices) == Rank, int>::type = 0>
  const T& operator()(Indices... indices) const
  {
    return view_(static_cast<int>(indices)...);
  }

  Array& operator=(const T& value)
  {
    std::fill(storage_.begin(), storage_.end(), value);
    return *this;
  }

  T* data()
  {
    return storage_.data();
  }

  const T* data() const
  {
    return storage_.data();
  }

  int rank() const
  {
    return static_cast<int>(Rank);
  }

  int extent(int r) const
  {
    return extents_[r];
  }

  std::size_t size() const
  {
    return storage_.size();
  }

  std::array<int, Rank> lbound() const
  {
    std::array<int, Rank> lower = {};
    return lower;
  }

  std::array<int, Rank> ubound() const
  {
    std::array<int, Rank> upper = {};

    for (std::size_t r = 0; r < Rank; r++) {
      upper[r] = extents_[r] - 1;
    }

    return upper;
  }

private:
  template <std::size_t... Indices>
  void update_view(std::index_sequence<Indices...>)
  {
    view_ = mdspan_type(storage_.data(), extents_[Indices]...);
  }

  std::array<int, Rank> extents_;
  std::vector<T>        storage_;
  mdspan_type           view_;
};

#endif

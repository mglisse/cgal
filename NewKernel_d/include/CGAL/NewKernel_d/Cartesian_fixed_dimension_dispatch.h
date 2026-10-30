// Copyright (c) 2026
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Marc Glisse

#ifndef CGAL_NEWKERNEL_D_CARTESIAN_FIXED_DIMENSION_DISPATCH_H
#define CGAL_NEWKERNEL_D_CARTESIAN_FIXED_DIMENSION_DISPATCH_H

#include <CGAL/basic.h>
#include <CGAL/NewKernel_d/functor_tags.h>
#include <CGAL/NewKernel_d/store_kernel.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/tools.h>
#include <CGAL/Default.h>
#include <CGAL/Dimension.h>
#include <array>
#include <optional>
#include <tuple>
#include <type_traits>

namespace CGAL {
template <class Dim> struct Epick_d;

namespace internal {
// Keep the dispatch range and predicate list in sync with Cartesian_static_filters.
#ifdef CGAL_NO_STATIC_FILTER_456
constexpr int fixed_epick_max_dimension = 3;
#else
constexpr int fixed_epick_max_dimension = 6;
#endif

// FIXME: move this to Cartesian_static_filters.h ?
template <class Tag, int D>
struct Has_fixed_epick_static_filter
    : std::bool_constant<(std::is_same_v<Tag, Orientation_of_points_tag> && D >= 2 && D <= fixed_epick_max_dimension) ||
                         (std::is_same_v<Tag, Side_of_oriented_sphere_tag> && (D == 2 || D == 3))>
{};

// Access policy keeps approximate-coordinate reads out of the lazy kernel.
struct Fixed_epick_point_access
{
  template <class K> static const K& kernel(const K& k) { return k; }
  template <class P> static const P& point(const P& p) { return p; }
};

template <class R_, class Base_predicate, class Predicate_tag, class Point_access>
struct Fixed_dimension_epick_predicate : private Store_kernel<R_>
{
  CGAL_FUNCTOR_INIT_STORE(Fixed_dimension_epick_predicate)
  using result_type = typename Base_predicate::result_type;

  // FIXME: Orientation_of_points also has an overload taking points directly in fixed dim (yeah, that was a bad idea).
  template <class Iter, class... Args> result_type operator()(Iter f, Iter e, const Args&... args) const {
    const auto& k = this->kernel();
    typename Get_functor<R_, Point_dimension_tag>::type pd(k);
    using Dimension = typename R_::Default_ambient_dimension;
    if constexpr(std::is_same_v<Dimension, Dynamic_dimension_tag>) {
      if(auto result = dispatch<2>(pd(*f), k, f, e, args...))
        return *result;
    } else {
      constexpr int D = Dimension::value;
      if constexpr(Has_fixed_epick_static_filter<Predicate_tag, D>::value)
        if(auto result = evaluate<D>(k, f, e, args...))
          return *result;
    }
    return Base_predicate(k)(f, e, args...);
  }

private:
  template <int D, class Point>
  static bool convert_point(const R_& k, const Point& source, typename Epick_d<Dimension_tag<D>>::Point_d& target) {
    const auto& source_kernel = Point_access::kernel(k);
    const auto& source_point = Point_access::point(source);
    using Source_kernel = std::decay_t<decltype(source_kernel)>;
    typename Get_functor<Source_kernel, Point_dimension_tag>::type pd(source_kernel);
    typename Get_functor<Source_kernel, Compute_point_cartesian_coordinate_tag>::type coordinate(source_kernel);
    if(pd(source_point) != D)
      return false;
    std::array<double, D> coords;
    for(int j = 0; j < D; ++j)
      if(!fit_in_double(coordinate(source_point, j), coords[j]))
        return false;
    target = typename Epick_d<Dimension_tag<D>>::Point_d(D, coords.begin(), coords.end());
    return true;
  }

  template <int D, class Iter, class... Args>
  static std::optional<result_type> evaluate(const R_& k, Iter f, Iter e, const Args&... args) {
    using Fixed_kernel = Epick_d<Dimension_tag<D>>;
    using Fixed_point = typename Fixed_kernel::Point_d;
    // Both predicates take a simplex; sphere additionally takes a query point.
    std::array<Fixed_point, D + 1> points;
    std::size_t i = 0;
    for(; f != e && i < points.size(); ++f, ++i)
      if(!convert_point<D>(k, *f, points[i]))
        return std::nullopt;
    if(f != e || i != points.size())
      return std::nullopt;

    // FIXME: having 1 path for orientation and 1 for sphere would be more readable.
    std::array<Fixed_point, sizeof...(Args)> extra;
    i = 0;
    if(!(convert_point<D>(k, args, extra[i++]) && ...))
      return std::nullopt;
    Fixed_kernel fixed(D);
    typename Get_functor<Fixed_kernel, Predicate_tag>::type predicate(fixed);
    return std::apply([&](const auto&... converted) { return predicate(points.begin(), points.end(), converted...); },
                      extra);
  }

  template <int D, class Iter, class... Args>
  static std::optional<result_type> dispatch(int dimension, const R_& k, Iter f, Iter e, const Args&... args) {
    if constexpr(Has_fixed_epick_static_filter<Predicate_tag, D>::value)
      if(dimension == D)
        return evaluate<D>(k, f, e, args...);
    if constexpr(D < fixed_epick_max_dimension)
      return dispatch<D + 1>(dimension, k, f, e, args...);
    return std::nullopt;
  }
};
} // namespace internal

template <class Base_, class Derived_ = Default, class Point_access = internal::Fixed_epick_point_access>
struct Cartesian_fixed_dimension_dispatch : public Base_
{
  constexpr Cartesian_fixed_dimension_dispatch() {}
  constexpr Cartesian_fixed_dimension_dispatch(int d)
      : Base_(d) {}
  using Self = Cartesian_fixed_dimension_dispatch<Base_, Derived_, Point_access>;
  using Derived = typename Default::Get<Derived_, Self>::type;
  template <class T, class = void> struct Functor : Inherit_functor<Base_, T>
  {};
  template <class D> struct Functor<Orientation_of_points_tag, D>
  {
    using Base_predicate = typename Get_functor<Base_, Orientation_of_points_tag>::type;
    using type =
        internal::Fixed_dimension_epick_predicate<Derived, Base_predicate, Orientation_of_points_tag, Point_access>;
  };
  template <class D> struct Functor<Side_of_oriented_sphere_tag, D>
  {
    using Base_predicate = typename Get_functor<Base_, Side_of_oriented_sphere_tag>::type;
    using type =
        internal::Fixed_dimension_epick_predicate<Derived, Base_predicate, Side_of_oriented_sphere_tag, Point_access>;
  };
};
} // namespace CGAL
#endif // CGAL_NEWKERNEL_D_CARTESIAN_FIXED_DIMENSION_DISPATCH_H

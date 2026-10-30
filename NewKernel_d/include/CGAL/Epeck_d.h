// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Marc Glisse

#ifndef CGAL_EPECK_D_H
#define CGAL_EPECK_D_H

#include <CGAL/disable_warnings.h>

#include <CGAL/NewKernel_d/Cartesian_base.h>
#include <CGAL/NewKernel_d/Wrapper/Cartesian_wrap.h>
#include <CGAL/NewKernel_d/Kernel_d_interface.h>
#include <CGAL/NewKernel_d/Lazy_cartesian.h>
#include <CGAL/NewKernel_d/KernelD_converter.h>
#include <CGAL/Number_types/internal/Exact_type_selector.h>
#include <CGAL/NewKernel_d/Types/Weighted_point.h>
#include <CGAL/Epick_d.h>

// For predicates where Epick_d has a static filter, we try an exact convertion to Epick_d first.
// Other predicates could also benefit (in particular because of the Mpzf fallback), but we don't handle them yet, it would increase compile time even more, and require the wrapper to be more generic.

namespace CGAL {
#define CGAL_KA Cartesian_base_d<Interval_nt_advanced,Dim>
#define CGAL_KE Cartesian_base_d<internal::Exact_field_selector<double>::Type, Dim>
template<class Dim> using Epeck_d_help1 = Lazy_cartesian<CGAL_KE, CGAL_KA, KernelD_converter<CGAL_KE, CGAL_KA>>;
#undef CGAL_KE
#undef CGAL_KA
#undef CGAL_BASE

namespace internal {
struct Fixed_epick_lazy_point_access
{
  template <class K> static decltype(auto) kernel(const K& k) { return k.approximate_kernel(); }
  template <class P> static decltype(auto) point(const P& p) {
    // Read the point approximation directly: coordinate access through the lazy
    // kernel would construct a lazy coordinate node.
    return CGAL::approx(p);
  }
};
} // namespace internal

// Both fixed and runtime dimensions use the same eligible-predicate shortcut.
#define CGAL_BASE \
  Cartesian_fixed_dimension_dispatch< \
    Epeck_d_help1<Dim>, \
    Default, \
    internal::Fixed_epick_lazy_point_access>
template <class Dim>
struct Epeck_d_help2
: CGAL_BASE
{
  constexpr Epeck_d_help2(){}
  constexpr Epeck_d_help2(int d):CGAL_BASE(d){}
};
#undef CGAL_BASE

#define CGAL_BASE \
  Kernel_d_interface< \
    Cartesian_wrap< \
      Epeck_d_help2<Dim>, \
      Epeck_d<Dim> > >
template<class Dim>
struct Epeck_d
: CGAL_BASE
{
  constexpr Epeck_d(){}
  constexpr Epeck_d(int d):CGAL_BASE(d){}
};
#undef CGAL_BASE
}

#include <CGAL/enable_warnings.h>

#endif

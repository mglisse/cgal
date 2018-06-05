// Copyright (c) 2001-2018  INRIA (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Marc Glisse

#ifndef CGAL_FILTERED_COMPUTE_H
#define CGAL_FILTERED_COMPUTE_H

#include <string>
#include <CGAL/config.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Uncertain.h>
#include <CGAL/Profile_counter.h>
#include <CGAL/NewKernel_d/store_kernel.h>
#include <CGAL/NT_converter.h>
#include <boost/preprocessor.hpp>

namespace CGAL {

// This template class is a wrapper that implements the filtering for any
// computation (dynamic filters with IA).

template <class FT, class EC, class AC, class C2E, class C2A, bool Protection = true>
class Filtered_compute
{
  typedef typename decay<typename AC::result_type>::type  Ares;
  typedef typename decay<typename EC::result_type>::type  Eres;

//TODO: pack (at least use a tuple)
//FIXME: is it better to store those, or just store enough to recreate them
//(i.e. possibly references to the kernels)?
  EC  ec;
  AC  ac;
  C2E c2e;
  C2A c2a;
  NT_converter<Ares,FT> a2b;
  NT_converter<Eres,FT> e2b;

public:

  typedef AC    Approximate_compute;
  typedef EC    Exact_compute;
  typedef C2E   To_exact_converter;
  typedef C2A   To_approximate_converter;

  // AK::FT must be Interval_nt and with EK::FT they must be NT_convertible to K::FT
  // It would be nicer if the precision test was embedded in to_double(Interval_nt).
  typedef FT result_type;

  Filtered_compute()
  {}

  template <class K>
  Filtered_compute(const K& k)
    : ec(k.exact_kernel()), ac(k.approximate_kernel()), c2e(k,k.exact_kernel()), c2a(k,k.approximate_kernel())
  {}

#ifdef CGAL_CXX11
  template <typename... Args>
  result_type
  operator()(Args&&... args) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    // Protection is outside the try block as VC8 has the CGAL_CFG_FPU_ROUNDING_MODE_UNWINDING_VC_BUG
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  // No forward here, the arguments may still be needed
	  Ares res = ac(c2a(args)...);
	  if (has_smaller_relative_precision(res, .00001))
	    return a2b(res);
	}
      catch (Uncertain_conversion_exception&) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return e2b(ec(c2e(std::forward<Args>(args))...));
  }
#else

#define CGAL_VAR(Z,N,C) C(a##N)
#define CGAL_CODE(Z,N,_) \
  template <BOOST_PP_ENUM_PARAMS(N,class A)> \
  result_type \
  operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, A, const& a)) const \
  { \
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp); \
    { \
      Protect_FPU_rounding<Protection> p; \
      try \
	{ \
	  Ares res = ac(BOOST_PP_ENUM(N,CGAL_VAR,c2a)); \
	  if (has_smaller_relative_precision(res, .00001)) \
	    return a2b(res); \
	} \
      catch (Uncertain_conversion_exception&) {} \
    } \
    CGAL_BRANCH_PROFILER_BRANCH(tmp); \
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST); \
    return e2b(ec(BOOST_PP_ENUM(N,CGAL_VAR,c2e))); \
  }
  BOOST_PP_REPEAT_FROM_TO(1, 10, CGAL_CODE, _ )
#undef CGAL_CODE
#undef CGAL_VAR

#endif
};

} //namespace CGAL

#endif // CGAL_FILTERED_COMPUTE_H

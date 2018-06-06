// Copyright (c) 2018
// INRIA (France)
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
// Author(s)     : Marc Glisse

#ifndef CGAL_KERNEL_D_CARTESIAN_FILTER_ALL_H
#define CGAL_KERNEL_D_CARTESIAN_FILTER_ALL_H

#include <CGAL/basic.h>
#include <CGAL/NewKernel_d/KernelD_converter.h>
#include <CGAL/NewKernel_d/Filtered_predicate2.h>
#include <CGAL/NewKernel_d/Filtered_compute.h>
#include <CGAL/NewKernel_d/Filtered_construct.h>
#include <boost/mpl/if.hpp>
#include <boost/mpl/and.hpp>

namespace CGAL {

// Duplicates the functionality of Cartesian_filter_K so we don't have to chain them. It would be nicer to make them orthogonal, but possibly not so useful?
template < typename Base_, typename AK_, typename EK_ >
struct Cartesian_filter_all : public Base_,
  private Store_kernel<AK_>, private Store_kernel2<EK_>
{
    CGAL_CONSTEXPR Cartesian_filter_all(){}
    CGAL_CONSTEXPR Cartesian_filter_all(int d):Base_(d){}
    //FIXME: or do we want an instance of AK and EK belonging to this kernel,
    //instead of a reference to external ones?
    CGAL_CONSTEXPR Cartesian_filter_all(AK_ const&a,EK_ const&b):Base_(),Store_kernel<AK_>(a),Store_kernel2<EK_>(b){}
    CGAL_CONSTEXPR Cartesian_filter_all(int d,AK_ const&a,EK_ const&b):Base_(d),Store_kernel<AK_>(a),Store_kernel2<EK_>(b){}
    typedef Base_ Kernel_base;
    typedef AK_ AK;
    typedef EK_ EK;
    typedef typename Store_kernel<AK_>::reference_type AK_rt;
    AK_rt approximate_kernel()const{return this->kernel();}
    typedef typename Store_kernel2<EK_>::reference2_type EK_rt;
    EK_rt exact_kernel()const{return this->kernel2();}

    // MSVC is too dumb to perform the empty base optimization.
    typedef boost::mpl::and_<
      internal::Do_not_store_kernel<Kernel_base>,
      internal::Do_not_store_kernel<AK>,
      internal::Do_not_store_kernel<EK> > Do_not_store_kernel;

    //TODO: C2A/C2E could be able to convert *this into this->kernel() or this->kernel2().
    typedef KernelD_converter<Kernel_base,AK> C2A;
    typedef KernelD_converter<Kernel_base,EK> C2E;
    typedef KernelD_converter<AK,Kernel_base> A2B;
    typedef KernelD_converter<EK,Kernel_base> E2B;

    // fix the types
    // TODO: only fix some types, based on some criterion?
    template<class T> struct Type : Get_type<Kernel_base,T> {};

    template<class T,class D=void,class=typename Get_functor_category<Cartesian_filter_all,T>::type> struct Functor :
      Inherit_functor<Kernel_base,T,D> {};
    template<class T,class D> struct Functor<T,D,Predicate_tag> {
      typedef typename Get_functor<AK, T>::type AP;
      typedef typename Get_functor<EK, T>::type EP;
      typedef Filtered_predicate2<EP,AP,C2E,C2A> type;
    };
    template<class T,class D> struct Functor<T,D,Compute_tag> {
      typedef typename Get_functor<AK, T>::type AP;
      typedef typename Get_functor<EK, T>::type EP;
      typedef typename Get_type<Kernel_base, FT_tag>::type FT;
      // This could use Filtered_construct instead.
      typedef Filtered_compute<FT,EP,AP,C2E,C2A> type;
    };
    template<class T,class D> struct Functor<T,D,Construct_tag> {
      typedef typename Get_functor<AK, T>::type AP;
      typedef typename Get_functor<EK, T>::type EP;
      typedef typename map_result_tag<T>::type result_tag;
      typedef typename Get_type<Cartesian_filter_all,result_tag>::type Res;
      typedef Filtered_construct<Res,EP,AP,C2E,C2A,E2B,A2B> type;
    };
    // FIXME: This uses sqrt, disable for now.
    // We should also disable it when the double computation is exact (e.g. Compute_point_coordinate for a usual kernel)
    template<class D> struct Functor<Point_of_sphere_tag,D,Construct_tag> :
      Inherit_functor<Kernel_base,Point_of_sphere_tag,D> {};
    //// FIXME: what about Misc_tag?
// TODO:
//    template<class T> struct Functor<T,No_filter_tag,Predicate_tag> :
//	    Kernel_base::template Functor<T,No_filter_tag> {};
// TODO:
// detect when Less_cartesian_coordinate doesn't need filtering
};

} //namespace CGAL

#endif // CGAL_KERNEL_D_CARTESIAN_FILTER_ALL_H

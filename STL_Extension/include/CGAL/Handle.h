// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// Author(s)     : Sylvain Pion

#ifndef CGAL_HANDLE_H
#define CGAL_HANDLE_H

#include <cstddef>
#include <CGAL/Handle_for.h>
#include <CGAL/assertions.h>

namespace CGAL {

class Rep
{
    friend class Handle;
  protected:
             Rep() { count = 1; }
    virtual ~Rep() {}

    int      count;
};

class Handle
{
  public:
    
    typedef std::ptrdiff_t Id_type ;
    
    Handle()
	: PTR(static_cast<Rep*>(0)) {}

    Handle(const Handle& x)
    {
      CGAL_precondition( x.PTR != static_cast<Rep*>(0) );
      PTR = x.PTR;
      PTR->count++;
    }

    ~Handle()
    {
      unlink();
    }

    Handle&
    operator=(const Handle& x)
    {
      CGAL_precondition( x.PTR != static_cast<Rep*>(0) );
      x.PTR->count++;
      unlink();
      PTR = x.PTR;
      return *this;
    }

#ifdef CGAL_CXX11
    Handle(Handle&& x)
    {
      CGAL_precondition( x.PTR != static_cast<Rep*>(0) );
      PTR = x.PTR;
      x.PTR = nullptr;
    }

    Handle&
    operator=(Handle&& x)
    {
      CGAL_precondition( x.PTR != static_cast<Rep*>(0) );
#if 1
      std::swap(PTR, x.PTR);
#else
      // Alternative to kill the old value faster
      CGAL_precondition( this != &x );
      unlink();
      PTR = x.PTR;
      x.PTR = nullptr;
#endif
      return *this;
    }
#endif

    int
    refs()  const { return PTR->count; }

    Id_type id() const { return PTR - static_cast<Rep*>(0); }

    bool identical(const Handle& h) const { return PTR == h.PTR; }

  protected:
    Rep* PTR;

    void unlink()
    {
	if ( PTR && (--PTR->count == 0))
	    delete PTR;
    }
};

//inline Handle::Id_type id(const Handle& x) { return x.id() ; }

inline bool identical(const Handle &h1, const Handle &h2) { return h1.identical(h2); }

} //namespace CGAL

#endif // CGAL_HANDLE_H

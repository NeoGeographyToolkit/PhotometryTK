// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __PHOTK_BBOX_ITERATOR_H__
#define __PHOTK_BBOX_ITERATOR_H__

#include <vw/Math/BBox.h>
#include <vw/Math/Vector.h>
#include <boost/iterator/iterator_facade.hpp>

namespace photk {

  template <class BBoxT, class RealT, size_t DimN>
  class BBoxIterator : public boost::iterator_facade<BBoxIterator<BBoxT,RealT,DimN>,
                                                     vw::Vector<RealT,DimN>,
                                                     boost::forward_traversal_tag> {
    friend class boost::iterator_core_access;

    vw::BBox<RealT,DimN> m_bbox;
    mutable vw::Vector<RealT,DimN> m_location;

    void increment() {
      m_location[0]++;
      for ( size_t i = 0; i < DimN - 1; i++ ) {
        if ( m_location[i] == m_bbox.max()[i] ) {
          m_location[i+1]++;
          m_location[i] = m_bbox.min()[i];
        }
      }
    }

    bool equal(BBoxIterator const& other) const {
      return  m_location == other.m_location;
    }

    vw::Vector<RealT,DimN>& dereference() const { return m_location; }
  public:
    BBoxIterator() :
      m_bbox(), m_location() {}

    BBoxIterator(vw::math::BBoxBase<BBoxT,RealT,DimN> const& box) :
      m_bbox(box.impl()), m_location(m_bbox.min()) {
    }

    BBoxIterator(vw::math::BBoxBase<BBoxT,RealT,DimN> const& box, vw::Vector<RealT,DimN> const& loc) :
      m_bbox(box.impl()), m_location(loc) {}
  };

  template <class BBoxT, class RealT, size_t DimN>
  class BBoxRange {
    vw::BBox<RealT,DimN> m_bbox;
  public:
    typedef BBoxIterator<BBoxT,RealT,DimN> iterator;
    typedef iterator const_iterator;

    BBoxRange(vw::math::BBoxBase<BBoxT,RealT,DimN> const& box) :
      m_bbox(box.impl()) {}

    const_iterator begin() const { return const_iterator(m_bbox); }
    const_iterator end() const {
      vw::Vector<RealT,DimN> limit = m_bbox.min();
      limit[DimN-1] = m_bbox.max()[DimN-1];
      return const_iterator(m_bbox,limit);
    }
  };

  template <class BBoxT, class RealT, size_t DimN>
  BBoxRange<BBoxT, RealT, DimN>
  bbox_range( vw::math::BBoxBase<BBoxT,RealT,DimN> const& box ) {
    return BBoxRange<BBoxT, RealT, DimN>( box.impl() );
  }
}

#endif//__PHOTK_BBOX_ITERATOR_H__

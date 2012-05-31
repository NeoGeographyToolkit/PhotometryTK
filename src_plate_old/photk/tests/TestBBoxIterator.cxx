// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <gtest/gtest.h>
#include <test/Helpers.h>
#include <boost/foreach.hpp>
#include <photk/BBoxIterator.h>

using namespace vw;
using namespace photk;

TEST(BBoxIterator, BoostForeach2i) {
  size_t element_count = 0;
  BBox2i box(150,-50,250,200);

  BOOST_FOREACH( Vector2i const& v, bbox_range(box) ) {
    EXPECT_TRUE( box.contains(v) ) << v;
    element_count++;
  }
  EXPECT_EQ( box.width()*box.height(), element_count );

  element_count = 0;
  BOOST_FOREACH( Vector2i v, bbox_range(box) ) {
    EXPECT_TRUE( box.contains(v) ) << v;
    element_count++;
  }
  EXPECT_EQ( box.width()*box.height(), element_count );
}

TEST(BBoxIterator, BoostForeach3i) {
  size_t element_count = 0;
  BBox3i box(40,-40,-15,50,50,50);

  BOOST_FOREACH( Vector3i const& v, bbox_range(box) ) {
    EXPECT_TRUE( box.contains(v) ) << v;
    element_count++;
  }
  EXPECT_EQ( box.width()*box.height()*box.depth(), element_count );
}

TEST(BBoxIterator, BoostForeach2u) {
  size_t element_count = 0;
  BBox<uint32,2> box(50,1627,88,99);

  typedef Vector<uint32,2> Vector2u;
  BOOST_FOREACH( Vector2u const& v, bbox_range(box) ) {
    EXPECT_TRUE( box.contains(v) ) << v;
    element_count++;
  }
  EXPECT_EQ( box.width()*box.height(), element_count );
}

TEST(BBoxIterator, VerifyPattern) {
  BBox2i box(150,-50,250,200);
  BBoxRange<BBox2i,int32,2>::iterator iter = bbox_range(box).begin();
  for ( int32 y = box.min().y(); y < box.max().y(); y++ ) {
    for ( int32 x = box.min().x(); x < box.max().x(); x++ ) {
      EXPECT_EQ( Vector2i(x,y), *iter );
      iter++;
    }
  }
}

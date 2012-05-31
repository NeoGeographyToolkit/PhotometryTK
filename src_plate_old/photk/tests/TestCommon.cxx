// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
#include <photk/Common.h>
#include <test/Helpers.h>
#include <boost/foreach.hpp>

using namespace vw;
using namespace photk;

TEST(Common, SplitSquareEqualArea) {
  for ( int32 level = 3; level < 20; level++ ) {
    int32 level_size = 0x1 << level;
    int32 k_divisions = level_size;
    if ( k_divisions > 61 )
      k_divisions = 61;
    for ( int32 k = 1; k < k_divisions; k++ ) {
      std::vector<BBox2i> regions =
        split_square_to_equal_area(level_size, k);
      ASSERT_EQ( k, regions.size() ) << "Level: " << level << " size: " << level_size << " k: " << k;
      int32 total_area = 0;
      BOOST_FOREACH( BBox2i const& box, regions )
        total_area += box.width()*box.height();
      EXPECT_EQ( level_size*level_size, total_area );
    }
  }

  // Verify the split_square_to_equal_area produce bboxes that overlap
  // and are contained by the input.
  std::vector<BBox2i> regions = split_square_to_equal_area(256,80);
  BBox2i global_box(0,0,256,256);
  for ( size_t i = 0; i < 80; i++ ) {
    EXPECT_TRUE( global_box.contains(regions[i]) );
    for ( size_t j = i+1; j < 80; j++ )
      EXPECT_FALSE( regions[i].intersects(regions[j]) );
  }

}

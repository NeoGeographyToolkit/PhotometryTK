// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <gtest/gtest.h>
#include <photk/PlateCommon.h>
#include <test/Helpers.h>

using namespace vw;
using namespace photk;

TEST( PlateCommon, TileCacheUnique ) {
  TileCache cache;
  ASSERT_EQ( 0, cache.size() );

  cache = TileCache(8,256);
  ASSERT_EQ( 8, cache.size() );

  for ( size_t i = 0; i < cache.size(); i++ ) {
    EXPECT_EQ( i, cache.index );
    TileCache::result_type& shallow_copy = cache.get();
    EXPECT_EQ( &cache.m_tiles[i](0,0), &shallow_copy(0,0) );
    EXPECT_EQ( i+1, cache.index );
    for ( size_t j = i+1; j < cache.size(); j++ )
      EXPECT_NE( &cache.m_tiles[i](0,0), &cache.m_tiles[j](0,0) );
  }
}

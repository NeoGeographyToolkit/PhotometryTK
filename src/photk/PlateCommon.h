// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

/// \file PlateCommon.h
///

#ifndef __PHOTK_PLATE_COMMON_H__
#define __PHOTK_PLATE_COMMON_H__

#include <map>
#include <list>
#include <vw/Image/ImageView.h>
#include <vw/Plate/PlateFile.h>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

namespace photk {

  typedef boost::shared_ptr<vw::platefile::PlateFile> PlatePtr;
  typedef std::list<vw::platefile::TileHeader> HeaderList;
  typedef boost::tuple<vw::uint32,vw::uint32> rowcol_t;

  typedef std::map<rowcol_t,HeaderList> composite_map_t;
  typedef std::pair<vw::uint32,vw::ImageView<vw::PixelGrayA<vw::float32> > > trans_view_t;
  typedef std::map<rowcol_t,std::list<trans_view_t> > cache_map_t;

  struct TileCache {
    typedef vw::ImageView<vw::PixelGrayA<vw::float32> > result_type;
    std::vector<result_type > m_tiles;
    size_t index;
    TileCache( size_t n, size_t tile_size );

    TileCache() : index(0) {}

    result_type& get();

    size_t size() const { return m_tiles.size(); }
  };

  // Helpful functions to deal with the consumption of tiles.
  vw::uint64 calc_cache_tile_count( PlatePtr plate );
  void cache_consume_tiles(PlatePtr plate, HeaderList const& headers,
                           cache_map_t& cmap, TileCache& cache );
  composite_map_t build_map(HeaderList const& headers);
}

#endif//__PHOTK_PLATE_COMMON_H__

// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <photk/PlateCommon.h>
#include <boost/foreach.hpp>

using namespace vw;
using namespace photk;

TileCache::TileCache( size_t n, size_t tile_size ) : m_tiles(n), index(0) {
  BOOST_FOREACH( result_type& tile, m_tiles )
    tile.set_size(tile_size,tile_size);
}

TileCache::result_type& TileCache::get() {
  VW_ASSERT( index < m_tiles.size(), ArgumentErr() << "Requesting tile outside cache, " << index << " of " << m_tiles.size() << "\n" );
  return m_tiles[index++];
}

vw::uint64 photk::calc_cache_tile_count( PlatePtr plate ) {
  const vw::uint64 TILE_BYTES = plate->default_tile_size() * plate->default_tile_size() * uint64(PixelNumBytes<PixelGrayA<float32> >::value);
  const vw::uint64 CACHE_BYTES = vw_settings().system_cache_size();
  return CACHE_BYTES/TILE_BYTES;
}

void photk::cache_consume_tiles(PlatePtr plate, HeaderList const& headers,
                         cache_map_t& cmap, TileCache& cache ) {
  cmap.clear(); // Clear to insure no redundant data. Redundant is
  // possible since we don't hash on transaction ID.
  platefile::Datastore::TileSearch tile_lookup;
  tile_lookup.reserve( headers.size() );
  std::copy(headers.begin(), headers.end(),
            std::back_inserter(tile_lookup));
  BOOST_FOREACH(const platefile::Tile& t, plate->batch_read(tile_lookup)) {
    trans_view_t image_data(t.hdr.transaction_id(), cache.get());
    boost::scoped_ptr<SrcImageResource> r(SrcMemoryImageResource::open(t.hdr.filetype(),&t.data->operator[](0), t.data->size()));
    read_image(image_data.second, *r);
    cmap[rowcol_t(t.hdr.row(),t.hdr.col())].push_back( image_data );
  }
}

composite_map_t photk::build_map(HeaderList const& headers) {
  composite_map_t cmap;
  BOOST_FOREACH( platefile::TileHeader const& hdr, headers )
    cmap[rowcol_t(hdr.row(),hdr.col())].push_back(hdr);
  return cmap;
}


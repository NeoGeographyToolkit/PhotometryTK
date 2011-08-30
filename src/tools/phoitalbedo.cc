// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// What's this file supposed to do ?
//
// (Pho)tometry (It)eration (Albedo) Update
//
// With Reflectance
// .... see docs
//
// With out Relectance
//      A = A + sum((I^k-T^k*A)*T^k*S^k)/sum((T^k*S^k)^2)
#include <photk/BBoxIterator.h>
#include <photk/Common.h>
#include <photk/Macros.h>
#include <photk/RemoteProjectFile.h>
#include <photk/AlbedoAccumulators.h>
#include <vw/Image.h>
#include <vw/Plate/PlateFile.h>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

using namespace vw;
using namespace vw::platefile;

// Custom types to cut down on my keystrokes
typedef boost::shared_ptr<PlateFile> PlatePtr;
typedef std::list<TileHeader> HeaderList;
typedef boost::tuple<uint32,uint32> rowcol_t;

typedef std::map<rowcol_t,HeaderList> composite_map_t;
typedef std::pair<uint32,ImageView<PixelGrayA<float32> > > trans_view_t;
typedef std::map<rowcol_t,std::list<trans_view_t> > cache_map_t;

struct TileCache {
  typedef ImageView<PixelGrayA<float32> > result_type;
  std::vector<result_type > m_tiles;
  size_t index;
  TileCache( size_t n, size_t tile_size ) : m_tiles(n) {
    BOOST_FOREACH( result_type& tile, m_tiles )
      tile.set_size(tile_size,tile_size);
  }

  TileCache() : index(0) {}

  result_type& get() {
    VW_ASSERT( index < m_tiles.size(), ArgumentErr() << "Requesting tile outside cache, " << index << " of " << m_tiles.size() << "\n" );
    return m_tiles[index++];
  }

  size_t size() const { return m_tiles.size(); }
};

struct Options : photk::BaseOptions {
  // Input
  Url ptk_url;
  int32 level;
  bool perform_2band;

  // For spawning multiple jobs
  int32 job_id, num_jobs;
  BBox2i work_area;

  // Commonly reused variables
  boost::scoped_ptr<photk::RemoteProjectFile> remote_ptk;
  photk::ProjectMeta project_info;
  PlatePtr drg_plate, albedo_plate, reflect_plate;
  std::vector<double> exposure_vec;
  size_t tile_size;
  TileCache tile_cache;
};

uint64 calc_cache_tile_count( PlatePtr plate ) {
  const uint64 TILE_BYTES = plate->default_tile_size() * plate->default_tile_size() * uint64(PixelNumBytes<PixelGrayA<float32> >::value);
  const uint64 CACHE_BYTES = vw_settings().system_cache_size();
  return CACHE_BYTES/TILE_BYTES;
}

void cache_consume_tiles(PlatePtr plate, HeaderList const& headers,
                         cache_map_t& cmap, TileCache& cache ) {
  cmap.clear(); // Clear to insure no redundant data. Redundant is
                // possible since we don't hash on transaction ID.
  Datastore::TileSearch tile_lookup;
  tile_lookup.reserve( headers.size() );
  std::copy(headers.begin(), headers.end(),
            std::back_inserter(tile_lookup));
  BOOST_FOREACH(const Tile& t, plate->batch_read(tile_lookup)) {
    trans_view_t image_data(t.hdr.transaction_id(), cache.get());
    boost::scoped_ptr<SrcImageResource> r(SrcMemoryImageResource::open(t.hdr.filetype(),&t.data->operator[](0), t.data->size()));
    read_image(image_data.second, *r);
    cmap[rowcol_t(t.hdr.row(),t.hdr.col())].push_back( image_data );
  }
}

composite_map_t build_map(HeaderList const& headers) {
  composite_map_t cmap;
  BOOST_FOREACH( TileHeader const& hdr, headers )
    cmap[rowcol_t(hdr.row(),hdr.col())].push_back(hdr);
  return cmap;
}

void initialize_albedo( HeaderList drg_headers,
                        ChannelAccumulator<MinMaxAccumulator<float32 > >& minmax_acc,
                        Options& opt ) {
  photk::AlbedoInitNRAccumulator<PixelGrayA<float32> > accum(opt.tile_size,opt.tile_size);

  composite_map_t drg_map = build_map( drg_headers );
  composite_map_t::const_iterator i = drg_map.begin(), end = drg_map.end();
  do {
    size_t size = 0;
    HeaderList sub_drg;

    for (; i != end; ++i) {
      if ( size + i->second.size() <= opt.tile_cache.size() ) {
        sub_drg.insert(sub_drg.end(),i->second.begin(),i->second.end());
        size += i->second.size();
      } else {
        break;
      }
    }

    // Load up tiles
    opt.tile_cache.index = 0;
    cache_map_t drg_cache;
    cache_consume_tiles( opt.drg_plate, sub_drg, drg_cache, opt.tile_cache );

    // Do the albedo work
    for ( cache_map_t::const_iterator map_it = drg_cache.begin();
          map_it != drg_cache.end(); map_it++ ) {

      BOOST_FOREACH( trans_view_t const& data, map_it->second )
        accum( data.second, opt.exposure_vec[data.first-1] );

      ImageView<PixelGrayA<float32> > result = accum.result();
      for_each_pixel(alpha_to_mask(result), minmax_acc);
      opt.albedo_plate->write_update( result, map_it->first.get<1>(),
                                      map_it->first.get<0>(), opt.level );
    }
  } while ( i != end );
}

void update_albedo( HeaderList drg_headers,
                    HeaderList albedo_headers,
                    ChannelAccumulator<MinMaxAccumulator<float32 > >& minmax_acc,
                    Options& opt ) {
  photk::AlbedoDeltaNRAccumulator<PixelGrayA<float32> > accum(opt.tile_size,opt.tile_size);
  opt.tile_cache.index = 0;

  // Load up all albedo, max 32 MB
  cache_map_t albedo_cache;
  cache_consume_tiles( opt.albedo_plate, albedo_headers, albedo_cache, opt.tile_cache );
  const size_t albedo_index = opt.tile_cache.index;

  // Start loading up chunks
  composite_map_t drg_map = build_map( drg_headers );
  composite_map_t::const_iterator i = drg_map.begin(), end = drg_map.end();
  do {
    size_t size = 0;
    HeaderList sub_drg;

    for (; i != end; ++i) {
      if ( size + i->second.size() <= opt.tile_cache.size() - albedo_index ) {
        sub_drg.insert(sub_drg.end(),i->second.begin(),i->second.end());
        size += i->second.size();
      } else {
        break;
      }
    }

    // Load up DRG tiles
    opt.tile_cache.index = albedo_index;
    cache_map_t drg_cache;
    cache_consume_tiles( opt.drg_plate, sub_drg, drg_cache, opt.tile_cache );

    // Do the albedo work
    for ( cache_map_t::const_iterator map_it = drg_cache.begin();
          map_it != drg_cache.end(); map_it++ ) {

      ImageView<PixelGrayA<float32> >& albedo =
        albedo_cache[map_it->first].back().second;

      BOOST_FOREACH( trans_view_t const& data, map_it->second )
        accum( data.second, albedo, opt.exposure_vec[data.first-1] );

      select_channel(albedo,0) += select_channel(accum.result(),0);
      for_each_pixel(alpha_to_mask(albedo), minmax_acc);
      opt.albedo_plate->write_update( albedo, map_it->first.get<1>(),
                                      map_it->first.get<0>(), opt.level );
    }

  } while ( i != end );
}

void handle_arguments( int argc, char *argv[], Options& opt ) {
  namespace po = boost::program_options;
  po::options_description general_options("");
  general_options.add_options()
    ("level,l", po::value(&opt.level)->default_value(-1), "Default is to process at lowest level.")
    ("job_id,j", po::value(&opt.job_id)->default_value(0), "")
    ("num_jobs,n", po::value(&opt.num_jobs)->default_value(1), "");
  general_options.add( photk::BaseOptionsDescription(opt) );

  po::options_description positional("");
  positional.add_options()
    ("ptk_url", po::value(&opt.ptk_url), "Input PTK url");

  po::positional_options_description positional_desc;
  positional_desc.add("ptk_url", 1);

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " <ptk-url>\n";

  po::variables_map vm =
    photk::check_command_line( argc, argv, opt, general_options,
                               positional, positional_desc, usage.str() );

  if ( opt.ptk_url == Url() )
    vw_throw( ArgumentErr() << "Missing project file url!\n"
              << usage.str() << general_options );

  // Open remote project file
  opt.remote_ptk.reset( new photk::RemoteProjectFile(opt.ptk_url) );
  opt.remote_ptk->get_project( opt.project_info );
  opt.remote_ptk->get_platefiles(opt.drg_plate,opt.albedo_plate,
                                 opt.reflect_plate);
  opt.tile_size = opt.albedo_plate->default_tile_size();
  opt.tile_cache = TileCache( calc_cache_tile_count( opt.albedo_plate ),
                              opt.tile_size );

  if ( opt.level < 0 )
    opt.level = opt.drg_plate->num_levels() - 1;
  if ( opt.level >= opt.drg_plate->num_levels() )
    vw_throw( ArgumentErr() << "Can't request level higher than available in DRG plate" );
  if ( opt.job_id >= opt.num_jobs ||
       opt.job_id < 0 || opt.num_jobs < 0 )
    vw_throw( ArgumentErr() << "Invalid job_id or num_jobs!");
}

int main( int argc, char *argv[] ) {

  Options opt;
  try {
    handle_arguments( argc, argv, opt );

    // Define work area
    opt.work_area =
      photk::split_square_to_equal_area( 0x1 << opt.level, opt.num_jobs )[opt.job_id];

    // Request current exposure values
    opt.exposure_vec.reserve( opt.project_info.num_cameras() );
    for ( int32 i = 0; i < opt.project_info.num_cameras(); i++ ) {
      photk::CameraMeta current_cam;
      opt.remote_ptk->get_camera( i, current_cam );
      opt.exposure_vec.push_back( current_cam.exposure_t() );
    }

    // Initialize pixval accumulator
    ChannelAccumulator<MinMaxAccumulator<float32 > > minmaxacc;

    // At a high level, request top tiles from DRG to figure out where
    // populated areas are. We'll never actually load these tiles.
    HeaderList l3_test_records =
      opt.drg_plate->search_by_region( opt.level - 3, opt.work_area / 8,
                                       TransactionRange(-1,-1) );

    if ( l3_test_records.size() == 0 )
      return 0;

    // Request transaction and write permissions
    opt.albedo_plate->transaction_begin("phoitalbedo [id="+boost::lexical_cast<std::string>(opt.job_id)+"]",-1);
    opt.albedo_plate->write_request();

    std::ostringstream job_prefix;
    job_prefix << "[Job " << opt.job_id << "/" << opt.num_jobs << "] ";
    vw_out() << job_prefix.str() << "Work area: " << opt.work_area << "\n";

    TerminalProgressCallback tpc("photometrytk", job_prefix.str());
    float tpc_inc = 1.0/float(l3_test_records.size());
    tpc.report_progress(0);
    BOOST_FOREACH( TileHeader const& l3_tile, l3_test_records ) {
      BBox2i query( l3_tile.col() * 8, l3_tile.row() * 8, 8, 8 );
      query.crop( opt.work_area );

      HeaderList drg_headers =
        opt.drg_plate->search_by_region( opt.level, query,
                                         TransactionRange(0,opt.project_info.num_cameras()+1) );
      HeaderList albedo_headers =
        opt.albedo_plate->search_by_region( opt.level, query,
                                            TransactionRange(-1,-1) );

      // Actually process jobs!
      if ( albedo_headers.empty() )
        initialize_albedo( drg_headers, minmaxacc, opt );
      else
        update_albedo( drg_headers, albedo_headers, minmaxacc, opt );

      tpc.report_incremental_progress(tpc_inc);
    }
    tpc.report_finished();

    // Finish up transaction writing
    opt.albedo_plate->write_complete();
    opt.albedo_plate->transaction_end(true);

    // Report the min max values we found
    try {
      opt.remote_ptk->add_pixvals(minmaxacc.minimum(),
                                  minmaxacc.maximum());
    } catch (ArgumentErr const& e ) {}

  } PHOTK_STANDARD_CATCHES;

  return 0;
}

// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

// What's this file supposed to do ?
//
// (Pho)tometry (It)eration Exposure (Time) Update
//
// With Reflectance
// .... see docs
//
// With out Relectance
//      T = T + sum((Ik-Tk*A)*A*Sk)/sum((A*Sk)^2)

#include <photk/Common.h>
#include <photk/Macros.h>
#include <photk/RemoteProjectFile.h>
#include <photk/TimeAccumulators.h>
#include <photk/PlateCommon.h>
#include <boost/foreach.hpp>

using namespace photk;
using namespace vw;
using namespace vw::platefile;

struct Options : photk::BaseOptions  {
  // Input
  Url ptk_url;

  // Settings
  int32 job_id, num_jobs, level;
  bool verbose;

  // Reused variables
  boost::scoped_ptr<photk::RemoteProjectFile> remote_ptk;
  photk::ProjectMeta project_info;
  photk::PlatePtr drg_plate, albedo_plate, reflect_plate;
  size_t tile_size;
  photk::TileCache tile_cache;
};

void handle_arguments( int argc, char *argv[], Options& opt ) {
  namespace po = boost::program_options;
  po::options_description general_options("");
  general_options.add_options()
    ("level,l", po::value(&opt.level)->default_value(-1), "Default is to process lowest level.")
    ("job_id,j", po::value(&opt.job_id)->default_value(0), "")
    ("num_jobs,n", po::value(&opt.num_jobs)->default_value(1), "")
    ("verbose", po::bool_switch(&opt.verbose), "Print out the changes happening to exposure time.");
  general_options.add( photk::BaseOptionsDescription(opt) );

    po::options_description positional("");
  positional.add_options()
    ("ptk_url",  po::value(&opt.ptk_url),  "Input PTK Url");

  po::positional_options_description positional_desc;
  positional_desc.add("ptk_url", 1);

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " <ptk-url>\n";

  po::variables_map vm =
    photk::check_command_line( argc, argv, opt, general_options,
                             positional, positional_desc, usage.str() );

  if ( opt.ptk_url == Url() )
    vw_throw( ArgumentErr() << "Missing project file url!\n"
              << usage.str() <<  general_options );

  // Open remote project file
  opt.remote_ptk.reset( new photk::RemoteProjectFile(opt.ptk_url) );
  opt.remote_ptk->get_project( opt.project_info );
  opt.remote_ptk->get_platefiles( opt.drg_plate, opt.albedo_plate,
                                  opt.reflect_plate );
  opt.tile_size = opt.albedo_plate->default_tile_size();
  opt.tile_cache = photk::TileCache( calc_cache_tile_count( opt.albedo_plate ),
                                     opt.tile_size );

  if ( opt.level < 0 )
    opt.level = opt.drg_plate->num_levels() - 1;
  if ( opt.level >= opt.drg_plate->num_levels() )
    vw_throw( ArgumentErr() << "Can't request level higher than available in DRG plate" );
  if ( opt.job_id >= opt.num_jobs ||
       opt.job_id < 0 || opt.num_jobs < 0 )
    vw_throw( ArgumentErr() << "Invalid job_id or num_jobs!");
}

int main( int argc, char* argv[] ) {
  Options opt;
  try {
    handle_arguments( argc, argv, opt );

    // Deciding what cameras
    int min_cam_idx =
      boost::numeric_cast<int>(float(opt.project_info.num_cameras()*opt.job_id)/
                               float(opt.num_jobs));
    int max_cam_idx =
      boost::numeric_cast<int>(float(opt.project_info.num_cameras()*(opt.job_id+1))/
                               float(opt.num_jobs));
    std::vector<TimeDeltaNRAccumulator> accumulators;
    std::vector<CameraMeta> camera_information;
    accumulators.reserve(max_cam_idx-min_cam_idx);
    camera_information.reserve(max_cam_idx-min_cam_idx);
    for( size_t i = 0; i < max_cam_idx-min_cam_idx; i++ ) {
      CameraMeta cam_info;
      opt.remote_ptk->get_camera(i+min_cam_idx,cam_info);

      // Fix any errors that might have worked in
      if ( std::isnan(cam_info.exposure_t()) )
        cam_info.set_exposure_t(1.0);

      camera_information.push_back(cam_info);
      accumulators.push_back(TimeDeltaNRAccumulator(camera_information.back().exposure_t()));
    }

    std::ostringstream job_prefix ;
    job_prefix << "[Job " << opt.job_id << "/" << opt.num_jobs << "] ";
    vw_out() << job_prefix.str() << "Camera range [" << min_cam_idx << " " << max_cam_idx << "]\n";

    size_t level_size = 0x1 << (opt.level - 3);
    HeaderList l3_test_records =
      opt.drg_plate->search_by_region( opt.level - 3,
                                       BBox2i(0,0,level_size,level_size),
                                       TransactionRange(-1,-1) );
    if ( l3_test_records.size() == 0 )
      vw_throw( LogicErr() << "There are no tile for which to do time correction with!" );

    // Iterate through geographical locations to start calculating the time deltas
    TerminalProgressCallback tpc("photometrytk", job_prefix.str());
    float tpc_inc = 1.0/float(l3_test_records.size());
    tpc.report_progress(0);
    BOOST_FOREACH( TileHeader const& l3_tile, l3_test_records ) {
      BBox2i query( l3_tile.col() * 8, l3_tile.row() * 8, 8, 8 );

      HeaderList albedo_headers =
        opt.albedo_plate->search_by_region( opt.level, query,
                                            TransactionRange(-1,-1) );
      cache_map_t albedo_cache;
      opt.tile_cache.index = 0;
      cache_consume_tiles( opt.albedo_plate, albedo_headers, albedo_cache, opt.tile_cache );
      const size_t albedo_index = opt.tile_cache.index;

      HeaderList drg_headers =
        opt.drg_plate->search_by_region( opt.level, query,
                                         TransactionRange(min_cam_idx+1,max_cam_idx) );
      const size_t drg_size = drg_headers.size();
      const size_t free_space = opt.tile_cache.size() - albedo_index;

      HeaderList::const_iterator drg_section_start = drg_headers.begin();
      HeaderList::const_iterator drg_section_end   = drg_headers.begin();
      for ( size_t load_index = 0; load_index < drg_size; load_index += free_space ) {
        // Reset the cache write location
        opt.tile_cache.index = albedo_index;

        size_t load_amount = free_space;
        if ( load_amount + load_index > drg_size )
          load_amount = drg_size - load_index;
        drg_section_start = drg_section_end;
        std::advance( drg_section_end, load_amount );
        HeaderList sub_drg_headers(drg_section_start, drg_section_end);
        std::vector<hdr_view_t> sub_drg_tiles;
        cache_consume_tiles( opt.drg_plate, sub_drg_headers,
                             sub_drg_tiles, opt.tile_cache );

        // Actually run the time accumulators on the retrieved data
        BOOST_FOREACH( hdr_view_t const& data, sub_drg_tiles ) {
          for_each_pixel( data.second,
                          albedo_cache[rowcol_t(data.first.row(),
                                                data.first.col())].front().second,
                          accumulators[data.first.transaction_id()-1-min_cam_idx] );
        }
      }

      tpc.report_incremental_progress(tpc_inc);
    }
    tpc.report_finished();

    // Modify and write the new camera informations
    for( size_t i = 0; i < accumulators.size(); i++ ) {
      if ( opt.verbose )
        vw_out() << "  Camera [" << i+min_cam_idx << "] : " << camera_information[i].exposure_t() << " + " << accumulators[i].value() << "\n";
      if ( !std::isnan(accumulators[i].value()) )
        camera_information[i].set_exposure_t( camera_information[i].exposure_t() +
                                              accumulators[i].value() );
      opt.remote_ptk->set_camera(i+min_cam_idx,
                                 camera_information[i]);
    }

    if ( opt.job_id == 0 ) {
      if ( opt.verbose )
        vw_out() << "  Iteration incremented to: " << opt.project_info.current_iteration() + 1 << "\n";
      opt.remote_ptk->set_iteration(opt.project_info.current_iteration()+1);
    }

  } PHOTK_STANDARD_CATCHES;

  return 0;
}

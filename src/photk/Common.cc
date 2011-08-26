// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <photk/Common.h>
#include <photk/config.h>
#include <vw/Math/BBox.h>
#include <vw/config.h>
#include <boost/filesystem/path.hpp>
using namespace vw;

namespace photk {

  BaseOptions::BaseOptions() {
#if defined(VW_HAS_BIGTIFF) && VW_HAS_BIGTIFF == 1
    gdal_options["COMPRESS"] = "LZW";
#else
    gdal_options["COMPRESS"] = "NONE";
    gdal_options["BIGTIFF"] = "NO";
#endif
    raster_tile_size =
      vw::Vector2i(vw::vw_settings().default_tile_size(),
                   vw::vw_settings().default_tile_size());
  }

  BaseOptionsDescription::BaseOptionsDescription( BaseOptions& opt ) {
    namespace po = boost::program_options;
    (*this).add_options()
      ("threads", po::value(&opt.num_threads)->default_value(0),
       "Select the number of processors (threads) to use.")
      ("no-bigtiff", "Tell GDAL to not create bigtiffs.")
      ("version,v", "Display the version of software.")
      ("help,h", "Display this help message");
  }

  boost::program_options::variables_map
  check_command_line( int argc, char *argv[], BaseOptions& opt,
                      boost::program_options::options_description const& public_options,
                      boost::program_options::options_description const& hidden_options,
                      boost::program_options::positional_options_description const& positional,
                      std::string const& help ) {
    namespace po = boost::program_options;
    po::variables_map vm;
    try {
      po::options_description all_options;
      all_options.add(public_options).add(hidden_options);
      po::store( po::command_line_parser( argc, argv ).options(all_options).positional(positional).run(), vm );
      po::notify( vm );
    } catch (po::error const& e) {
      vw::vw_throw( vw::ArgumentErr() << "Error parsing input:\n"
                    << e.what() << "\n" << help << "\n" << public_options );
    }
    // We really don't want to use BIGTIFF unless we have to. It's
    // hard to find viewers for bigtiff.
    if ( vm.count("no-bigtiff") ) {
      opt.gdal_options["BIGTIFF"] = "NO";
    } else {
      opt.gdal_options["BIGTIFF"] = "IF_SAFER";
    }
    if ( vm.count("help") )
      vw::vw_throw( vw::ArgumentErr() << help << "\n" << public_options );
    if ( vm.count("version") )
      vw::vw_throw( vw::ArgumentErr() << PHOTK_PACKAGE_STRING << "\n\n"
                    << "Built against:\n  " << VW_PACKAGE_STRING << "\n  BOOST "
                    << PHOTK_BOOST_VERSION << "\n");
    if ( opt.num_threads != 0 ) {
      vw::vw_out() << "\t--> Setting number of processing threads to: "
                   << opt.num_threads << std::endl;
      vw::vw_settings().set_default_num_threads(opt.num_threads);
    }

    return vm;
  }

  bool has_cam_extension( std::string input ) {
    boost::filesystem::path ipath( input );
    std::string ext = ipath.extension();
    if ( ext == ".cahvor" || ext == ".cahv" ||
         ext == ".pin" || ext == ".pinhole" ||
         ext == ".tsai" || ext == ".cmod" ||
         ext == ".cahvore" )
      return true;
    return false;
  }

  // This algorithm is lifted from, "Cutting circles and squares into
  // equal pieces" by Bose et. al.
  std::vector<vw::BBox2i>
  split_square_to_equal_area( int32 const& side_length,
                              size_t const& k ) {
    using namespace vw;
    std::vector<BBox2i> result;
    result.reserve( k );

    VW_ASSERT( uint64(side_length)*uint64(side_length) >= uint64(k),
               ArgumentErr() << "Requesting more regions then are pixels in square. K = " << k << ", Side Length = " << side_length );
    VW_ASSERT( side_length >= k,
               ArgumentErr() << "K must be smaller or equal to side_length." );

    int32 a = int32(floor(sqrt(float(k))));
    int32 r = k - a*a;
    int32 r_p_h_num = 0, r_pp_h_num = 0;
    if ( 0 <= r && r <= a ) {
      r_p_h_num = r;
      r_pp_h_num = a - r;
    } else if ( a + 1 <= r && r <= 2*a ) {
      r_p_h_num = r - a;
      r_pp_h_num = 2*a - r + 1;
    } else {
      vw_throw( LogicErr() << "Split square has hit undiscovered error." );
    }
    int32 region_split =
      int32(float(side_length * (a+1) * r_p_h_num)/float(k));
    int32 r_p_h = r_p_h_num == 0 ? 0 : ceil(float(region_split) / float(r_p_h_num));
    int32 r_pp_h = r_pp_h_num == 0 ? 0 : int32(ceil(float(side_length - region_split) / float(r_pp_h_num)));
    int32 r_p_v = int32(ceil(float(side_length)/float(a + 1)));
    int32 r_pp_v = int32(ceil(float(side_length)/float(a)));

    int32 x,y;
    if ( r_p_h != 0 ) {
      // Adding boxes in the R' region
      for (y = 0; y < side_length - r_p_v; y += r_p_v) {
        for (x = 0; x < region_split - r_p_h; x += r_p_h)
          result.push_back(BBox2i(x,y,r_p_h,r_p_v));
        result.push_back(BBox2i(x,y,region_split-x,r_p_v));
      }
      for (x = 0; x < region_split - r_p_h; x += r_p_h )
        result.push_back(BBox2i(x,y,r_p_h,side_length-y));
      result.push_back(BBox2i(x,y,region_split-x,side_length-y));
    }

    // Adding boxes in the R'' region
    if ( r_pp_h != 0 ) {
      for (y = 0; y < side_length - r_pp_v; y += r_pp_v) {
        for (x = region_split; x < side_length-r_pp_h; x += r_pp_h)
          result.push_back(BBox2i(x,y,r_pp_h,r_pp_v));
        result.push_back(BBox2i(x,y,side_length-x,r_pp_v));
      }
      for (x = region_split; x < side_length-r_pp_h; x += r_pp_h)
        result.push_back(BBox2i(x,y,r_pp_h,side_length-y));
      result.push_back(BBox2i(x,y,side_length-x,side_length-y));
    }

    return result;
  }

}

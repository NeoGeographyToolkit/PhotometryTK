// __BEGIN_LICENSE__
//  Copyright (c) 2009-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NGT platform is licensed under the Apache License, Version 2.0 (the
//  "License"); you may not use this file except in compliance with the
//  License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__

/// \file dem_mosaic.cc
///

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <limits>
using namespace std;

#include <vw/FileIO.h>
#include <vw/Image.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace vw;

#include <photk/Misc.h>
#include <photk/Weights.h>
using namespace photometry;

#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

// A tool to mosaic and blend DEMs.
// To do: Make it work for other projections except longlat.

class DemMosaicView: public ImageViewBase<DemMosaicView>{
  int m_cols, m_rows;
  vector<ModelParams> & m_modelParamsArray;
  cartography::GeoReference m_out_georef;
  vector<double> m_nodata_values;
  double m_out_nodata_value;

public:
  DemMosaicView(int cols, int rows, vector<ModelParams> & modelParamsArray,
                cartography::GeoReference const& out_georef,
                vector<double> const& nodata_values, double out_nodata_value):
    m_cols(cols), m_rows(rows), m_modelParamsArray(modelParamsArray),
    m_out_georef(out_georef), m_nodata_values(nodata_values),
    m_out_nodata_value(out_nodata_value){}
  
  typedef float pixel_type;
  typedef pixel_type result_type;
  typedef PixelMask<float> masked_pixel_type;
  typedef ProceduralPixelAccessor<DemMosaicView> pixel_accessor;
  
  inline int32 cols() const { return m_cols; }
  inline int32 rows() const { return m_rows; }
  inline int32 planes() const { return 1; }
  
  inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

  inline pixel_type operator()( double/*i*/, double/*j*/, int32/*p*/ = 0 ) const {
    vw_throw(NoImplErr() << "DemMosaicView::operator()(...) is not implemented");
    return pixel_type();
  }

  typedef CropView<ImageView<pixel_type> > prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {

    ImageView<pixel_type> tile(bbox.width(), bbox.height());
    ImageView<float>      weights(bbox.width(), bbox.height());
    fill( tile, m_out_nodata_value );
    fill( weights, 0.0 );

    for (int dem_iter = 0; dem_iter < (int)m_modelParamsArray.size(); dem_iter++){

      std::string curr_file = m_modelParamsArray[dem_iter].inputFilename;
      cartography::GeoReference georef;
      bool is_good = cartography::read_georeference( georef, curr_file );
      if (!is_good)
        vw_throw(ArgumentErr() << "DEM: " << curr_file
                 << " does not have a georeference.\n");

      double nodata_value = m_nodata_values[dem_iter];
      
      DiskImageView<float> curr_disk_dem(curr_file);

      // The tile corners as pixels in curr_dem
      Vector2 b = georef.lonlat_to_pixel
        (m_out_georef.pixel_to_lonlat(bbox.min()));
      Vector2 e = georef.lonlat_to_pixel
        (m_out_georef.pixel_to_lonlat(bbox.max()));
      if (b[0] > e[0]) std::swap(b[0], e[0]);
      if (b[1] > e[1]) std::swap(b[1], e[1]);
      b = floor(b); e = ceil(e);
      BBox2i curr_box(b[0], b[1], e[0] - b[0], e[1] - b[1]);
      curr_box.expand(BilinearInterpolation::pixel_buffer + 1);
      curr_box.crop(bounding_box(curr_disk_dem));
      if (curr_box.empty()) continue;
            
      // Read the whole chunk into memory, to speed things up
      ImageView<float> curr_dem = crop(curr_disk_dem, curr_box);
      ImageViewRef< PixelMask<float> > interp_dem
        = interpolate(create_mask(curr_dem, nodata_value),
                      BilinearInterpolation(), ConstantEdgeExtension());

      for (int c = 0; c < bbox.width(); c++){
        for (int r = 0; r < bbox.height(); r++){
          Vector2 out_pix(c +  bbox.min().x(), r +  bbox.min().y());
          Vector2 in_pix = georef.lonlat_to_pixel
            (m_out_georef.pixel_to_lonlat(out_pix));
          double x = in_pix[0] - curr_box.min().x();
          double y = in_pix[1] - curr_box.min().y();
          // below must use x <= cols()-1 as x is double
          if (x >= 0 && x <= interp_dem.cols()-1 &&
              y >= 0 && y <= interp_dem.rows()-1 ){
            PixelMask<float> pval = interp_dem(x, y);
            if (! is_valid(pval)) continue;
            double val = pval.child();
            double wt = ComputeLineWeightsHV(in_pix,
                                             m_modelParamsArray[dem_iter]);
            if (wt <= 0) continue;
            if ( tile(c, r) == m_out_nodata_value || isnan(tile(c, r)) )
              tile(c, r) = 0;
            tile(c, r) += wt*val;
            weights(c, r) += wt;
            
          }
        }
      }
      
    } // end iterating over DEMs

    // Divide by the weights
    for (int c = 0; c < bbox.width(); c++){
      for (int r = 0; r < bbox.height(); r++){
        if ( weights(c, r) > 0 )
          tile(c, r) /= weights(c, r);        
      }
    }
    
    return prerasterize_type(tile, -bbox.min().x(), -bbox.min().y(),
                             cols(), rows() );
  }

  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};

int main( int argc, char *argv[] ) {

  for (int s = 0; s < argc; s++) vw_out() << argv[s] << " ";
  vw_out() << endl;
  
  try{
    
    string dem_list_file, out_dem_file;
    double mpp = 0.0;
    int num_threads = 0;
    po::options_description general_options("Options");
    general_options.add_options()
      ("dem-list-file,l", po::value<string>(&dem_list_file),
       "List of DEM files to mosaic, one per line.")
      ("mpp", po::value<double>(&mpp),
       "Output DEM resolution in meters per pixel.")
      ("output-dem,o", po::value<string>(&out_dem_file),
       "Output DEM file name.")
      ("threads", po::value<int>(&num_threads),
       "Number of threads to use.")
      ("help,h", "Display this help message.");
    ;

    // Parse options
    po::options_description hidden_options("");
    po::options_description options("Allowed Options");
    options.add(general_options).add(hidden_options);
    po::positional_options_description p;
    ostringstream usage;
    usage << "A tool for mosaicking DEMs without seams.\n\n";
    usage << general_options << endl;
    po::variables_map vm;
    po::store( po::command_line_parser( argc, argv ).options(options)
               .positional(p).run(), vm );
    po::notify( vm );
    if  (vm.count("help") || argc <= 1){
      vw_out() << usage.str() << "\n";
      return 1;
    }

    // Error checking
    if (dem_list_file == "")
      vw_throw(ArgumentErr() << "No list of DEMs was specified.\n");
    if (mpp <= 0.0)
      vw_throw(ArgumentErr() << "The output resolution in meters per "
               << "pixel must be set and positive.\n");
    if (out_dem_file == "")
      vw_throw(ArgumentErr() << "No output DEM file was specified.\n");
    if (num_threads == 0)
      vw_throw(ArgumentErr() << "The number of threads must be set and "
               << "positive.\n");

    // Read the DEMs to mosaic
    vector<string> dem_files;
    ifstream is(dem_list_file.c_str());
    string file;
    while (is >> file) dem_files.push_back(file);
    if (dem_files.empty())
      vw_throw(ArgumentErr() << "No DEM files to mosaic.\n");
    is.close();
    
    // Read nodata from first DEM
    double out_nodata_value = numeric_limits<double>::quiet_NaN();
    {
      DiskImageResourceGDAL in_rsrc(dem_files[0]);
      if ( in_rsrc.has_nodata_read() ) out_nodata_value = in_rsrc.nodata_read();
    }
    
    // Find the lon-lat bounding box of all DEMs
    double big = numeric_limits<double>::max();
    Vector4 ll_bbox(big, -big, big, -big);
    for (int dem_iter = 0; dem_iter < (int)dem_files.size(); dem_iter++){
      std::string curr_file = dem_files[dem_iter];
      Vector4 corners = getImageCorners(curr_file);
      ll_bbox[0] = std::min(ll_bbox[0], corners[0]);
      ll_bbox[1] = std::max(ll_bbox[1], corners[1]);
      ll_bbox[2] = std::min(ll_bbox[2], corners[2]);
      ll_bbox[3] = std::max(ll_bbox[3], corners[3]);
    }

    // Convert from meters/pixel to pixels/degree, and form the
    // georef. We will use the geographic projection only, so this may
    // not work well around poles.
    cartography::GeoReference georef;
    bool is_good = read_georeference(georef, dem_files[0]);
    if (!is_good){
      std::cerr << "No georeference found in " << dem_files[0] << std::endl;
      exit(1);
    }
    double spacing = 360.0*mpp/( 2*M_PI*georef.datum().semi_major_axis() );
    georef.set_geographic(); // wipe original georef
    Matrix<double,3,3> transform;
    transform.set_identity();
    transform(0, 0) = spacing;
    transform(1, 1) = -spacing;
    georef.set_transform(transform);
    Vector2 beg_pix = georef.lonlat_to_pixel(Vector2(ll_bbox[0], ll_bbox[3]));
    georef = crop(georef, beg_pix[0], beg_pix[1]);

    // Image size
    Vector2 end_pix = georef.lonlat_to_pixel(Vector2(ll_bbox[1], ll_bbox[2]));
    int cols = (int)ceil(end_pix[0]);
    int rows = (int)ceil(end_pix[1]);

    // Compute the weights, and store the no-data values
    vector<ModelParams> modelParamsArray;
    vector<double> nodata_values;
    for (int dem_iter = 0; dem_iter < (int)dem_files.size(); dem_iter++){
      modelParamsArray.push_back(ModelParams());
      modelParamsArray[dem_iter].inputFilename = dem_files[dem_iter];
      ComputeDEMCenterLines(modelParamsArray[dem_iter]);

      double curr_nodata_value = out_nodata_value;
      DiskImageResourceGDAL in_rsrc(dem_files[dem_iter]);
      if ( in_rsrc.has_nodata_read() ) curr_nodata_value = in_rsrc.nodata_read();
      nodata_values.push_back(curr_nodata_value);
    }

    // Form the mosaic and write it to disk
    vw_out() << "Writing: " << out_dem_file << std::endl;
    vw::vw_settings().set_default_num_threads(num_threads);
    DiskImageResourceGDAL::Options gdal_options;
    gdal_options["COMPRESS"] = "LZW";
    ImageViewRef<float> out_dem = DemMosaicView(cols, rows, modelParamsArray,
                                                georef, nodata_values, out_nodata_value);
    DiskImageResourceGDAL rsrc(out_dem_file, out_dem.format(), Vector2i(256, 256),
                               gdal_options);
    if (!isnan(out_nodata_value))
      rsrc.set_nodata_write(out_nodata_value);
    write_georeference(rsrc, georef);
    block_write_image(rsrc, out_dem, TerminalProgressCallback("{Core}",
                                                              "Processing:"));
    
  } catch ( const exception& e ) {
    cerr << "\n\nError: " << e.what() <<  endl;
    return 1;                                     
  }

  return 0;
}

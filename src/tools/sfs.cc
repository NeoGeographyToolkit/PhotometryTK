//__BEGIN_LICENSE__
//  Copyright (c) 2009-2012, United States Government as represented by the
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

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
using namespace std;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
#include <photk/Photometry.h>

// ISIS
#include "Camera.h"
#include "CameraFactory.h"
#include "TProjection.h"
#include "ProjectionFactory.h"
#include "ProcessRubberSheet.h"
#include "IException.h"
#include "Pvl.h"
#include "IString.h"
#include "PushFrameCameraDetectorMap.h"
#include "Target.h"
#include "CameraDistortionMap.h"
#include "CameraFocalPlaneMap.h"
#include "SurfacePoint.h"
#include "Latitude.h"
#include "Longitude.h"

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace photometry;

// Compute mean and standard deviation of an image
template <class ImageT>
void compute_image_stats(ImageT const& I, double & mean, double & stdev){
  mean  = 0;
  stdev = 0;
  double sum = 0.0, sum2 = 0.0, count = 0.0;
  for (int col = 0; col < I.cols(); col++){
    for (int row = 0; row < I.rows(); row++){
      if (!is_valid(I(col, row))) continue;
      count++;
      double val = I(col, row).child();
      sum += val;
      sum2 += val*val;
    }
  }

  if (count > 0){
    mean = sum/count;
    stdev = sqrt( sum2/count - mean*mean );
  }
  
}

// Write an image together with georeference and nodata value.
template<class ImgT>
void write_image_nodata(std::string imgfile,
                        ImgT const& img, GeoReference geo, float nodata_val){
  
  std::cout << "Writing: " << imgfile << std::endl;
  DiskImageResourceGDAL::Options gdal_options;
  gdal_options["COMPRESS"] = "LZW";
  DiskImageResourceGDAL rsrc(imgfile, img.format(), Vector2i(256, 256),
                             gdal_options);
  rsrc.set_nodata_write(nodata_val);
  write_georeference(rsrc, geo);
  write_image(rsrc, img,
              TerminalProgressCallback("{Core}", "Processing:"));
}

double calc_max_abs(ImageView<PixelGray<float> >& DEM,
                    float DEMnodata){
  double mx = 0.0;
  for (int col = 0; col < DEM.cols(); col++){
    for (int row = 0; row < DEM.rows(); row++){
      if (DEM(col, row) == DEMnodata) continue;
      mx = std::max(mx, (double)fabs(DEM(col, row)));
    }
  }

  return mx;
}

void calc_min_max(ImageView<PixelGray<float> >& DEM,
                  float DEMnodata,
                  double & mn, double & mx){
  mn = DBL_MAX;
  mx = -DBL_MAX;
  for (int col = 0; col < DEM.cols(); col++){
    for (int row = 0; row < DEM.rows(); row++){
      if (DEM(col, row) == DEMnodata) continue;
      mx = std::max(mx, (double)DEM(col, row));
      mn = std::min(mn, (double)DEM(col, row));
    }
  }

}

// Utilities
string print_vec(vector<double> const& vec){
  ostringstream oss;
  oss.precision(18);
  oss << "(";
  for (int s = 0; s < (int)vec.size()-1; s++){
    oss << vec[s] << ", ";
  }
  if (!vec.empty()) oss << vec[vec.size()-1];
  oss << ")";
  return oss.str();
}

vector<double> diff(vector<double> const& vec1, vector<double> const& vec2){
  vector<double> ans = vec1;
  for (int i = 0; i < (int)vec1.size(); i++){
    ans[i] -= vec2[i];
  }
  return ans;
}

vector<double> prod(double val, vector<double> const& vec){
  vector<double> ans = vec;
  for (int i = 0; i < (int)vec.size(); i++){
    ans[i] *= val;
  }
  return ans;
}

vector<double> div(vector<double> const& vec, double val){
  vector<double> ans = vec;
  for (int i = 0; i < (int)vec.size(); i++){
    ans[i] /= val;
  }
  return ans;
}

vector<double> normalize(vector<double> const& vec){
  double len = 0;
  for (int i = 0; i < (int)vec.size(); i++){
    len += vec[i]*vec[i];
  }
  len = sqrt(len);
  vector<double> ans = vec;
  for (int i = 0; i < (int)vec.size(); i++){
    ans[i] /= len;
  }
  return ans;
}

// Abstract base class
class IsisInterface {
public:
  IsisInterface( string const& file );
  virtual ~IsisInterface(); // Can't be declared here since we have
  // incomplete types from Isis.

  virtual string type() = 0;
  static IsisInterface* open( string const& filename );

  // Standard Methods
  //------------------------------------------------------

  // These are the standard camera request; IsisInterface allows for
  // them to be customized for the type of camera so that they are
  // fast and not too full of conditionals.

  // Note: The smallest pixel is (0, 0), not (1, 1)!
  virtual vector<double>
  point_to_pixel( vector<double> const& point ) const = 0;
  virtual vector<double>
  pixel_to_vector( vector<double> const& pix ) const = 0;
  virtual vector<double>
  camera_center() const = 0;

protected:
  // Standard Variables
  //------------------------------------------------------
  Isis::Pvl * m_label;
  Isis::Camera * m_camera;
  Isis::Cube * m_cube;

};

// ISIS frame camera class
class IsisInterfaceFrame : public IsisInterface {
  
public:
  IsisInterfaceFrame( string const& filename );
  virtual ~IsisInterfaceFrame(){}
  
  virtual string type()  { return "Frame"; }
  
  // Standard Methods
  //-------------------------------------------------
  
  // Note: The smallest pixel is (0, 0), not (1, 1)!
  virtual vector<double>
  point_to_pixel( vector<double> const& point ) const;
  virtual vector<double>
  pixel_to_vector( vector<double> const& pix ) const;
  virtual vector<double>
  camera_center() const;
  
protected:
  
  // Custom Variables
  Isis::CameraDistortionMap *m_distortmap;
  Isis::CameraFocalPlaneMap *m_focalmap;
  Isis::CameraDetectorMap   *m_detectmap;
  mutable Isis::AlphaCube   m_alphacube;
  
  vector<double> m_center;
};

// ISIS linescan camera class
class IsisInterfaceLineScan : public IsisInterface {

public:
  IsisInterfaceLineScan( std::string const& filename );

  virtual ~IsisInterfaceLineScan() {}

  virtual std::string type()  { return "LineScan"; }

  // Standard Methods
  //-------------------------------------------------

  // Note: The smallest pixel is (0, 0), not (1, 1)!
  virtual vector<double>
  point_to_pixel( vector<double> const& point ) const;
  virtual vector<double>
  pixel_to_vector( vector<double> const& pix ) const;
  virtual vector<double>
  camera_center() const;

protected:

  // Custom Variables
  Isis::CameraDistortionMap *m_distortmap;
  Isis::CameraFocalPlaneMap *m_focalmap;
  Isis::CameraDetectorMap   *m_detectmap;
  mutable Isis::AlphaCube   m_alphacube; // Doesn't use const
};

// IsisInterface base class implementation

IsisInterface::IsisInterface( string const& file ) {
  // Opening labels and camera
  Isis::FileName ifilename( QString::fromStdString(file) );
  m_label = new Isis::Pvl();
  m_label->read( ifilename.expanded() );

  // Opening Isis::Camera
  m_cube = new Isis::Cube(QString::fromStdString(file));
  m_camera = Isis::CameraFactory::Create( *m_cube );
}

IsisInterface::~IsisInterface(){
  delete m_label;
  delete m_camera;
  delete m_cube;
}

IsisInterface* IsisInterface::open( string const& filename ) {
  IsisInterface* result;

  // Opening Labels (This should be done somehow though labels)
  Isis::FileName ifilename( QString::fromStdString(filename) );
  Isis::Pvl label;
  label.read( ifilename.expanded() );

  Isis::Cube tempCube(QString::fromStdString(filename));
  Isis::Camera* camera = Isis::CameraFactory::Create( tempCube );

  switch ( camera->GetCameraType() ) {
  case 0:
    // Framing Camera
    if ( camera->HasProjection() ){
      throw string("Don't support Isis Camera Type with projection at this moment");
      //result = new IsisInterfaceMapFrame( filename );
    }
    else{
      result = new IsisInterfaceFrame( filename );
    }
    break;
  case 2:
    // Linescan Camera
    if ( camera->HasProjection() ){
      throw string("Don't support Projected Linescan Isis Camera Type at this moment");
    }
    else{
      result = new IsisInterfaceLineScan( filename );
    }
    break;
  default:
    throw string("Don't support unknown Isis Camera Type at this moment");
  }
  return result;
}

// IsisInterfaceFrame class implementation

IsisInterfaceFrame::IsisInterfaceFrame( string const& filename ) :
  IsisInterface(filename), m_alphacube( *m_cube ) {

  // Gutting Isis::Camera
  m_distortmap = m_camera->DistortionMap();
  m_focalmap   = m_camera->FocalPlaneMap();
  m_detectmap  = m_camera->DetectorMap();

  // Calculating Center (just once)
  m_center.resize(3);
  m_camera->instrumentPosition(&m_center[0]);
  for (int i = 0; i < (int)m_center.size(); i++){
    m_center[i] *= 1000; // convert to meters
  }
}

vector<double>
IsisInterfaceFrame::point_to_pixel( vector<double> const& point ) const{

  // Note: The smallest pixel is (0, 0), not (1, 1)!

  vector<double> look = normalize( diff(point, m_center) );

  vector<double> lookB_copy(3);
  std::copy( look.begin(), look.end(), lookB_copy.begin() );
  lookB_copy = m_camera->bodyRotation()->J2000Vector(lookB_copy);
  lookB_copy = m_camera->instrumentRotation()->ReferenceVector(lookB_copy);
  std::copy( lookB_copy.begin(), lookB_copy.end(), look.begin() );
  look = div(prod(m_camera->FocalLength(), look), fabs(look[2]) );

  // Back Projecting
  m_distortmap->SetUndistortedFocalPlane( look[0], look[1] );
  m_focalmap->SetFocalPlane( m_distortmap->FocalPlaneX(),
                             m_distortmap->FocalPlaneY() );
  m_detectmap->SetDetector( m_focalmap->DetectorSample(),
                            m_focalmap->DetectorLine() );

  vector<double> ans(2);
  ans[0] = m_alphacube.BetaSample( m_detectmap->ParentSample() ) - 1;
  ans[1] = m_alphacube.BetaLine( m_detectmap->ParentLine() ) - 1;
  return ans;
}

vector<double> 
IsisInterfaceFrame::pixel_to_vector( vector<double> const& px ) const {

  // Note: The smallest pixel is (0, 0), not (1, 1)!

  m_detectmap->SetParent( m_alphacube.AlphaSample(px[0]+1),
                          m_alphacube.AlphaLine(px[1]+1));
  m_focalmap->SetDetector( m_detectmap->DetectorSample(),
                           m_detectmap->DetectorLine() );
  m_distortmap->SetFocalPlane( m_focalmap->FocalPlaneX(),
                               m_focalmap->FocalPlaneY() );

  vector<double> look(3);
  look[0] = m_distortmap->UndistortedFocalPlaneX();
  look[1] = m_distortmap->UndistortedFocalPlaneY();
  look[2] = m_distortmap->UndistortedFocalPlaneZ();

  look = normalize( look );
  look = m_camera->instrumentRotation()->J2000Vector(look);
  look = m_camera->bodyRotation()->ReferenceVector(look);
  
  return look;
}

vector<double> IsisInterfaceFrame::camera_center() const{
  return m_center;
}

// IsisInterfaceLineScan class implementation
IsisInterfaceLineScan::IsisInterfaceLineScan( std::string const& filename ):
  IsisInterface(filename), m_alphacube( *m_cube ) {

  // Gutting Isis::Camera
  m_distortmap = m_camera->DistortionMap();
  m_focalmap   = m_camera->FocalPlaneMap();
  m_detectmap  = m_camera->DetectorMap();
}

vector<double>
IsisInterfaceLineScan::point_to_pixel( vector<double> const& point ) const{
  
  // Note: The smallest pixel is (0, 0), not (1, 1)!
  using namespace Isis;
  SurfacePoint surfP(Displacement(point[0], Displacement::Meters),
                     Displacement(point[1], Displacement::Meters),
                     Displacement(point[2], Displacement::Meters));
  
  vector<double> pix(2);
  if(m_camera->SetUniversalGround(surfP.GetLatitude().degrees(),
                                  surfP.GetLongitude().degrees(),
                                  surfP.GetLocalRadius().meters()
                                  )) {
    pix[0] = m_camera->Sample()-1;
    pix[1] = m_camera->Line()-1;

    return pix;
  }
  
  throw string("Could not project point into the camera.");
  return pix;
}

vector<double> 
IsisInterfaceLineScan::pixel_to_vector( vector<double> const& px ) const {

  // Note: The smallest pixel is (0, 0), not (1, 1)!

  m_detectmap->SetParent( m_alphacube.AlphaSample(px[0]+1),
                          m_alphacube.AlphaLine(px[1]+1));
  m_focalmap->SetDetector( m_detectmap->DetectorSample(),
                           m_detectmap->DetectorLine() );
  m_distortmap->SetFocalPlane( m_focalmap->FocalPlaneX(),
                               m_focalmap->FocalPlaneY() );

  vector<double> look(3);
  look[0] = m_distortmap->UndistortedFocalPlaneX();
  look[1] = m_distortmap->UndistortedFocalPlaneY();
  look[2] = m_distortmap->UndistortedFocalPlaneZ();

  look = normalize( look );
  look = m_camera->instrumentRotation()->J2000Vector(look);
  look = m_camera->bodyRotation()->ReferenceVector(look);

  return look;
}

vector<double> IsisInterfaceLineScan::camera_center() const{
  throw string("Computation of camera center is not implemented")
    + " for linescan cameras.";
}

// Convolve a DEM with exp(-sigma*x^2). The input DEM must
// have its invalid pixels masked.

template <class ImageT>
class BlurDEM:
  public ImageViewBase< BlurDEM<ImageT> > {
  ImageT m_img;
  int m_search_dist; // half of size of window to convolve with
  ImageView<double> m_gauss_kernel;
  typedef typename ImageT::pixel_type PixelT;

public:

  typedef PixelT pixel_type;
  typedef PixelT result_type;
  typedef ProceduralPixelAccessor<BlurDEM> pixel_accessor;

  BlurDEM( ImageViewBase<ImageT> const& img,
                     double blur_sigma) :
    m_img(img.impl()) {
    VW_ASSERT(blur_sigma > 0,
              ArgumentErr() << "Expecting positive sigma.");

    // Cut the gaussian exp(-sigma*x^2) where its value is 'scale'.
    double scale = 0.001;
    m_search_dist = (int)ceil(sqrt(-log(scale)/blur_sigma));
    std::cout << "Blur search distance is " << m_search_dist << std::endl;

    // The gaussian kernel
    int h = m_search_dist;
    m_gauss_kernel.set_size(2*h+1, 2*h+1);
    for (int c = 0; c < m_gauss_kernel.cols(); c++){
      for (int r = 0; r < m_gauss_kernel.rows(); r++){
        double r2 = double(c-h)*(c-h) + double(r-h)*(r-h);
        m_gauss_kernel(c, r) = exp(-blur_sigma*r2);
      }
    }
  
  }

  inline int32 cols() const { return m_img.cols(); }
  inline int32 rows() const { return m_img.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return pixel_accessor(*this); }

  inline result_type operator()( size_t i, size_t j, size_t p=0 ) const {
    vw_throw( NoImplErr() << "BlurDEM: operator() not implemented.\n" );
  }
  
  
  typedef CropView< ImageView<PixelT> > prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {

    // Crop into an expanded box as to have enough pixels to do
    // the blurring at every pixel in the current box.
    int h = m_search_dist; // shorten
    BBox2i biased_box = bbox;
    biased_box.expand(h+1);
    biased_box.crop(bounding_box(m_img));
    ImageView<PixelT> img( crop( m_img, biased_box ) );
    ImageView<PixelT> filled_img = copy(img);

    int nc = img.cols(), nr = img.rows(); // shorten
    for (int row = 0; row < nr; row++){
      for (int col = 0; col < nc; col++){
        PixelT V; V.validate();
        double sum = 0.0;
        for (int c = std::max(col-h, 0); c <= std::min(col+h, nc-1); c++){
          for (int r = std::max(row-h, 0); r <= std::min(row+h, nr-1); r++){
            if (!is_valid(img(c, r))) continue;
            double wt = m_gauss_kernel(c-col+h, r-row+h);
              V   += wt*img(c, r);
              sum += wt;
          }
        }
        if (sum > 0) filled_img(col, row) = V/sum;
        
      }
    } 
      
    return prerasterize_type(filled_img,
                             -biased_box.min().x(), -biased_box.min().y(),
                             cols(), rows());
      
  }
  template <class ImgT>
  inline void rasterize( ImgT const& img, BBox2i const& bbox ) const {
    vw::rasterize( prerasterize(bbox), img, bbox );
  }
  
};

template <class ImgT>
BlurDEM<ImgT>
blur_dem( ImageViewBase<ImgT> const& img, double blur_sigma) {
  typedef BlurDEM<ImgT> result_type;
  return result_type( img.impl(), blur_sigma);
}

template<class ImageT>
void calc_image_and_reflectance_at_pt(ImageView<PixelGray<float> > const& DEM,
                                      ImageT const& img,
                                      GeoReference const& DEMGeo,
                                      float DEMnodata,
                                      ModelParams const & img_params,
                                      GlobalParams const& globalParams,
                                      int col, int row, 
                                      IsisInterface const* cam,
                                      PixelMask<float> & imgval,
                                      PixelMask<PixelGray<float> > & refl
                                      ){

  refl.invalidate();
  imgval.invalidate();

  double h = DEM(col, row);
  if (h == DEMnodata) return;

  ImageViewRef<typename ImageT::pixel_type> interp_img
    = interpolate(img,
                  BicubicInterpolation(), ZeroEdgeExtension());
  
  Vector2 lonlat = DEMGeo.pixel_to_lonlat(Vector2(col, row));
  Vector3 llh(lonlat[0], lonlat[1], h);
  Vector3 xyz = DEMGeo.datum().geodetic_to_cartesian(llh);

  // Switch from VW vector to vector<double>
  vector<double> xyz2(3);
  for (int k = 0; k < (int)xyz.size(); k++) xyz2[k] = xyz[k];

  vector<double> cpix = cam->point_to_pixel(xyz2);
  
  if (cpix[0] < 1 || cpix[0] > interp_img.cols() - 2 ||
      cpix[1] < 1 || cpix[1] > interp_img.rows() - 2){
    imgval.invalidate();
  }else{
    imgval = interp_img(cpix[0], cpix[1]);
  }

  bool savePhaseAngle = false;
  float outputReflectanceVal, phaseAngleVal;
  computeReflectanceAtPixel(col, row, DEM, DEMGeo, DEMnodata, img_params,
                            globalParams, savePhaseAngle,
                            outputReflectanceVal,
                            phaseAngleVal);
  refl = outputReflectanceVal;
  if (refl.child() <= 0){
    //imgval.invalidate();
    refl.invalidate();
  }
  
}

void dump_iter(int iter,
               double a, double b,
               GeoReference const& DEMGeo,
               ImageView<PixelMask<float> > const& I,
               ImageView<PixelMask<PixelGray<float> > > const& R
               ){

  // The image intensity at every DEM point (I), and the reflectance
  // at every DEM point (R) are floating-point quantities.  Towards
  // the end of optimization we will have I is about a*R + b.

  // Here, we will clamp I to [mean(I) - sigma*stdev(I), mean(I) + sigma*stddev(I)],
  // and cast to int. We will clamp a*R + b to the same interval for easy comparison.

  double imgmean, imgstdev;
  compute_image_stats(I, imgmean, imgstdev);
  
  double sigma = 5.0, mn = imgmean - sigma*imgstdev, mx = imgmean + sigma*imgstdev;
  ostringstream os; os << "_iter" << iter << "_int.tif";
  float nodata = 0;
  ImageView<uint8> V(I.cols(), I.rows());

  for (int col = 0; col < V.cols(); col++){
    for (int row = 0; row < V.rows(); row++){
      int v;
      if (!is_valid(I(col, row)))
        v = 0;
      else
        v = (int)round(255*(I(col, row).child() - mn)/(mx - mn));
      v = std::max(0, v);
      v = std::min(255, v);
      V(col, row) = uint8(v);
    }
  }
  write_image_nodata("demimg" + os.str(), V, DEMGeo, nodata);

  for (int col = 0; col < V.cols(); col++){
    for (int row = 0; row < V.rows(); row++){
      int v;
      if (!is_valid(R(col, row)))
        v = 0;
      else
        v = (int)round(255*((a*R(col, row).child()+b) - mn)/(mx - mn));
      v = std::max(0, v);
      v = std::min(255, v);
      V(col, row) = uint8(v);
    }
  }
  write_image_nodata("demrefl" + os.str(), V, DEMGeo, nodata);
  
  
}

double calc_cost_fun(ImageView<PixelGray<float> >& DEM,
                     GeoReference const& DEMGeo,
                     float DEMnodata,
                     ModelParams const& img_params,
                     GlobalParams const& globalParams,
                     IsisInterface const* cam,
                     ImageView<float> & img,
                     double a, double b,
                     ImageView<PixelMask<float> > & I,
                     ImageView<PixelMask<PixelGray<float> > > & R
                     ){

  // Compute the cost function.
  
  I.set_size(DEM.cols(), DEM.rows());
  R.set_size(DEM.cols(), DEM.rows());
  for (int col = 0; col < DEM.cols(); col++){
    for (int row = 0; row < DEM.rows(); row++){
      R(col, row).invalidate();
      I(col, row).invalidate();
    }
  }
  
  double cost = 0.0;
  
  // We will keep the border pixels fixed
  for (int col = 1; col < DEM.cols()-1; col++){
    for (int row = 1; row < DEM.rows()-1; row++){
      
      if (DEM(col, row) == DEMnodata) continue;
      
      bool bad_value = false;
      double cost2 = 0.0;
      for (int c = col-1; c <= col+1; c++){
        for (int r = row-1; r <= row+1; r++){
          
          PixelMask<float> imgval;
          PixelMask<PixelGray<float> > refl;
          calc_image_and_reflectance_at_pt(DEM, img,  
                                           DEMGeo, DEMnodata,  
                                           img_params, globalParams,  
                                           c, r, cam,  
                                           imgval, refl
                                           );
          
          if (!is_valid(imgval) || !is_valid(refl)){
            bad_value = true;
            continue;
          }
          if (c == col && r == row){
            I(c, r) = imgval;
            R(c, r) = refl; //a*refl + b;
          }
          
          double v = imgval.child() - a*refl.child() - b;
          cost2 += v*v;
        }
      }
      if (!bad_value){
        cost += cost2;
      }
      
    }
  }

  return cost;
}

void calc_update(ImageView<PixelGray<float> >& DEM,
                 GeoReference const& DEMGeo,
                 float DEMnodata,
                 ModelParams const& img_params,
                 GlobalParams const& globalParams,
                 IsisInterface const* cam,
                 ImageView<float> & img,
                 double a, double b, double deltah,
                 ImageView<PixelGray<float> > & update
                 ){
  
  // Calculate the update via gradient descent (the gradient is
  // computed numerically).
  
  update.set_size(DEM.cols(), DEM.rows());
  for (int col = 0; col < DEM.cols(); col++){
    for (int row = 0; row < DEM.rows(); row++){
      update(col, row) = 0;
    }
  }
  
  // We will keep the border pixels fixed
  for (int col = 1; col < DEM.cols()-1; col++){
    for (int row = 1; row < DEM.rows()-1; row++){
      
      if (DEM(col, row) == DEMnodata) continue;
      
      bool bad_value = false;
      
      double cost = 0.0;
      for (int c = col-1; c <= col+1; c++){
        for (int r = row-1; r <= row+1; r++){
          
          PixelMask<float> imgval;
          PixelMask<PixelGray<float> > refl;
          calc_image_and_reflectance_at_pt(DEM, img,  
                                           DEMGeo, DEMnodata,  
                                           img_params, globalParams,  
                                           c, r, cam,  
                                           imgval, refl
                                           );
          
          if (!is_valid(imgval) || !is_valid(refl)){
            bad_value = true;
            continue;
          }
          double v = imgval.child() - a*refl.child() - b;
          cost += v*v;
        }
      }
      
      // Now bias the DEM at current location
      DEM(col, row) += deltah;
      
      double cost2 = 0.0;
      for (int c = col-1; c <= col+1; c++){
        for (int r = row-1; r <= row+1; r++){
          
          PixelMask<float> imgval;
          PixelMask<PixelGray<float> > refl;
          calc_image_and_reflectance_at_pt(DEM, img,  
                                           DEMGeo, DEMnodata,  
                                           img_params, globalParams,  
                                           c, r, cam,  
                                           imgval, refl
                                           );
          
          if (!is_valid(imgval) || !is_valid(refl)){
            bad_value = true;
            continue;
          }
          double v = imgval.child() - a*refl.child() - b;
          cost2 += v*v;
        }
      }
      
      // Must bias the DEM back
      DEM(col, row) -= deltah;
      
      if (!bad_value){
        update(col, row) = (cost2-cost)/deltah;
      }
      
    }
  }
  
}

void calc_image_and_reflectance(ImageView<PixelGray<float> >& DEM,
                                GeoReference const& DEMGeo,
                                float DEMnodata,
                                ModelParams const& img_params,
                                GlobalParams const& globalParams,
                                int num_iter
                                ){

  std::string imgfile = img_params.inputFilename;

  boost::shared_ptr<IsisInterface> cam(IsisInterface::open(imgfile));

  boost::shared_ptr<DiskImageResource>
    img_rsrc( DiskImageResource::open(imgfile));
  ImageView<float> img = copy(DiskImageView<float>( img_rsrc ));

  // Calculate a and b in the camera transfer function model
  // I = a*R + b.
  double a = 1.0, b = 0.0;
  ImageView<PixelMask<float> > I;
  ImageView<PixelMask<PixelGray<float> > > R;
  calc_cost_fun(DEM, DEMGeo, DEMnodata,  
                img_params,  globalParams, cam.get(), img, a, b,
                I, R
                );
  double imgmean, imgstdev, refmean, refstdev;
  compute_image_stats(I, imgmean, imgstdev);
  compute_image_stats(R, refmean, refstdev);
  a = imgstdev/refstdev;
  b = imgmean - a*refmean;
  
  // The cost fun is sum (I-a*R-b)^2. Each time the height of a DEM
  // point is changed, I and R change only at that point and its
  // immediate neighbors.

  // Will vary height by this much to compute the numerical gradient
  double max_dem = calc_max_abs(DEM, DEMnodata);
  double deltah = 0.00001*max_dem;
  std::cout << "deltah is " << deltah << std::endl;
    
  // Initialize the update
  ImageView<PixelGray<float> > update;
  calc_update(DEM, DEMGeo, DEMnodata,  img_params,  
              globalParams, cam.get(), img, a, b, deltah, update
              );

  // Determine how to scale the updates, so they are not too big or
  // too small.
  double max_up  = calc_max_abs(update, 0.0);
  double mn, mx;
  calc_min_max(DEM, DEMnodata, mn, mx);
  double pct = 0.05;
  double factor = pct*(mx - mn)/max_up;
  if (factor == 0){
    // Handle the case when the input DEM is constant
    factor = 2.0/max_up; // 2 meter update
  }
    
  double cost = calc_cost_fun(DEM, DEMGeo, DEMnodata,  
                              img_params,  globalParams, cam.get(), img,  a, b,
                              I, R
                              );

  vector<double> Cost;
  Cost.push_back(cost);
  ofstream cs("cost.txt");
  cs << cost << std::endl;

  dump_iter(0, a, b, DEMGeo, I, R);

  for (int i = 1; i <= num_iter; i++){
      
    std::cout << "iteration: " << i << std::endl;

    double sigma = 1.5;
    ImageView<PixelGray<float> > blurredDEM  
      = apply_mask(blur_dem
                   (create_mask(DEM, DEMnodata), sigma),
                   DEMnodata);
    DEM = copy(blurredDEM);
      
    // Cost after blur
    cost = calc_cost_fun(DEM, DEMGeo, DEMnodata,  
                         img_params,  globalParams, cam.get(), img, a, b,
                         I, R
                         );
    Cost.push_back(cost);
    cs << cost << std::endl;
      
    calc_update(DEM, DEMGeo, DEMnodata,  img_params,  
                globalParams, cam.get(), img, a, b, deltah, update
                );

    double max_up = 0.0;
    for (int col = 1; col < DEM.cols()-1; col++){
      for (int row = 1; row < DEM.rows()-1; row++){
        if (DEM(col, row) == DEMnodata) continue;
        double u = factor*update(col, row);
        DEM(col, row) -= u;
        max_up = std::max(max_up, u);
      }
    }

    std::cout << "max update is " << max_up << std::endl;
      
    // Cost after update
    cost = calc_cost_fun(DEM, DEMGeo, DEMnodata,  
                         img_params,  globalParams, cam.get(), img,  a, b,
                         I, R
                         );
    Cost.push_back(cost);
    cs << cost << std::endl;

    if (i%10 == 0) dump_iter(i, a, b, DEMGeo, I, R);
    
  }

  cs.close();

  dump_iter(num_iter, a, b, DEMGeo, I, R);
}

int main( int argc, char *argv[] ) {

  // Do shape-from shading. See the INSTALL file about how to compile
  // and run this program.

  // Important: The step size for gradient descent needs tweaking for
  // best performance.  Ideally it would be adaptive, growing or
  // shrinking depending on how the cost function is decreasing.
  
  // To do: Remove the unneeded initial blur
  
  if (argc < 5){
    std::cerr << "Usage: " << argv[0] << " image.cub DEM.tif metadata_dir num_iter" << std::endl;
    exit(1);
  }
  
  vector<string> images;
  images.push_back(argv[1]);
  
  std::string DEMFile = argv[2];
  std::string meta_dir = argv[3];
  int num_iter = atoi(argv[4]);
  
  ImageView<PixelGray<float> > DEMTile = copy(DiskImageView<PixelGray<float> >(DEMFile));
  GeoReference DEMGeo;
  read_georeference(DEMGeo, DEMFile);
  float DEMnodata = -32768;
  if (readNoDEMDataVal(DEMFile, DEMnodata)){
    std::cout << "Found DEM nodata value: " << DEMnodata << std::endl;
  }

  GlobalParams globalParams;
  globalParams.sunPosFile = meta_dir + "/sunpos.txt";
  globalParams.spacecraftPosFile = meta_dir + "/spacecraftpos.txt";
  globalParams.reflectanceType = LUNAR_LAMBERT;
  globalParams.phaseCoeffC1    = 1.383488;
  globalParams.phaseCoeffC2    = 0.501149;
  //PrintGlobalParams(globalParams);


  std::map<std::string, Vector3> sunPositions;
  std::map<std::string, Vector3> spacecraftPositions;
  ReadSunOrSpacecraftPosition(globalParams.sunPosFile, // Input
                              sunPositions             // Output
                              );
  ReadSunOrSpacecraftPosition(globalParams.spacecraftPosFile, // Input
                              spacecraftPositions             // Output
                              );
  
  std::vector<ModelParams> modelParamsArray;
  std::cout.precision(18);
  modelParamsArray.resize(images.size());
  for (int k = 0; k < (int)images.size(); k++){
    std::string prefix = getFirstElevenCharsFromFileName(images[k]);

    modelParamsArray[k].inputFilename = images[k];
      
    if ( sunPositions.find(prefix) == sunPositions.end()){
      std::cerr << "Could not find the sun position for the image file: "
                << images[k] << std::endl;
      exit(1);
    }
    // Go from kilometers to meters
    modelParamsArray[k].sunPosition = 1000*sunPositions[prefix];
    std::cout << "sun position: " <<  modelParamsArray[k].inputFilename
              << ' ' <<  modelParamsArray[k].sunPosition << std::endl;

    if (spacecraftPositions.find(prefix) == spacecraftPositions.end()){
      std::cerr << "Could not find the spacecraft position for the DRG file: "
                << images[k] << std::endl;
      exit(1);
    }
    // Go from kilometers to meters
    modelParamsArray[k].spacecraftPosition = 1000*spacecraftPositions[prefix];
    std::cout << "spacecraft position: " <<  modelParamsArray[k].inputFilename
              << ' ' <<  modelParamsArray[k].spacecraftPosition << std::endl;
  }
  
#if 0
  // If to start with a constant initial guess
  for (int col = 0; col < DEMTile.cols(); col++){
    for (int row = 0; row < DEMTile.rows(); row++){
      DEMTile(col, row) = mean;
    }
  }
#endif
  
  calc_image_and_reflectance(DEMTile, DEMGeo,  
                             DEMnodata,  
                             modelParamsArray[0],  
                             globalParams,
                             num_iter
                             );
  return 0;
}


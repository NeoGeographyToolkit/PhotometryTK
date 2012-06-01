// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Misc.h

#ifndef __VW_PHOTOMETRY_MISC_H__
#define __VW_PHOTOMETRY_MISC_H__

#include <sys/stat.h>
#include <math.h>
#include <time.h>
#include <vector>

#include <vw/Image/ImageViewBase.h>
#include <vw/Cartography/GeoReference.h>
#include <photk/Reconstruct.h>

namespace photometry {

  void upsample_uint8_image(std::string output_file, std::string input_file, int upsampleFactor);
  
  //upsamples a DEM by an upsampleFactor
  void upsample_image(std::string output_file, std::string input_file, int upsampleFactor);

  //subsamples a geo referenced tiff image by two
  void subsample_image(std::string output_file, std::string input_file);

  // Given two images and two georeferences, this function picks a set
  // of matching pixel samples between the two images.  It rejects
  // pixels that are not valid, and it should probably also reject
  // pixels that are near saturation (though it does not yet!).
  template <class ViewT>
  std::vector<vw::Vector4> sample_images(vw::ImageViewBase<ViewT> const& image1,
                                         vw::ImageViewBase<ViewT> const& image2,
                                         vw::cartography::GeoReference const& geo1,
                                         vw::cartography::GeoReference const& geo2,
                                         int num_samples,
                                         std::string const& DEM_file,
                                         std::vector<vw::Vector3> *normalArray,
                                         std::vector<vw::Vector3> *xyzArray );


  /// Erases a file suffix if one exists and returns the base string
  static std::string prefix_from_filename(std::string const& filename) {
    std::string result = filename;
    int index = result.rfind(".");
    if (index != -1)
      result.erase(index, result.size());
    return result;
  }

  /// Erases a file suffix if one exists and returns the base string less3 characters
  static std::string prefix_less3_from_filename(std::string const& filename) {
    std::string result = filename;
    int index = result.rfind(".");
    if (index != -1)
      result.erase(index-3, result.size()+3);
    return result;
  }

  /// Erases a file suffix if one exists and returns the base string less3 characters
  static std::string suffix_from_filename(std::string const& filename) {
    std::string result = filename;
    int index = result.rfind("/");
    if (index != -1)
      result.erase(0, index);
    return result;
  }

  //reads the tiff DEM into a 3D coordinate
  //pos is a vw::Vector2 of pixel coordinates, GR is georeference
  template <class ViewT>
  vw::Vector3 pixel_to_cart (vw::Vector2 pos, vw::ImageViewBase<ViewT> const& img,
                             vw::cartography::GeoReference GR);


  std::vector<std::string> parse_command_arguments(int argc, char *argv[] );

  void getTileCornersWithoutPadding(// Inputs
                                    int numCols, int numRows,
                                    vw::cartography::GeoReference const& geoRef,
                                    double tileSize, int pixelPadding,
                                    // Outputs
                                    double & min_x, double & max_x,
                                    double & min_y, double & max_y
                                    );
  
  void applyPaddingToTileCorners(// Inputs
                                 vw::cartography::GeoReference const& geoRef,
                                 int pixelPadding,
                                 double min_x, double max_x,
                                 double min_y, double max_y,
                                 // Outputs
                                 double & min_x_padded, double & max_x_padded,
                                 double & min_y_padded, double & max_y_padded);
    
  void readDEMTilesIntersectingBox(// Inputs
                                   double noDEMDataValue,
                                   vw::Vector2 boxNW, vw::Vector2 boxSE,
                                   std::vector<std::string> const& DEMTiles,
                                   // Outputs
                                   vw::ImageView<vw::PixelGray<float> > & combinedDEM,
                                   vw::cartography::GeoReference    & combinedDEM_geo);

  void listTifsInDir(const std::string & dirName,
                     std::vector<std::string> & tifsInDir
                     );

  void writeSunAndSpacecraftPosition(std::string prefix, std::string sunFile, std::string spacecraftFile,
                                     vw::Vector3 sunPosition, vw::Vector3 spacecraftPosition);
    
  std::string getFirstElevenCharsFromFileName(std::string fileName);
    
  void indexFilesByKey(std::string dirName, std::map<std::string, std::string> & index);

  void enforceUint8Img(std::string imgName);

  bool readNoDEMDataVal(std::string DEMFile, float & noDEMDataValue);
    
  void maskPixels(std::string imgFile, std::string maskFile, double shadowThresh, std::string outDir);

  void ReadPhaseCoeffsFromFile(std::string phaseDir, GlobalParams& settings);
  void AppendPhaseCoeffsToFile(const GlobalParams& settings);

  float getShadowThresh(const GlobalParams& settings, float exposureRefl);

  void resampleImage(std::string initFilename, std::string outputFilename, int factor);

  bool boxesOverlap(const vw::Vector4 & box1Corners, const vw::Vector4 & box2Corners);

  vw::Vector4 ComputeGeoBoundary(vw::cartography::GeoReference Geo, int width, int height);

  vw::Vector4 getImageCorners(std::string imageFile);

  void listTifsInDirOverlappingWithBox(const std::string & dirName,
                                       vw::Vector4 & boxCorners,
                                       const std::string & outputListName);
      
  void createAlbedoTilesOverlappingWithDRG(double tileSize, int pixelPadding,
                                           std::string imageFile, vw::Vector4 const& simulationBox,
                                           std::vector<ImageRecord> const& drgRecords,
                                           std::string blankTilesList,  std::string blankTilesDir,
                                           std::string DEMTilesList,    std::string meanDEMDir,
                                           std::string albedoTilesList, std::string albedoDir
                                           );

  std::vector<int> GetInputIndices( std::vector<std::string> inputFiles, std::vector<std::string> DRGFiles);
  std::vector<int> makeOverlapList(const std::vector<ModelParams>& drgFiles,
                                   const std::string& currFile);

  std::vector<int> makeOverlapList(const std::vector<ImageRecord>& drgRecords,
                                   const std::string& currFile);
      
  std::vector<std::vector<int> > makeOverlapList(const std::vector<std::string>& inputFiles,
                                                 const std::vector<ModelParams>& DRGFiles);
    
  void printOverlapList(std::vector<std::vector<int> > overlapIndices);

  vw::Vector4 parseSimBox(std::string simulationBoxStr);

  void extractSimBox(char * line, vw::Vector4 & simulationBox);

  int ReadConfigFile(char *config_filename, struct GlobalParams & settings);

  void PrintGlobalParams(GlobalParams& settings);

  bool readImagesFile(std::vector<ImageRecord>& images,
                      const std::string& imagesListName);

  void list_DRG_in_box_and_all_DEM(bool useTiles, bool useReflectance,
                                   std::string allDRGIndex, std::string allDEMIndex,
                                   vw::Vector4 simulationBox, 
                                   std::string DRGDir,  std::string DEMDir, 
                                   std::string DRGInBoxList
                                   );

    
  template <class pixelInType, class pixelOutType>
  bool getSubImageWithMargin(// Inputs
                             vw::Vector2 begLonLat, vw::Vector2 endLonLat,
                             std::string imageFile,
                             // Outputs
                             vw::ImageView<pixelOutType>   & subImage,
                             vw::cartography::GeoReference & sub_geo){
    
    // Read from disk only the portion of image whose pixel lon lat
    // values are in between begLonLat and endLonLat.  Adjust the
    // georeference accordingly. This is necessary to save on memory
    // for large images.
    
    // IMPORTANT: We assume that the input image is of pixelInType, and
    // we will convert the portion we read into pixelOutType.  Also we
    // read a few more pixels on each side, to help later with
    // interpolation.
  
    int extra = 2;
    vw::DiskImageView<pixelInType> input_img(imageFile);
    //std::cout << "Reading: " << imageFile << std::endl;
    
    vw::cartography::GeoReference input_geo;
    read_georeference(input_geo, imageFile);

    vw::Vector2 begPix = input_geo.lonlat_to_pixel(begLonLat);
    vw::Vector2 endPix = input_geo.lonlat_to_pixel(endLonLat);

    int beg_col = std::max(0,                (int)floor(begPix(0)) - extra);
    int end_col = std::min(input_img.cols(), (int)ceil(endPix(0))  + extra);
    int beg_row = std::max(0,                (int)floor(begPix(1)) - extra);
    int end_row = std::min(input_img.rows(), (int)ceil(endPix(1))  + extra);
  
    if (beg_col >= end_col || beg_row >= end_row) return false;
  
    subImage.set_size(end_col - beg_col, end_row - beg_row);
    for (int row = beg_row; row < end_row; row++){
      for (int col = beg_col; col < end_col; col++){
        // Note the pixelInType to pixelOutType conversion below
        subImage(col - beg_col, row - beg_row) = input_img(col, row); 
      }
    }
  
    // In sub_geo, the (0, 0) pixel will where (beg_col, beg_row) is in input_geo.
    sub_geo = vw::cartography::crop(input_geo, beg_col, beg_row);

    //system("echo getSubImageWithMargin top is $(top -u $(whoami) -b -n 1|grep lt-reconstruct)");
  
    return true;
  }

  template<class T>
  void dumpImageToFile(T & img, std::string fileName){
    
    printf("dumping %s\n", fileName.c_str());
    
    std::ofstream fs(fileName.c_str());
    fs.precision(20);
    
    for (int k = 0 ; k < img.rows(); ++k) {
      for (int l = 0; l < img.cols(); ++l) {
        fs << (double)img(l, k) << " ";
      }
      fs << std::endl;
    }
    fs.close();
  }

    
} // end photometry

namespace photometry {

  template <class ImageT>
  bool isGoodPixel(ImageT const& maskImg, vw::cartography::GeoReference maskGeo, 
                   vw::Vector2 lonlat, double t){
  
    vw::Vector2 maskPix = maskGeo.lonlat_to_pixel(lonlat);
    float x = maskPix(0), y = maskPix(1);
  
    int x0 = (int)floor(x), x1 = (int)ceil(x);
    int y0 = (int)floor(y), y1 = (int)ceil(y);
    bool isGood = ((x >= 0) && (x <= maskImg.cols()-1) &&
                   (y >= 0) && (y <= maskImg.rows()-1) && 
                   is_valid(maskImg(x0, y0)) && (float)maskImg(x0, y0) >= t &&
                   is_valid(maskImg(x0, y1)) && (float)maskImg(x0, y1) >= t &&
                   is_valid(maskImg(x1, y0)) && (float)maskImg(x1, y0) >= t &&
                   is_valid(maskImg(x1, y1)) && (float)maskImg(x1, y1) >= t
                   );
  
    return isGood;
  }

  // Mask the pixels in a given image using pixels from a mask image.
  // All the process is lazy, the full images are never stored in memory.

  template <class ImageT>
  class MaskImage : public vw::ImageViewBase<MaskImage<ImageT> > {

    typedef typename boost::mpl::if_<vw::IsFloatingPointIndexable<ImageT>, double, vw::int32>::type offset_type;

    ImageT m_inputImg;
    const vw::ImageViewRef<vw::PixelMask<vw::PixelGray<vw::uint8> > > & m_maskImg;
    double m_shadowThresh;
    vw::cartography::GeoReference m_imgGeo;
    vw::cartography::GeoReference m_maskGeo;
    
  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef const pixel_type result_type;
    typedef vw::ProceduralPixelAccessor<MaskImage> pixel_accessor;

    MaskImage(vw::ImageViewRef<vw::PixelMask<vw::PixelGray<vw::uint8> > > const& inputImg,
              ImageT const& inputImg_prerasterized,
              vw::ImageViewRef<vw::PixelMask<vw::PixelGray<vw::uint8> > > const& maskImg,
              double shadowThresh, 
              vw::cartography::GeoReference const& imgGeo,
              vw::cartography::GeoReference const& maskGeo):
      m_inputImg(inputImg_prerasterized),
      m_maskImg(maskImg),
      m_shadowThresh(shadowThresh),
      m_imgGeo(imgGeo),
      m_maskGeo(maskGeo){}
    
    inline vw::int32 cols() const { return m_inputImg.cols(); }
    inline vw::int32 rows() const { return m_inputImg.rows(); }
    inline vw::int32 planes() const { return m_inputImg.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(*this); }

    inline result_type operator()( offset_type i, offset_type j, vw::int32 p=0 ) const {

      if  ((float)m_inputImg(i, j) == 0) result_type(0.0);
      
      // Note: Below we take into account that the lonlat of the image
      // may be shifted by 360 degrees from the lonlat of the mask.
      // This is a bug fix.
      vw::Vector2 lonlat = m_imgGeo.pixel_to_lonlat(vw::Vector2(i, j));
      double t = m_shadowThresh;
      if (!isGoodPixel(m_maskImg, m_maskGeo, lonlat, t)                   &&
          !isGoodPixel(m_maskImg, m_maskGeo, lonlat + vw::Vector2(360, 0), t) && 
          !isGoodPixel(m_maskImg, m_maskGeo, lonlat - vw::Vector2(360, 0), t)
          ){
        return result_type(0.0);
      }else{
        return m_inputImg(i, j);
      }
      
    }
    
    /// \cond INTERNAL
    typedef MaskImage<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( vw::BBox2i const& bbox ) const {
      return prerasterize_type(m_inputImg, m_inputImg.prerasterize(bbox), m_maskImg, m_shadowThresh, 
                               m_imgGeo, m_maskGeo);
    }
    template <class DestT> inline void rasterize( DestT const& dest, vw::BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  // --------------------------------------------------------------------------
  // Functional API. See the description of the MaskImage for more info.
  // --------------------------------------------------------------------------
  template <class ImageT>
  MaskImage<ImageT>
  mask_image(vw::ImageViewBase<ImageT> const& inputImg,
             vw::ImageViewBase<ImageT> const& maskImg,
             double shadowThresh,
             vw::cartography::GeoReference const& imgGeo,
             vw::cartography::GeoReference const& maskGeo) {
    return MaskImage<ImageT>(inputImg.impl(), inputImg.impl(), maskImg.impl(), shadowThresh, imgGeo, maskGeo);
  }

} // namespace photometry

namespace vw{
  /// \cond INTERNAL
  // Type traits
  template <class ImageT>
  struct IsFloatingPointIndexable< photometry::MaskImage<ImageT> > : public IsFloatingPointIndexable<ImageT> {};
  /// \endcond
}

#endif//__VW_PHOTOMETRY_MISC_H__

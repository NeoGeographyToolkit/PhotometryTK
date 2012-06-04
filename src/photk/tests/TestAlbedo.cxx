// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <string>
#include <vector>
#include <gtest/gtest.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <photk/Photometry.h>
#include <test/Helpers.h>

using namespace vw;
using namespace photometry;
using namespace std;

bool compareDirsWithTol(string dir1, string dir2){

  // Given two output directories of running albedo, see if the albedo tiles
  // in those directories are similar enough
  
  int tol = 2; // image values go from 0 to 255, use 2 as the tolerance
  
  vector<string> tifsInDir;
  listTifsInDir(dir1 + "/albedo", tifsInDir);
  for (int s = 0; s < tifsInDir.size(); s++){
    string file1 = tifsInDir[s];
    string file2 = dir2 + "/albedo" + suffix_from_filename(file1);
    
    std::cout << "Comparing: " << file1  << ' ' << file2 << std::endl;

    ifstream h1(file1.c_str());
    if (!h1) {
      std::cerr << "Missing file: " << file1 << std::endl;
      return false;
    }
    ifstream h2(file2.c_str());
    if (!h2) {
      std::cerr << "Missing file: " << file2 << std::endl;
      return false;
    }

    DiskImageView<PixelMask<PixelGray<uint8> > >  img1(file1);
    DiskImageView<PixelMask<PixelGray<uint8> > >  img2(file2);
    if (img1.rows() != img2.rows() || img1.cols() != img2.cols()){
      std::cout << "The following images have different sizes: " << file1  << ' ' << file2 << std::endl;
      return false;
    }

    for (int k = 0 ; k < img1.rows(); ++k) {
      for (int l = 0; l < img1.cols(); ++l) {
        if ( abs( (double)img1(l, k) - (double)img2(l, k) ) > tol ) return false;
      }
    }
    
  }

  return true;
}

class TestAlbedo : public ::testing::Test {
protected:
  TestAlbedo() {}
  virtual void SetUp() {}
};

TEST_F( TestAlbedo, FullCycle ) {

  // Make symbolic links in the current directory to the script and data necessary to run this test.
  string files[] = {
    "cmpdirs.sh", "meta",
    "DIM_input_1280mpp_masked", "DIM_input_2560mpp",
    "DEM_tiles_sub64",
    "photometry_settings_1.txt", "albedo_gold_1",
    "photometry_settings_2.txt", "albedo_gold_2"
  };
  
  for (int s = 0 ; s < sizeof(files)/sizeof(string); s++){
    string cmd = "ln -s ../../../../src/photk/tests/" + files[s] + " .";
    std::cout << cmd << std::endl;
    system(cmd.c_str());
  }
  string cmd = "ln -s ../../../../src/tools/reconstruct.sh .";
  std::cout << cmd << std::endl;
  int flag = system(cmd.c_str());
  
  // Run test 1
  cmd="./reconstruct.sh photometry_settings_1.txt curr_1 > output1.txt";
  std::cout << cmd << std::endl;
  system(cmd.c_str());
  flag = compareDirsWithTol("albedo_gold_1", "albedo_curr_1");
  EXPECT_EQ(flag, 1);
  
  // Run test 2
  cmd="./reconstruct.sh photometry_settings_2.txt curr_2 > output2.txt";
  std::cout << cmd << std::endl;
  system(cmd.c_str());
  flag = compareDirsWithTol("albedo_gold_2", "albedo_curr_2");
  EXPECT_EQ(flag, 1);
}


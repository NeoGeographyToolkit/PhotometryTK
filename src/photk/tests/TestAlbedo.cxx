// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <string>
#include <gtest/gtest.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <photk/Photometry.h>
#include <test/Helpers.h>

using namespace vw;
using namespace photometry;
using namespace std;

template <class ChannelT>
class FullIterationTest : public ::testing::Test {
protected:
  FullIterationTest() {}
  virtual void SetUp() {}
};

typedef FullIterationTest<uint8> FullU8IterationTest;

TEST_F( FullU8IterationTest, FullCycle ) {

  system("pwd");
  
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
  // Run a directory comparison and check the status
  cmd="cmpdirs.sh albedo_gold_1 albedo_curr_1";
  std::cout << cmd << std::endl;
  flag = system(cmd.c_str());
  EXPECT_EQ(flag, 0);
  
  // Run test 2
  cmd="./reconstruct.sh photometry_settings_2.txt curr_2 > output2.txt";
  std::cout << cmd << std::endl;
  system(cmd.c_str());
  // Run a directory comparison and check the status
  cmd="cmpdirs.sh albedo_gold_2 albedo_curr_2";
  std::cout << cmd << std::endl;
  flag = system(cmd.c_str());
  EXPECT_EQ(flag, 0);

}


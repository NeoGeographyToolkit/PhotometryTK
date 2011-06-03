// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Macros.h
///

#ifndef __PHOTK_TOOLS_MACROS_H__
#define __PHOTK_TOOLS_MACROS_H__

#define PHOTK_STANDARD_CATCHES                              \
  catch ( const vw::ArgumentErr& e ) {                      \
    vw_out() << e.what() << std::endl;                      \
    return 1;                                               \
  } catch ( const vw::Exception& e ) {                      \
    std::cerr << "\n\nVW Error: " << e.what() << std::endl; \
    return 1;                                               \
  } catch ( const std::bad_alloc& e ) {                     \
    std::cerr << "\n\nError: Ran out of Memory!" << std::endl; \
    return 1;                                               \
  } catch ( const std::exception& e ) {                     \
    std::cerr << "\n\nError: " << e.what() <<  std::endl;   \
    return 1;                                               \
  }

#endif//__PHOTK_TOOLS_MACROS_H__

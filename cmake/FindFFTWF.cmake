# - Find FFTWF
# Find the native FFTWF includes and library
#
#  FFTWF_INCLUDES    - where to find fftw3.h
#  FFTWF_LIBRARIES   - List of libraries when using FFTWF.
#  FFTWF_FOUND       - True if FFTWF found.

if (FFTWF_INCLUDES)
  # Already in cache, be silent
  set (FFTWF_FIND_QUIETLY TRUE)
endif (FFTWF_INCLUDES)

find_path (FFTWF_INCLUDES fftw3.h)

find_library (FFTWF_LIBRARIES NAMES fftw3f)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTWF DEFAULT_MSG FFTWF_LIBRARIES FFTWF_INCLUDES)

mark_as_advanced (FFTWF_LIBRARIES FFTWF_INCLUDES)


program testing_ground
    implicit none
#ifdef _WIN32
      print *,'Windows'
#else
      print *,'Linux'
#endif
end program testing_ground

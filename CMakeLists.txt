INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
)

LINK_DIRECTORIES(
  ${PROJECT_BINARY_DIR}/lemon
)

# Uncomment (and adjust) the following two lines. 'myprog' is the name
# of the final executable ('.exe' will automatically be added to the
# name on Windows) and 'myprog-main.cc' is the source code it is
# compiled from. You can add more source files separated by
# whitespaces. Moreover, you can add multiple similar blocks if you
# want to build more than one executables.

ADD_EXECUTABLE(robust_energy Robust_Network_Energy/main.cpp)
#target_sources(robust Robust_Network_Energy/robust_energy.cpp )
#TARGET_LINK_LIBRARIES(energy lemon)
TARGET_LINK_LIBRARIES(robust_energy lemon)

find_package(Boost 1.48.0 QUIET COMPONENTS unit_test_framework timer)

if(NOT Boost_FOUND)
	message("Boost libraries not found. C++ tests will not be build.")
	return()
endif(NOT Boost_FOUND)



# use C++14 because of a 'auto' parameter in a lambda
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wall -std=c++14 -pedantic")

add_definitions(-DBOOST_TEST_DYN_LINK)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR}) # for config.hpp
include_directories(${Boost_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/LAPACKE/include)


configure_file(config.hpp.in config.hpp @ONLY)


add_executable(tools_tests main.cpp tools_tests.cpp)
target_link_libraries(
	tools_tests tmglib ${LAPACK_LIBRARIES} ${Boost_LIBRARIES}
)
add_test(NAME tools_tests COMMAND tools_tests)


add_executable(xLASRTI_tests main.cpp xLASRTI_tests.cpp)
target_link_libraries(
	xLASRTI_tests tmglib ${LAPACK_LIBRARIES} ${Boost_LIBRARIES}
)
add_test(NAME xLASRTI_tests COMMAND xLASRTI_tests)


add_executable(xLASRTR_tests main.cpp xLASRTR_tests.cpp)
target_link_libraries(
	xLASRTR_tests tmglib ${LAPACK_LIBRARIES} ${Boost_LIBRARIES}
)
add_test(NAME xLASRTR_tests COMMAND xLASRTR_tests)


add_executable(xGGQRCS_tests main.cpp xGGQRCS_tests.cpp)
target_link_libraries(
	xGGQRCS_tests tmglib ${LAPACK_LIBRARIES} ${Boost_LIBRARIES}
)
add_test(NAME xGGQRCS_tests COMMAND xGGQRCS_tests)


add_executable(xUNCSD2BY1_tests main.cpp xUNCSD2BY1_tests.cpp)
target_link_libraries(
	xUNCSD2BY1_tests tmglib ${LAPACK_LIBRARIES} ${Boost_LIBRARIES}
)
add_test(NAME xUNCSD2BY1_tests COMMAND xUNCSD2BY1_tests)

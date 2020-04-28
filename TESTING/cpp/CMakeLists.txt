find_package(Boost 1.48.0 QUIET COMPONENTS unit_test_framework timer)

if(NOT Boost_FOUND)
	message("Boost libraries not found. C++ tests will not be build.")
	return()
endif(NOT Boost_FOUND)



# use C++14 because of a 'auto' parameter in a lambda
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wall -std=c++14 -pedantic")

add_definitions(-DBOOST_TEST_DYN_LINK)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/LAPACKE/include)


if(BUILD_SINGLE)
	add_executable(SGGQRCS_tests main.cpp sggqrcs.cpp)
	target_link_libraries(
		SGGQRCS_tests tmglib ${LAPACK_LIBRARIES} ${Boost_LIBRARIES}
	)
	add_test(NAME SGGQRCS_tests COMMAND SGGQRCS_tests)

	add_executable(SLASRTI_tests main.cpp slasrti.cpp)
	target_link_libraries(
		SLASRTI_tests tmglib ${LAPACK_LIBRARIES} ${Boost_LIBRARIES}
	)
	add_test(NAME SLASRTI_tests COMMAND SLASRTI_tests)
endif()

if(BUILD_DOUBLE)
	add_executable(DGGQRCS_tests main.cpp dggqrcs.cpp)
	target_link_libraries(
		DGGQRCS_tests tmglib ${LAPACK_LIBRARIES} ${Boost_LIBRARIES}
	)
	add_test(NAME DGGQRCS_tests COMMAND DGGQRCS_tests)
endif()

if(BUILD_COMPLEX)
	add_executable(CGGQRCS_tests main.cpp cggqrcs.cpp)
	target_link_libraries(
		CGGQRCS_tests tmglib ${LAPACK_LIBRARIES} ${Boost_LIBRARIES}
	)
	add_test(NAME CGGQRCS_tests COMMAND CGGQRCS_tests)
endif()

if(BUILD_COMPLEX16)
	add_executable(ZGGQRCS_tests main.cpp zggqrcs.cpp)
	target_link_libraries(
		ZGGQRCS_tests tmglib ${LAPACK_LIBRARIES} ${Boost_LIBRARIES}
	)
	add_test(NAME ZGGQRCS_tests COMMAND ZGGQRCS_tests)
endif()

find_package(Boost 1.48.0 QUIET COMPONENTS unit_test_framework)

if(NOT Boost_FOUND)
	message("Boost libraries not found. C++ tests will not be build.")
	return()
endif(NOT Boost_FOUND)


configure_file(config.hpp.in config.hpp @ONLY)


# use C++14 instead of C++11 because of a 'auto' parameter in a lambda
add_compile_options(-Wextra -Wall -std=c++14 -pedantic)


add_library(test_main OBJECT main.cpp)
target_compile_definitions(test_main PUBLIC -DBOOST_TEST_DYN_LINK)
# target_include_directories does not seem to do anything
include_directories(${Boost_INCLUDE_DIR})



function(add_cxx_test target source_file)
	add_executable(
		${target}
		$<TARGET_OBJECTS:test_main>
		${source_file}
	)

	target_include_directories(
		${target}
		PRIVATE ${Boost_INCLUDE_DIR}
		PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
		PRIVATE ${CMAKE_CURRENT_BINARY_DIR} # for config.hpp
		PRIVATE ${CMAKE_SOURCE_DIR}/LAPACKE/include
	)
	target_link_libraries(
		${target}
		PRIVATE tmglib ${LAPACK_LIBRARIES} ${Boost_LIBRARIES}
	)

	add_test(NAME ${target} COMMAND ${target})
endfunction()


add_cxx_test(tools_tests tools_tests.cpp)
add_cxx_test(xGGQRCS_tests xGGQRCS_tests.cpp)
add_cxx_test(xLASRTI_tests xLASRTI_tests.cpp)
add_cxx_test(xLASRTR_tests xLASRTR_tests.cpp)
add_cxx_test(xUNCSD2BY1_tests xUNCSD2BY1_tests.cpp)


#
# Benchmarks
#

add_executable(gsvd_benchmarks gsvd_benchmarks.cpp)
target_include_directories(
	gsvd_benchmarks
	PRIVATE ${Boost_INCLUDE_DIR}
	PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
	PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
	PRIVATE ${CMAKE_SOURCE_DIR}/LAPACKE/include
)
target_link_libraries(
	gsvd_benchmarks
	PRIVATE ${Boost_LIBRARIES}
	PRIVATE tmglib
	PRIVATE ${LAPACK_LIBRARIES}
)

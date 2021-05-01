// Copyright 2021 Christoph Conrads

// This code implements several matrix transposition approaches and compares
// their performance.

#include <cassert>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <chrono>
#include <limits>
#include <random>
#include <vector>


struct CpuClock {
	using duration = std::chrono::duration<std::intmax_t, std::nano>;
	using time_point = std::chrono::time_point<CpuClock>;

	static time_point now() {
		auto tp = timespec{-1, -1};
		auto ret = clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tp);

		if(ret != 0) {
			std::perror("clock_gettime");
			std::_Exit(EXIT_FAILURE);
		}

		auto second2nanosecond = std::intmax_t{1000000000};
		auto t =
			second2nanosecond * std::intmax_t{tp.tv_sec}
			+ std::intmax_t{tp.tv_nsec}
		;

		return time_point(duration(t));
	}
};

template<typename Number>
void transpose_simple(
	std::size_t m, std::size_t n,
	const Number* a, std::size_t lda,
	Number* b, std::size_t ldb)
	__attribute__((noinline));

template<typename Number>
void transpose_simple(
	std::size_t m, std::size_t n,
	const Number* a, std::size_t lda,
	Number* b, std::size_t ldb) {
	assert(lda >= m);
	assert(ldb >= n);

	for(auto j = std::size_t{0}; j < n; ++j) {
		for(auto i = std::size_t{0}; i < m; ++i) {
			b[i * ldb + j] = a[j * lda + i];
		}
	}
}


template<typename Number>
void transpose_blocked(
	std::size_t m, std::size_t n,
	const Number* a, std::size_t lda,
	Number* b, std::size_t ldb)
	__attribute__((noinline));

template<typename Number>
void transpose_blocked(
	std::size_t m, std::size_t n,
	const Number* u, std::size_t ldu,
	Number* v, std::size_t ldv) {
	assert(ldu >= m);
	assert(ldv >= n);

	constexpr auto B = std::size_t{64} / sizeof(Number);

	for(auto j = std::size_t{0}; j < n + (B - 1); j += B) {
		for(auto i = std::size_t{0}; i < m + (B - 1); i += B) {
			for(auto l = j; l < std::min(j + B, n); ++l) {
				for(auto k = i; k < std::min(i + B, m); ++k) {
					v[k * ldv + l] = u[l * ldu + k];
				}
			}
		}
	}
}




int main() {
	using Number = float;
	using Real = float;

	auto m = std::size_t{27};
	auto n = std::size_t{69};

	constexpr auto nan = std::numeric_limits<Real>::quiet_NaN();
	auto lda = m;
	auto a = std::vector<Number>(lda * n, nan);
	auto ldb = n;
	auto b = std::vector<Number>(ldb * m, nan);
	auto gen = std::mt19937_64();
	auto dist = std::uniform_real_distribution<Real>(-1024, +1024);

	std::generate(a.begin(), a.end(), [&gen, &dist] () { return dist(gen); });

	auto num_iterations = 1u << 16;
	auto t0 = CpuClock::now();
	for(auto it = std::size_t{0}; it < num_iterations; ++it) {
		__builtin___clear_cache(a.data(), a.data() + a.size());
		__builtin___clear_cache(b.data(), b.data() + b.size());
		transpose_simple(m, n, a.data(), lda, b.data(), ldb);
	}
	auto t1 = CpuClock::now();
	auto t_iter = (t1 - t0) / num_iterations;

	for(auto j = std::size_t{0}; j < n; ++j) {
		for(auto i = std::size_t{0}; i < m; ++i) {
			assert(a[j * lda + i] == b[i * ldb + j]);
		}
	}

	std::printf("CPU time (sec): %8.2e\n", t_iter.count() / 1e9);
	auto sum = std::accumulate(b.cbegin(), b.cend(), Number{0});
	std::printf("dummy: %8.2e\n", sum);
}

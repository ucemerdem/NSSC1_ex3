// g++ main.cpp -std=c++17 -O3 -march=native -ffast-math -o solver_ref
// ./solver_ref reference 10 10000

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include <omp.h>

/// for mode selection for the makefile
#ifndef SCHEDULE_MODE
#define SCHEDULE_MODE static //deafult is 'static'
#endif
/// for output in the terminal
#ifndef OUTPUT_MODE
#define OUTPUT_MODE 0 
#endif

namespace program_options {

struct Options {
  std::string name;
  size_t N;
  size_t iters;
  void print() const {
    std::cout << "######## START JACOBI SOLVER ########" << std::endl << std::endl;
    std::printf("    name: %s\n", name.c_str());
    std::printf("       N: %zu\n", N);
    std::printf("   iters: %zu\n", iters);
    std::cout << std::endl << std::endl;
  }
};

auto parse(int argc, char *argv[]) {
  if (argc != 4)
    throw std::runtime_error("unexpected number of arguments");
  Options opts;
  opts.name = argv[1];
  if (std::sscanf(argv[2], "%zu", &opts.N) != 1 && opts.N >= 2)
    throw std::runtime_error("invalid parameter for N");
  if (std::sscanf(argv[3], "%zu", &opts.iters) != 1 && opts.iters != 0)
    throw std::runtime_error("invalid parameter for iters");
  return opts;
}

} // namespace program_options

int main(int argc, char *argv[]) try {

  // parse args
  auto opts = program_options::parse(argc, argv);
  opts.print();

  // initial guess (0.0) with fixed values in west (-100) and east (100)
  auto init = [N = opts.N, W = -100.0, E = 100.0]() -> auto {
    std::vector<double> res(N * N);
    for (size_t j = 0; j < N; ++j)
      for (size_t i = 0; i < N; ++i) {
        res[i + j * N] = 0.0;
        if (i % N == 0)
          res[i + j * N] = W;
        if (i % N == N - 1)
          res[i + j * N] = E;
      }
    return res;
  };

  // solver update
  auto jacobi_iter = [N = opts.N](const auto &xold, auto &xnew,
                                  bool residual = false) {
    auto h = 1.0 / (N - 1);
    auto h2 = h * h;

    int threads = 0;
    #pragma omp parallel
      #pragma omp single  //only for getting thread number for later output
      {
        threads = omp_get_num_threads();
      }

    ///////// PARALELLIZE HERE ///////////
    #pragma omp parallel for schedule(SCHEDULE_MODE)
    // all interior points
      for (size_t j = 1; j < N - 1; ++j) {
        // if(residual == false) {
        //   std::cout << "Hello from Thread " << omp_get_thread_num() << std::endl;
        // }
        for (size_t i = 1; i < N - 1; ++i) {
          auto w = xold[(i - 1) + (j)*N];
          auto e = xold[(i + 1) + (j)*N];
          auto n = xold[(i) + (j + 1) * N];
          auto s = xold[(i) + (j - 1) * N];
          auto c = xold[(i) + (j)*N];
          if (!residual)
            xnew[i + j * N] = (- (-1.0 / h2) * (w + e + n + s)) * h2 / 4.0;
          else
            xnew[i + j * N] = (-1.0 / h2) * (w + e + n + s - 4.0 * c);
        }
      }
      // isolating south boundary
      {
        size_t j = 0;
        for (size_t i = 1; i < N - 1; ++i) {
          auto w = xold[(i - 1) + (j)*N];
          auto e = xold[(i + 1) + (j)*N];
          auto n = xold[(i) + (j + 1) * N];
          auto s = n;
          auto c = xold[(i) + (j)*N];
          if (!residual)
            xnew[i + j * N] = (- (-1.0 / h2) * (w + e + n + s)) * h2 / 4.0;
          else
            xnew[i + j * N] = (-1.0 / h2) * (w + e + n + s - 4 * c);
        }
      }
      // isolating north boundary
      {
        size_t j = N - 1;
        for (size_t i = 1; i < N - 1; ++i) {
          auto w = xold[(i - 1) + (j)*N];
          auto e = xold[(i + 1) + (j)*N];
          auto s = xold[(i) + (j - 1) * N];
          auto n = s;
          auto c = xold[(i) + (j)*N];
          if (!residual)
            xnew[i + j * N] = (- (-1.0 / h2) * (w + e + n + s)) * h2 / 4.0;
          else
            xnew[i + j * N] = (-1.0 / h2) * (w + e + n + s - 4 * c);
        }
      }
      return threads;
  };

  // write vector to csv
  auto write = [N = opts.N, name = opts.name](const auto &x) -> auto {
    std::ofstream csv;
    csv.open(name + ".csv");
    for (size_t j = 0; j < N; ++j) {
      for (size_t i = 0; i < N - 1; ++i) {
        csv << x[i + j * N] << " ";
      }
      csv << x[(N - 1) + j * N];
      csv << "\n";
    }
    csv.close();
  };

  // 2 norm
  auto norm2 = [N = opts.N](const auto &vec) -> auto {
    double sum = 0.0;
    for (size_t j = 0; j < N; ++j)
      for (size_t i = 1; i < (N - 1); ++i)
        sum += vec[i + j * N] * vec[i + j * N];

    return std::sqrt(sum);
  };

  // Inf norm
  auto normInf = [N = opts.N](const auto &vec) -> auto {
    double max = 0.0;
    for (size_t j = 0; j < N; ++j)
      for (size_t i = 1; i < (N - 1); ++i)
        max = std::fabs(vec[i + j * N]) > max ? std::fabs(vec[i + j * N]) : max;
    return max;
  };

  //////// MAIN-MAIN --> start time here ////////
  double start = omp_get_wtime();

  auto x1 = init();
  auto x2 = x1;
  int threads = 0;
  for (size_t iter = 0; iter <= opts.iters; ++iter) {
    threads = jacobi_iter(x1, x2);
    std::swap(x1, x2);
  }
  double end = omp_get_wtime();
  double duration = end - start;
  ///////////////////////////////////////////////

  // write(b);

  //write(x2);
  jacobi_iter(x1, x2, true);
  int mode = OUTPUT_MODE;
  std::cout << "  norm2 = " << norm2(x2) << std::endl;
  std::cout << "normInf = " << normInf(x2) << std::endl << std::endl;
  
  if (mode == 0) {
    std::cout << "Using scheduling mode: static" << std::endl;
  } else if (mode == 1) {
    std::cout << "Using scheduling mode: static,1" << std::endl;
  } else if (mode == 2) {
    std::cout << "Using scheduling mode: dynamic" << std::endl;
  }
  std::cout << "Threads: " << threads << std::endl;
  std::cout << "--> DURATION: " << duration << std::endl << std::endl;
  std::cout << "#####################################" << std::endl;

  return EXIT_SUCCESS;
  
} catch (std::exception &e) {
  std::cout << e.what() << std::endl;
  return EXIT_FAILURE;
}

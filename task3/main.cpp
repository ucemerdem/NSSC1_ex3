#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <random>   // std::mt19937 (mersenne twister)
#include <omp.h>

// forward declaration
double process_MCI(int seed, double (*func)(double), double xmin, double xmax, int N);


//////// master function ////////////

auto master_function(int value, double (*func)(double), double xmin, double xmax, int N){
    std::mt19937 rng(value); // first RNG
    std::cout << std::endl << "##### START MONTE CARLO INTEGRATION #####" << std::endl << std::endl;
    

    double final_approx = 0.0;
    int threads_num = 0;

    #pragma omp parallel reduction(+:final_approx, threads_num)
    {
        threads_num += 1;   // for getting total thread number

        #pragma omp critical 
        {
            // for each thread: own random number
            int seed = rng();
            
            // call MCI function to calc for each thread
            double thread_approx = process_MCI(seed, func, xmin, xmax, N);
            std::cout << "-- Approx from Thread #" << omp_get_thread_num() << ": " << thread_approx << std::endl << std::endl;
            final_approx += thread_approx;
        }
        
    }
    final_approx = final_approx/double(threads_num);
    std::cout << "Threads used: " << threads_num  << std::endl;
    std::cout << "Total samples: " << N << std::endl << std::endl;
    std::cout << "FINAL APPROXIMATION: " << final_approx << std::endl;
            
};

////// math functions /////

double SINX(double x) {
    return sin(x);
}

double COS2XINV(double x) {
    return cos(1/x) * cos(1/x);
}

double X4M5(double x) {
    return 5 * x * x * x * x;
}

//// monte carlo integral ////

double process_MCI(int seed, double (*func)(double), double xmin, double xmax, int N) {
    // define rng and range
    std::mt19937 rng2(seed);
    std::uniform_real_distribution<double> x_i(static_cast<double>(xmin), static_cast<double>(xmax));

    double prae_factor = (xmax - xmin) / N;
    double sum = 0;
    
    // for each N: pick random value & calc integral
    for(int c = 0; c < N; c++) {
        double value = x_i(rng2);
        //std::cout << "--- loop " << c << ": " << value << std::endl;
        double I = prae_factor * func(value);
        sum += I;
    } 
    
    // return sum: whole integral
    return sum;


}; 


/////////// MAIN ///////////////

// use: 'export OMP_NUM_THREADS=number' (1, 5, ...) before executing!


int main(int argc, char* argv[]) {
    
    if (argc != 5) {
        std::cerr << "Error: Number of arguments incorrect!\n";
        return 1;
    }

    //// parsing //// 
    std::string func_name = argv[1];       
    double xmin = std::stod(argv[2]); 
    double xmax = std::stod(argv[3]); 
    int samples = std::stod(argv[4]); 

    double (*func)(double) = nullptr;
    if (func_name == "SINX") {
        func = SINX;
    } else if (func_name == "COS2XINV") {
        func = COS2XINV;
    } else if (func_name == "X4M5") {
        func = X4M5;
    } else {
        std::cerr << "Error: Unknown function name '" << func_name << "'\n";
        return 1;
    }
    /////////////////

    /// MCI start ///
    double start = omp_get_wtime();
    master_function(42, func, xmin, xmax, samples);
    double end = omp_get_wtime();

    double duration = end - start;
    std::cout << "RUNTIME: " << duration << " s" << std::endl << std::endl;

    return 0;
}


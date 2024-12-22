#include <iostream>
#include <vector>
#include <random>
#include <omp.h>
#include <limits>
#include <chrono>

int calculate_minmax(const std::vector<long long>& vec) {

    std::cout << "Input threads number for omp: ";

    int min_val = std::numeric_limits<int>::max();
    int max_val = std::numeric_limits<int>::min();

    auto start = std::chrono::high_resolution_clock::now();
/*
#pragma omp parallel for reduction(min:min_val) reduction(max:max_val)
    for (int i = 0; i < vec.size(); ++i) {
        if (vec[i] < min_val) { min_val = vec[i]; }
        if (vec[i] > max_val) { max_val = vec[i]; }
    }
    */
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;


    std::cout << "Min: " << min_val << ", Max: " << max_val << std::endl;
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;
    return 0;
}

double simple_min_max(const std::vector<int>& vec, int& thread_num) {
    omp_set_num_threads(thread_num);

    auto start = std::chrono::high_resolution_clock::now();

    int min_val = std::numeric_limits<int>::max();
    int max_val = std::numeric_limits<int>::min();

    #pragma omp parallel
    {
        int local_min = std::numeric_limits<int>::max();
        int local_max = std::numeric_limits<int>::min();

        #pragma omp for
        for (int i = 0; i < vec.size(); ++i) {
            if (vec[i] < local_min) {
                local_min = vec[i];
            }
            if (vec[i] > local_max) {
                local_max = vec[i];
            }
        }

        #pragma omp critical
        {
            if (local_min < min_val) {
                min_val = local_min;
            }
            if (local_max > max_val) {
                max_val = local_max;
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    return duration.count();
}

double simple_min_max_2(const std::vector<long long>& vec, int& thread_num) {
    omp_set_num_threads(thread_num);

    auto start = std::chrono::high_resolution_clock::now();

    long long min_val = std::numeric_limits<long long>::max();
    long long max_val = std::numeric_limits<long long>::min();

#pragma omp parallel
    {
        long long local_min = std::numeric_limits<long long>::max();
        long long local_max = std::numeric_limits<long long>::min();

#pragma omp for
        for (int i = 0; i < vec.size(); ++i) {
            if (vec[i] < local_min) {
                local_min = vec[i];
            }
            if (vec[i] > local_max) {
                local_max = vec[i];
            }
        }

#pragma omp critical
        {
            if (local_min < min_val) {
                min_val = local_min;
            }
            if (local_max > max_val) {
                max_val = local_max;
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    return duration.count();
}

int main() {
    std::vector<int> thread_experiments = { 1, 2, 4, 8, 16, 24, 32, 64, 128, 256, 512 };
    std::vector<long long> vector_size_experiments = { 1'000'000,  1'000'000'000, 2'000'000'000 };
    int runs_count = 5;

    std::random_device rd; // Источник случайных чисел
    std::mt19937 gen(rd()); // Генератор случайных чисел Mersenne Twister
    std::uniform_int_distribution<long long> dist(INT_MIN, INT_MAX);

    std::cout << "Starting omp..." << std::endl;

    for (int ii = 0; ii < vector_size_experiments.size(); ii++)
    {
        long long vector_size_experiment = vector_size_experiments[ii];
        std::vector<long long> vec(vector_size_experiment);

        for (long long iii = 0; iii < vector_size_experiment; ++iii) {
            vec[iii] = dist(gen);
        }

        for (int i = 0; i < thread_experiments.size(); i++)
        {
            int current_thread_experiment = thread_experiments[i];

            double total_execution_time = 0;
            for (int j = 0; j < runs_count; j++)
            {
                total_execution_time += simple_min_max_2(vec, current_thread_experiment);
            }
            double avg_exexution_time = total_execution_time / runs_count;
            std::cout << "Vector size:" << vector_size_experiment << "; Threads: " << current_thread_experiment << "; Execution time " << avg_exexution_time << ";" << std::endl;
        }
    }

    std::cout << "Waiting for exit...";
    int temp;
    std::cin >> temp;
    return 0;
}
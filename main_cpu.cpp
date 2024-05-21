#include <iostream>
#include <boost/program_options.hpp>
#include <cmath>
#include <memory>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <chrono>

namespace opt = boost::program_options;


double linearInterpolation(double x, double x1, double y1, double x2, double y2) {
    return y1 + ((x - x1) * (y2 - y1) / (x2 - x1));
}


void initMatrix(std::unique_ptr<double[]> &arr ,int N){
    arr[0] = 10.0;
    arr[N-1] = 20.0;
    arr[(N-1)*N + (N-1)] = 30.0;
    arr[(N-1)*N] = 20.0;
    for (size_t i = 1; i < N-1; i++){
        arr[0*N+i] = linearInterpolation(i,0.0,arr[0],N-1,arr[N-1]);
        arr[i*N+0] = linearInterpolation(i,0.0,arr[0],N-1,arr[(N-1)*N]);
        arr[i*N+(N-1)] = linearInterpolation(i,0.0,arr[N-1],N-1,arr[(N-1)*N + (N-1)]);
        arr[(N-1)*N+i] = linearInterpolation(i,0.0,arr[(N-1)*N],N-1,arr[(N-1)*N + (N-1)]);
    }
}


int main(int argc, char const *argv[]){

    opt::options_description desc("опции");
    desc.add_options()
        ("accuracy",opt::value<double>()->default_value(1e-6),"точность")
        ("cellsCount",opt::value<int>()->default_value(256),"размер матрицы")
        ("iterCount",opt::value<int>()->default_value(1000000),"количество операций")
        ("help","помощь")
    ;
    opt::variables_map vm;
    opt::store(opt::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }
    opt::notify(vm);

    int N = vm["cellsCount"].as<int>();
    double accuracy = vm["accuracy"].as<double>();
    int countIter = vm["iterCount"].as<int>();
    double error = 1.0;
    int iter = 0;

    std::unique_ptr<double[]> A(new double[N*N]);
    std::unique_ptr<double[]> Anew(new double[N*N]);
    initMatrix(A,N);
    initMatrix(Anew,N);
    double* prevmatrix = Anew.get();
    double* curmatrix = A.get();


    auto start = std::chrono::high_resolution_clock::now();
    while (iter < countIter && iter<10000000 && error > accuracy){
        #pragma acc parallel loop independent collapse(2) vector vector_length(16) gang num_gangs(40)
        for (size_t i = 1; i < N-1; i++){
            for (size_t j = 1; j < N-1; j++){
                curmatrix[i*N+j]  = (prevmatrix[i*N+j+1] + prevmatrix[i*N+j-1] + prevmatrix[(i-1)*N+j] + prevmatrix[(i+1)*N+j])/4;
            }
        }
        if ((iter+1) %10 == 0){
            error = 0.0;
            #pragma acc parallel loop independent collapse(2) vector vector_length(16) gang num_gangs(40) reduction(max:error)
            for (size_t i = 1; i < N-1; i++){   
                for (size_t j = 1; j < N-1; j++){
                    error = fmax(error,fabs(curmatrix[i*N+j]-prevmatrix[i*N+j]));
                }
            }
        }

        double* temp = prevmatrix;
        prevmatrix = curmatrix;
        curmatrix = temp;   
        iter++;
    }
    auto end = std::chrono::high_resolution_clock::now();
    double time_s = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count())/1000;
    std::cout<<"time: " << time_s<<" error: "<<error << " iterarion: " << iter<<std::endl;
    if (N == 13 || N == 10){
        for (size_t i = 0; i < N; i++){
            for (size_t j = 0; j < N; j++){
                std::cout << A[i*N+j] << ' ';
                }
            std::cout << std::endl;
        }
    }
    A = nullptr;
    Anew = nullptr;
    return 0;
}
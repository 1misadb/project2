#include "geometry.h"
#include "json.hpp"
#include <fstream>
#include <random>
#include <iostream>
#include <cmath>

using json = nlohmann::json;

static Paths64 load_shape(const std::string& file){
    std::ifstream f(file);
    json j; f >> j;
    Paths64 shape;
    if(j.empty()) return shape;
    const auto& paths = j[0];
    for(const auto& path : paths){
        Path64 p;
        for(const auto& pt : path){
            double x_d = pt[0].get<double>();
            double y_d = pt[1].get<double>();
            long long x = static_cast<long long>(std::llround(x_d*1000.0));
            long long y = static_cast<long long>(std::llround(y_d*1000.0));
            p.push_back({x,y});
        }
        shape.push_back(std::move(p));
    }
    return shape;
}

int main(){
    std::vector<std::string> files = {"../part5.json","../part6.json","../part7.json","../part8.json"};
    std::vector<Paths64> shapes;
    for(const auto& f:files){
        shapes.push_back(load_shape(f));
    }
    if(shapes.size() < 2){
        std::cerr << "need at least two shapes" << std::endl;
        return 1;
    }
    std::mt19937 rng(123);
    std::uniform_int_distribution<int> dist(0, static_cast<int>(shapes.size()-1));
    for(int i=0;i<10;++i){
        int a = dist(rng);
        int b = dist(rng);
        while(b==a) b = dist(rng);
        std::vector<Paths64> others{shapes[b]};
        bool cpu = overlap(shapes[a], others[0]);
        auto gpu = overlapBatchGPU(shapes[a], others);
        if(gpu.size() != 1 || gpu[0] != cpu){
            std::cerr << "mismatch on pair " << a << "," << b << std::endl;
            return 1;
        }
    }
    std::cout << "OK" << std::endl;
    return 0;
}

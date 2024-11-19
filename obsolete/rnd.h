#pragma once

#include <random>

std::mt19937_64 mt;

std::normal_distribution<double> dist_gaussian;
double gaussian(){ return dist_gaussian(mt); }

std::uniform_real_distribution<> dist_01;
double uniform(){ return dist_01(mt); }


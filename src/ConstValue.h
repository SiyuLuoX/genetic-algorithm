#pragma once

#include <cstdio>
#include <climits>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <cmath>
//圆周率
#define M_PI 3.14159265358979323846 
//自然常数e
#define M_E  2.71828182845904523536

//x的维度
#define DIM 8
//ACKLEY函数的
#define A 20.0
#define B 0.20
#define C 2.0*M_PI

//基因互换采用的方法
#define CROSSOVER_TYPE ("Signal Point")
//发生交叉的几率，此处为100%
#define CROSSOVER_RATE 100
#define MUTATION_SIZE 0.01
//发生变异的几率，此处为10%
#define MUTATION_RATE 10
//父母基因的选择方法
#define PARENTS_SELECTION ("Roullete method")
//种群的个数
#define POPULATION_SIZE 200
//初始化方式随机
#define INITIALIZATION_TYPE ("Aleatory")
//迭代次数
#define CALLS_LIMIT 50000
#define STOP_CONDITION (" fitness calls OR optimum achived")

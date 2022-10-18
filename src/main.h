#include "Chromosomes.h"

//定义函数指针ParentSelectionFunction
typedef Chrom* (*ParentSelectionFunction)(Chrom*, Chrom*);
//定义函数指针CrossoverFunction
typedef void (*CrossoverFunction)(Chrom*, Chrom*);
//定义函数指针MutationFunction
typedef void (*MutationFunction)(Chrom*);



Chrom* generateInitialPopulation(Chrom* chroms);

void printInformations(Chrom* chroms);

double ackleysFunction(double* x, double c1, double c2, double c3);

Chrom* createGeneration(Chrom* chroms, ParentSelectionFunction f, CrossoverFunction c, MutationFunction m);

Chrom* roulleteSelection(Chrom* parent, Chrom* chroms);

void averageBetweenParentsCrossover(Chrom* parent, Chrom* chroms);

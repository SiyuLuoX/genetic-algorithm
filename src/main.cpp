#include "Chromosomes.h"
#include "Roullete.h"
#include "main.h"
#include <iostream>
#include <fstream>
using namespace std;
//0-1的随机浮点数
#define RAND_DOUBLE(interval) ((double) rand() * interval / ((double) RAND_MAX))


//生成随机的染色体
Chrom* generateInitialPopulation(Chrom* chroms)
{
	int i, j;
	chroms = new Chrom[POPULATION_SIZE];
	/* 初始化伪随机数生成器 */
	srand((unsigned int)time(NULL));

	for (i = 0; i < POPULATION_SIZE; i++) {
		for (j = 0; j < DIM; j++) {
			//0-5之间的随机数
			chroms[i].genes[j] = RAND_DOUBLE(5);
			//如果是偶数则变为负数
			if (rand() % 2 == 0)
				chroms[i].genes[j] *= (-1.0);
		}
	}

	return chroms;
}

void printInformations(Chrom* chroms) {
	printf("\nCrossover: %s.\n", CROSSOVER_TYPE);
	printf("Crossover rate: %d%s.\n", CROSSOVER_RATE, "%");
	printf("Mutation: +/- %lf.\n", MUTATION_SIZE);
	printf("Mutation rate: %d%s.\n", MUTATION_RATE, "%");
	printf("Parents selection: %s.\n", PARENTS_SELECTION);
	printf("Population size: %d.\n", POPULATION_SIZE);
	printf("Initialization: %s.\n", INITIALIZATION_TYPE);
	printf("Stop condition: %d %s.\n", CALLS_LIMIT, STOP_CONDITION);
	printf("Dimensions: %d.\n", DIM);
	printf("\n");
}

//计算fitness也就是ACKLEY函数的值
double ackleysFunction(double* x, double a, double b, double c)
{
	double sum1 = 0;
	double sum2 = 0;

	for (int i = 0; i < DIM; i++) {
		sum1 += x[i] * x[i];
		sum2 += cos(c * x[i]);
	}
	return (-1.0) * a * exp((-1.0) * b * sqrt((1.0 / DIM) * sum1)) - exp((1.0 / DIM) * sum2) + M_E + a;
}


Chrom* createGeneration(Chrom* chroms, ParentSelectionFunction f, CrossoverFunction c, MutationFunction m)
{
	Chrom* parent = new Chrom[POPULATION_SIZE];
	//selecting parent
	parent = f(parent, chroms);
	//crossover
	c(parent, chroms);
	//mutating offspring
	m(chroms);
	delete []parent;
	return chroms;
}


//轮盘赌(别人写的)
Chrom* roulleteSelection(Chrom* parent, Chrom* chroms)
{
	int i, j;
	double last = 0.0;
	double r;
	double sum = 0;
	double bigger = 0;
	Roullete* roullete = NULL;
	//求
	for (i = 0; i < POPULATION_SIZE; i++) {
		sum += fabs(chroms[i].fitness);
		bigger = (fabs(chroms[i].fitness) > bigger) ? fabs(chroms[i].fitness) : bigger;
	}
	roullete = new Roullete[POPULATION_SIZE + 1];
	roullete[0].value = 0.0;
	//适应度越好(低，轮盘赌空间越大
	for (i = 0; i < POPULATION_SIZE; i++) {
		roullete[i + 1].value = last + bigger / sum - (fabs(chroms[i].fitness) / sum);
		roullete[i].chrom = chroms[i];
		last += bigger / sum - (fabs(chroms[i].fitness) / sum);
	}

	for (i = 0; i < POPULATION_SIZE; i++) {
		r = RAND_DOUBLE(last);
		for (j = 0; j < POPULATION_SIZE; j++) {
			if (r >= roullete[j].value && r < roullete[j + 1].value) {
				memcpy(&(parent[i]), &(roullete[j].chrom), sizeof(Chrom));
				break;
			}
		}
	}
	delete[] roullete;
	return parent;
}


//轮盘赌
Chrom* roulleteSelection02(Chrom* parent, Chrom* chroms)
{
	int i, j;
	double last = 0.0;
	double r;
	double sum = 0;
	double bigger = 0;
	Roullete* roullete = NULL;
	//求
	for (i = 0; i < POPULATION_SIZE; i++) {
		sum += 1 / chroms[i].fitness;
	}
	roullete = new Roullete[POPULATION_SIZE + 1];
	roullete[0].value = 0.0;
	//适应度越好(低，轮盘赌空间越大
	for (i = 0; i < POPULATION_SIZE; i++) {
		roullete[i + 1].value = last + (1 / (chroms[i].fitness) / sum);
		roullete[i].chrom = chroms[i];
		last += (1 / chroms[i].fitness) / sum;
	}

	for (i = 0; i < POPULATION_SIZE; i++) {
		r = RAND_DOUBLE(last);
		for (j = 0; j < POPULATION_SIZE; j++) {
			if (r >= roullete[j].value && r < roullete[j + 1].value) {
				memcpy(&(parent[i]), &(roullete[j].chrom), sizeof(Chrom));
				break;
			}
		}
	}
	delete[] roullete;
	return parent;
}


//轮盘赌
Chrom* roulleteSelection03(Chrom* parent, Chrom* chroms)
{
	int i, j;
	double last = 0.0;
	double r;
	double sum = 0;
	double bigger = 0;
	Roullete* roullete = NULL;
	double bestFitness=INT_MAX;
	//求
	for (i = 0; i < POPULATION_SIZE; i++) {
		sum += 1/chroms[i].fitness;
		if (chroms[i].fitness<bestFitness){
			bestFitness = chroms[i].fitness;
		}
	}
	roullete = new Roullete[POPULATION_SIZE + 1];
	roullete[0].value = 0.0;
	//适应度越好(低，轮盘赌空间越大
	for (i = 0; i < POPULATION_SIZE; i++) {
		roullete[i + 1].value = last + (1/chroms[i].fitness) / sum;
		roullete[i].chrom = chroms[i];
		last += (1/chroms[i].fitness) / sum;
	}

	//复制其中最好的个体
	for (j = 0; j < POPULATION_SIZE; j++) {
		if (roullete[j].chrom.fitness=bestFitness) {
			memcpy(&(parent[0]), &(roullete[j].chrom), sizeof(Chrom));
			break;
		}
	}
	
	for (i = 1; i < POPULATION_SIZE; i++) {
		r = RAND_DOUBLE(last);
		for (j = 0; j < POPULATION_SIZE; j++) {
			if (r >= roullete[j].value && r < roullete[j + 1].value) {
				memcpy(&(parent[i]), &(roullete[j].chrom), sizeof(Chrom));
				break;
			}
		}
	}
	delete[] roullete;
	return parent;
}


//取双亲中的连号基因再求平均值
void averageBetweenParentsCrossover(Chrom* parent, Chrom* chroms) {
	int i, k;

	for (i = 0; i < POPULATION_SIZE; i++) {
		if (RAND_DOUBLE(100) <= CROSSOVER_RATE) {
			//assign new value to the offspring
			for (k = 0; k < DIM; k++)
				chroms[i].genes[k] = (parent[i].genes[k] + parent[(i + 1) % POPULATION_SIZE].genes[k]) / 2.0;
		}
		//crossover didnt happened, both chromossomes pass to the next generation without modifications
		else {
			for (k = 0; k < DIM; k++) {
				chroms[i].genes[k] = parent[i].genes[k];
			}
		}
	}
}


//单点杂交
void singlePointCrossover(Chrom* parent, Chrom* chroms) {
	int i, k;
	for (i = 0; i < POPULATION_SIZE; i++) {
		//闭区间[1,6]
		srand(unsigned(time(NULL)));
		int a = rand() % (DIM - 2) + 1;
		//百分百发生基因交叉
		if (RAND_DOUBLE(100) <= CROSSOVER_RATE) {
			for (k = 0; k < DIM; k++) {
				if (k <= a) {
					chroms[i].genes[k] = parent[i].genes[k];
				}else {
					chroms[i].genes[k] = parent[(i + 1) % POPULATION_SIZE].genes[k];
				}
			}
		}
		//不发生变异则直接复制
		else {
			for (k = 0; k < DIM; k++) {
				chroms[i].genes[k] = parent[i].genes[k];
			}
		}
	}
}


//基因发生变异
void sumMutation(Chrom* chroms) {
	int i, k;

	for (i = 0; i < POPULATION_SIZE; i++) {
		for (k = 0; k < DIM; k++) {
			if (RAND_DOUBLE(100) <= MUTATION_RATE)
				chroms[i].genes[k] += rand() % 2 == 0 ? MUTATION_SIZE : (-1.0) * MUTATION_SIZE;
		}
	}
}


//变异函数02
void sumMutation02(Chrom* chroms) {
	int i, k;

	for (i = 0; i < POPULATION_SIZE; i++) {
		for (k = 0; k < DIM; k++) {
			if (RAND_DOUBLE(100) <= MUTATION_RATE)
				chroms[i].genes[k] += rand() % 2 == 0 ? RAND_DOUBLE(1): (-1.0) * RAND_DOUBLE(1);
		}
	}
}

int main(int argc, char* argv[]) {
	int counter = 0; 
	double best = INT_MAX;
	double globalBest = INT_MAX;
	int i, j;
	Chrom globalChromBest;
	Chrom chromBest;
	Chrom* chroms = NULL;
	// 以写模式打开文件
	ofstream outfile;
	outfile.open("data.txt", ios::out | ios::trunc);

	chroms = generateInitialPopulation(chroms);

	printInformations(chroms);

	while (++counter <= CALLS_LIMIT) {
		//		printf("Generation %d\n", counter);
		//		printGeneration(chroms);

		for (i = 0; i < POPULATION_SIZE; i++) {
			chroms[i].fitness = ackleysFunction(chroms[i].genes, A, B, C);
			if (fabs(chroms[i].fitness) < fabs(best)) {
				best = chroms[i].fitness;
				chromBest = chroms[i];
			}

			//全局最小值,函数下界为0
			if (fabs(chroms[i].fitness) == 0) {
				printf("Chromossome genes: ");
				for (j = 0; j < DIM; j++)
					printf("%lf ", chroms[i].genes[j]);
				printf("- Iterations: %d - Ackely's Function = %.8lf\n", counter, chroms[i].fitness);
				delete []chroms;
				return 0;
			}
		}

		if (fabs(best)<fabs(globalBest)){
			globalBest = best;
			globalChromBest = chromBest;
		}

		//每一千轮输出最好结果
		if (counter % 1000 == 0) {
			printf("Best result at generation %d: %.8lf\n", counter, best);
			// 向文件写入用户输入的数据
			outfile << counter<<","<<best << endl;
			best = INT_MAX;
		}

		//开始新的一代
		//chroms = createGeneration(chroms, roulleteSelection, averageBetweenParentsCrossover, sumMutation);
		chroms = createGeneration(chroms, roulleteSelection02, averageBetweenParentsCrossover, sumMutation02);
	}

	// 关闭打开的文件
	outfile.close();

	printf("Best overall result: %.8lf\n**Genes**\n", globalBest);
	for (i = 0; i < DIM; i++)
		printf("%lf\n", globalChromBest.genes[i]);
	delete chroms;
	return 0;
}
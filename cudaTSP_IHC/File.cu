#include <cuda_runtime.h>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <thread>
#include <execution>
#include <concurrencysal.h>
#include <ppl.h>
#include <string>
#include <device_launch_parameters.h>

using namespace concurrency;
using namespace std;


const float boxsize = 10.f;
const int npoints = 1024; //GPU parallel broken past 1025 points.
const int to_possibilites = (npoints - 2) * (npoints - 1) / 2;
const int distarr_n = npoints * (npoints - 1) / 2;
FILE *fp = NULL;
FILE *gnupipe = NULL;

int variables[to_possibilites][2];
int variablesi[to_possibilites];
int variablesk[to_possibilites];
float results[to_possibilites];

struct point {
	float x, y;
	int number;
};
point List[npoints];
float Listx[npoints];
float Listy[npoints];
float distmatrix[npoints][npoints];
float vnnmatrix[npoints][npoints];

float distarray[npoints * npoints];
float bestdist = FLT_MAX;

int Path[npoints + 1];

/*
void init() {
	distmatrix = (float**)malloc(npoints * sizeof(float*));
	vnnmatrix = (float**)malloc(npoints * sizeof(float*));
	for (int j = 0; j < npoints; j++) {
		distmatrix[j] = (float*)malloc(npoints * sizeof(float));
		vnnmatrix[j] = (float*)malloc(npoints * sizeof(float));
	}
}
*/

void PrintMatrix() {
	for (int i = 0; i < npoints; i++) {
		for (int j = 0; j < npoints; j++) {
			cout << distmatrix[i][j] << "\t \t";
		}
		cout << "\n";
	}
	cout << "\n";
}

void RandomCoordinate(point * a, int number) {
	float x = boxsize * (float)rand() / (RAND_MAX);
	float y = boxsize * (float)rand() / (RAND_MAX);
	Listx[number] = x;
	Listy[number] = y;
	a->x = x;
	a->y = y;
	a->number = number;
}

float CalcDist(point point1, point point2) {
	return sqrt(pow(point2.x - point1.x, 2) + pow(point2.y - point1.y, 2));
}

float GetDist(int x, int y) {
	if (x < y) return distarray[npoints*x - (x + 1) * (x + 2) / 2 + y];
	if (x > y) return distarray[npoints*y - (y + 1) * (y + 2) / 2 + x];
	else return FLT_MAX;
}

float CalcDistPath2(int* Path2) {
	float total = 0;
	for (int i = 0; i < npoints; i++) {
		total += distmatrix[Path[i]][Path[i + 1]];
	}
	return total;
}

float CalcDistPath3(int* Path2) {
	float total = 0;
	for (int i = 0; i < npoints; i++) {
		total += distmatrix[Path[i]][Path[i + 1]];
	}
	return total;
}

void FillVariables() {
	int j = 0;
	for (int i = 1; i < npoints - 1; i++) {
		for (int k = i + 1; k < npoints; k++) {
			variables[j][0] = i;
			variables[j][1] = k;
			variablesi[j] = i;
			variablesk[j] = k;
			j++;
		}
	}
}

int* GetNumPath(point* Path) {
	int NumPath[npoints + 1];
	for (int i = 0; i < npoints + 1; i++) {
		NumPath[i] = Path[i].number;
	}
	return NumPath;
}

void PopulateList() {
	for (int i = 0; i < npoints; i++) {
		RandomCoordinate(&List[i],i);
	}
}

void MakeDistMatrix() {
	int f = 0;
	for (int i = 0; i < npoints; i++) {
		for (int j = 0; j < npoints; j++) {
			if (i == j) {
				distmatrix[i][j] = FLT_MAX;
				distarray[f] = FLT_MAX;
				f++;
			}	
			else {
				float dist = CalcDist(List[i], List[j]);
				distmatrix[i][j] = dist;
				distarray[f] = dist;
				f++;	
			}
		}
	}
}

void SetCol(int col) {
	for (int i = 0; i < npoints; i++) {
		vnnmatrix[i][col] = FLT_MAX;
	}
}

void two_opt(int i, int k) {
	while (i < k) {
		int temp = Path[i];
		Path[i] = Path[k];
		Path[k] = temp;
		i++;
		k--;
	}
}

void Graph() {
	char* GnuCommands[] = { "set title \"Start\"", "plot 'data.tmp' pointtype 7 with linespoints" };
	fp = fopen("data.tmp", "w");
	gnupipe = _popen("gnuplot -persistent", "w");

	for (int i = 0; i < npoints; i++) {
		fprintf(fp, "%f %f\n", List[i].x, List[i].y);
	}

	for (int i = 0; i < 2; i++) {
		fprintf(gnupipe, "%s\n", GnuCommands[i]);
	}
}

void VNNGraph() {
	string GnuCommands2[] = { "set title \"VNN " + to_string(bestdist) + "\"", "plot 'data2.tmp' pointtype 7 with linespoints"};
	fp = fopen("data2.tmp", "w");
	gnupipe = _popen("gnuplot -persistent", "w");

	for (int i = 0; i < npoints + 1; i++) {
		fprintf(fp, "%f %f\n", List[Path[i]].x, List[Path[i]].y);
	}

	for (int i = 0; i < 2; i++) {
		fprintf(gnupipe, "%s\n", GnuCommands2[i].c_str());
	}
}

void IHCGraph() {
	string GnuCommands3[] = { "set title \"IHC " + to_string(bestdist) + "\"", "plot 'data3.tmp' pointtype 7 with linespoints" };
	fp = fopen("data3.tmp", "w");
	gnupipe = _popen("gnuplot -persistent", "w");

	for (int i = 0; i < npoints + 1; i++) {
		fprintf(fp, "%f %f\n", List[Path[i]].x, List[Path[i]].y);
	}

	for (int i = 0; i < 2; i++) {
		fprintf(gnupipe, "%s\n", GnuCommands3[i].c_str());
	}
}


void vnn(point points[]) {
	Path[0] = List[0].number;

	memcpy(vnnmatrix, distmatrix, npoints * npoints * sizeof(float));

	SetCol(0);
	for (int i = 0; i < npoints; i++) {
		float* min = min_element(vnnmatrix[Path[i]], vnnmatrix[Path[i]] + npoints);
		int index = min - &vnnmatrix[Path[i]][0];

		Path[i + 1] = List[index].number;
		SetCol(index);

	}
	Path[npoints] = List[0].number;

}


void IHC() {
	int j = 0;
	while (j < 1000) {
		int BestNewPath[npoints + 1];
		float BestNewDist = FLT_MAX;
		int bestvalue[2];

		for (int i = 1; i < npoints - 1; i++) {
			for (int k = i + 1; k < npoints; k++) {
				float dist = bestdist - distmatrix[Path[i - 1]][Path[i]] - distmatrix[Path[k]][Path[k + 1]];
				dist += distmatrix[Path[i - 1]][Path[k]] + distmatrix[Path[i]][Path[k + 1]];

				if (dist < bestdist) {
					bestvalue[0] = i;
					bestvalue[1] = k;

					BestNewDist = dist;
				}
			}
		}
		if (BestNewDist != FLT_MAX) {
			two_opt(bestvalue[0], bestvalue[1]);
			bestdist = BestNewDist;
			
		}
		else {
			cout << " done" << " " << j;
			break;

		}
		j++;
	}
}

void IHC_parallel() {
	int j = 0;
	while (j < 1000) {
		int BestNewPath[npoints + 1];
		float BestNewDist = FLT_MAX;
		int bestvalue[2];

		parallel_for_each(begin(variables), end(variables), [&](int value[]) {
			float dist = bestdist - distmatrix[Path[value[0] - 1]][Path[value[0]]] - distmatrix[Path[value[1]]][Path[value[1] + 1]];
			dist +=	distmatrix[Path[value[0] - 1]][Path[value[1]]] + distmatrix[Path[value[0]]][Path[value[1] + 1]];

			if (dist < bestdist) {
				bestvalue[0] = value[0];
				bestvalue[1] = value[1];

				BestNewDist = dist;
			}

		});

		if (BestNewDist != FLT_MAX) {
			two_opt(bestvalue[0], bestvalue[1]);
			bestdist = BestNewDist;
		}
		else {
			cout << " done" << " " << j;
			break;
		}
		j++;
	}
}

__global__ void GPU_parallel(int* Path, int N, float* results, float bestdist, float* distmatrix, int* variablesi, int* variablesk) {
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	int index = row * N + col;

	results[index] = bestdist - distmatrix[Path[variablesi[index] -1] * N + Path[variablesi[index]]] - distmatrix[Path[variablesk[index]] * N + Path[variablesk[index] + 1]];
	results[index] += distmatrix[Path[variablesi[index] - 1] * N + Path[variablesk[index]]] + distmatrix[Path[variablesi[index]] * N + Path[variablesk[index] + 1]];
}

void IHC_CUDA() {
	


	size_t bytes_matrix = npoints * npoints * sizeof(float);
	size_t bytes_list = npoints * sizeof(float);
	size_t bytes_variables = to_possibilites * sizeof(int);
	size_t bytes_path = (npoints + 1) * sizeof(int);
	size_t bytes_float = sizeof(float);

	float* d_matrix, * d_listx, * d_listy, *d_results;
	int* d_variablesi, * d_variablesk;
	float d_bestdist;
	cudaMalloc(&d_matrix, bytes_matrix);
	cudaMalloc(&d_listx, bytes_list);
	cudaMalloc(&d_listy, bytes_list);
	cudaMalloc(&d_variablesi, bytes_variables);
	cudaMalloc(&d_variablesk, bytes_variables);
	cudaMalloc(&d_results, bytes_variables);
	

	int THREADS = npoints - 1;

	int BLOCKS = to_possibilites / THREADS;

	dim3 threads(THREADS, THREADS);
	dim3 blocks(BLOCKS, BLOCKS);

	cudaMemcpy(d_matrix, distarray, bytes_matrix, cudaMemcpyHostToDevice);
	cudaMemcpy(d_listx, Listx, bytes_list, cudaMemcpyHostToDevice);
	cudaMemcpy(d_listy, Listy, bytes_list, cudaMemcpyHostToDevice);
	cudaMemcpy(d_variablesi, variablesi, bytes_variables, cudaMemcpyHostToDevice);
	cudaMemcpy(d_variablesk, variablesk, bytes_variables, cudaMemcpyHostToDevice);
	cudaMemcpy(d_results, results, bytes_variables, cudaMemcpyHostToDevice);

	int* d_Path;
	cudaMalloc(&d_Path, bytes_path);

	int j = 0;
	while (j < 500) {
		int BestNewPath[npoints + 1];
		float BestNewDist = FLT_MAX;
		int bestvaluei;
		int bestvaluek;

		
		//const clock_t begin2 = clock();
		
		cudaMemcpy(d_Path, Path, bytes_path, cudaMemcpyHostToDevice);
		
		GPU_parallel<<<BLOCKS, THREADS>>>(d_Path, npoints, d_results, bestdist, d_matrix, d_variablesi, d_variablesk);
		
		cudaMemcpy(results, d_results, bytes_variables, cudaMemcpyDeviceToHost);
		
		//cout << "\nTime taken for " << npoints << " in linear: " << float(clock() - begin2) / CLOCKS_PER_SEC;

		float* bestvalue = std::min_element(std::begin(results), std::end(results));
		int bestindex =  std::distance(std::begin(results), bestvalue);

		if (bestvalue[0] < bestdist) {
			two_opt(variablesi[bestindex], variablesk[bestindex]);
			bestdist = bestvalue[0];
		}
		else {
			cout << " done" << " " << j;
			break;
		}

		j++;
	}

	cudaFree(d_matrix);
	cudaFree(d_listx);
	cudaFree(d_listy);
	cudaFree(d_variablesi);
	cudaFree(d_variablesk);
	cudaFree(d_results);

	cudaFree(d_Path);
}





int main() {
	
	//init();
	PopulateList();
	
	
	MakeDistMatrix();


	
	//Graph();
	//bestdist = CalcDistPath(Path);
	vnn(List);
	//VNNGraph();
	FillVariables();
	
	// LINEAR
	bestdist = CalcDistPath3(Path);
	const clock_t begin = clock();
	IHC();
	cout << "\nTime taken for " << npoints << " in linear: " << float(clock() - begin) / CLOCKS_PER_SEC;
	IHCGraph();
	
	
	// PARALLEL CPU
	vnn(List);
	bestdist = CalcDistPath3(Path);
	const clock_t begin2 = clock();
	IHC_parallel();
	cout << "\nTime taken for " << npoints << " in CPU parallel: " << float(clock() - begin2) / CLOCKS_PER_SEC;
	IHCGraph();
	

	
	// PRALLEL GPU
	vnn(List);
	bestdist = CalcDistPath3(Path);

	const clock_t begin3 = clock();
	IHC_CUDA();
	cout << "\nTime taken for " << npoints << " in GPU parallel: " << float(clock() - begin3) / CLOCKS_PER_SEC;
	IHCGraph();
	
}
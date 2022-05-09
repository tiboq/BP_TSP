#include <cuda_runtime.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
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
#include <cublas_v2.h>
#ifdef __INTELLISENSE__
#include "intellisense_cuda_intrinsics.h"
#endif


using namespace concurrency;
using namespace std;


const float boxsize = 10.f;
const int npoints = 280; //GPU parallel broken past 1025 points.
const int to_possibilites = (npoints - 2) * (npoints - 1) / 2;
const int distarr_n = npoints * (npoints - 1) / 2;
const float approxdist = 0.765 * sqrt(npoints * boxsize * boxsize);
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
cublasHandle_t cublasHandle;
cublasStatus_t cublasStatus;

vector<float> distances;

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

void SetCoordinate(point * a, int number, float x, float y) {
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
	//for (int i = npoints - 2; i > 0; i--) {
		for (int k = i + 1; k < npoints; k++) {
		//for (int k = npoints -1; k > i; k--) {
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

void DistGraph() {

	string GnuCommands4[] = { "set title \"IHC " + to_string(bestdist) + "\"", "plot 'data4.tmp' pointtype 7 with linespoints\n replot " + to_string(approxdist) };
	fp = fopen("data4.tmp", "w");
	gnupipe = _popen("gnuplot -persistent", "w");

	for (int i = 0; i < size(distances); i++) {
		fprintf(fp, "%i %f\n", i, distances[i]);
	}
	//fprintf(fp, "%s\n", "plot 20");


	for (int i = 0; i < 2; i++) {
		fprintf(gnupipe, "%s\n", GnuCommands4[i].c_str());
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

int GetThreadcount() {
	if (to_possibilites < 1025) {
		return (npoints - 1);
	} 
	else {
		for (int i = 1024; i > 1; i--) {
			if ((to_possibilites % i) == 0) return i;
		}
	}
}



void IHC() {
	int j = 0;
	while (j < 10000) {
		int BestNewPath[npoints + 1];
		float BestNewDist = FLT_MAX;
		int bestvalue[2];

		//for (int i = 1; i < npoints - 1; i++) {
		for (int i = npoints - 2; i > 1; i--) {
			//for (int k = i + 1; k < npoints; k++) {
			for (int k = npoints -1; k > i + 1; k--) {
				float dist = bestdist - distmatrix[Path[i - 1]][Path[i]] - distmatrix[Path[k]][Path[k + 1]];
				dist += distmatrix[Path[i - 1]][Path[k]] + distmatrix[Path[i]][Path[k + 1]];

				if ((dist < bestdist) && (bestdist - dist) > 0.00002f) { //prevent rounding error
					bestvalue[0] = i;
					bestvalue[1] = k;

					BestNewDist = dist;

					///*
					two_opt(bestvalue[0], bestvalue[1]);
					bestdist = BestNewDist;
					distances.push_back(BestNewDist);

					i = 1;
					k = 2;
					//*/
					
				}
			}
		}
		if (BestNewDist != FLT_MAX) {
			/*
			two_opt(bestvalue[0], bestvalue[1]);
			bestdist = BestNewDist;
			distances.push_back(BestNewDist);
			*/
		}
		else {
			cout << "linear done" << " " << j;
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

			if ((dist < bestdist) && (bestdist - dist) > 0.00002f) {
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
			cout << "CPU done" << " " << j;
			break;
		}
		j++;
	}
}



__global__ void GPU_parallel(int* Path, int N, float* results, float bestdist, float* distmatrix, int* variablesi, int* variablesk) {
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	int index = row * N + col;

	float dist = bestdist - distmatrix[Path[variablesi[index] - 1] * N + Path[variablesi[index]]] - distmatrix[Path[variablesk[index]] * N + Path[variablesk[index] + 1]];
	dist += distmatrix[Path[variablesi[index] - 1] * N + Path[variablesk[index]]] + distmatrix[Path[variablesi[index]] * N + Path[variablesk[index] + 1]];
	results[index] = dist;

}

void IHC_CUDA() {
	const int THREADS = GetThreadcount();
	const int BLOCKS = to_possibilites / THREADS;

	size_t bytes_matrix = npoints * npoints * sizeof(float);
	size_t bytes_list = npoints * sizeof(float);
	size_t bytes_variables = to_possibilites * sizeof(int);
	size_t bytes_path = (npoints + 1) * sizeof(int);
	size_t bytes_results2 = BLOCKS * sizeof(float);

	float* d_matrix, * d_listx, * d_listy, *d_results;
	int* d_variablesi, * d_variablesk;
	cudaMalloc(&d_matrix, bytes_matrix);
	cudaMalloc(&d_listx, bytes_list);
	cudaMalloc(&d_listy, bytes_list);
	cudaMalloc(&d_variablesi, bytes_variables);
	cudaMalloc(&d_variablesk, bytes_variables);
	cudaMalloc(&d_results, bytes_variables);
	
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

	
	//const clock_t begin3 = clock();

	int j = 0;
	while (j < 1000) {

		
		cudaMemcpy(d_Path, Path, bytes_path, cudaMemcpyHostToDevice);

		GPU_parallel<<<BLOCKS, THREADS>>>(d_Path, npoints, d_results, bestdist, d_matrix, d_variablesi, d_variablesk);	

		int index = 0;
		cublasIsamin(cublasHandle, to_possibilites, d_results, 1, &index);
		index--;

		float dist = bestdist - distarray[Path[variablesi[index] - 1] * npoints + Path[variablesi[index]]] - distarray[Path[variablesk[index]] * npoints + Path[variablesk[index] + 1]];
		dist += distarray[Path[variablesi[index] - 1] * npoints + Path[variablesk[index]]] + distarray[Path[variablesi[index]] * npoints + Path[variablesk[index] + 1]];

		if ((dist < bestdist) && (bestdist - dist) > 0.00002f) {
			two_opt(variablesi[index], variablesk[index]);
			bestdist = dist;
			//distances.push_back(dist);
		}
		else {
			cout << "GPU done" << " " << j;
			break;
		}
		
		j++;
	}

	//cout << "\nTime taken for " << npoints << " in GPU parallel: " << float(clock() - begin3) / CLOCKS_PER_SEC;
	cudaFree(d_matrix);
	cudaFree(d_listx);
	cudaFree(d_listy);
	cudaFree(d_variablesi);
	cudaFree(d_variablesk);
	cudaFree(d_results);
	cudaFree(d_Path);
}

void readTSPLib() {
	string line;
	ifstream myfile("TSPlib/a280.tsp");
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::bind1st(std::not_equal_to<char>(), ' ')));
			string delimiter = " ";
			int index = std::stoi(line.substr(0, line.find(delimiter))) - 1;

			string line2 = line.substr(line.find(delimiter), line.size() - 1);
			line2.erase(line2.begin(), std::find_if(line2.begin(), line2.end(), std::bind1st(std::not_equal_to<char>(), ' ')));
			float x = std::stof(line2.substr(0, line2.find(delimiter)));

			string line3 = line2.substr(line2.find(delimiter), line2.size() - 1);
			line3.erase(line3.begin(), std::find_if(line3.begin(), line3.end(), std::bind1st(std::not_equal_to<char>(), ' ')));
			float y = std::stof(line3.substr(0, line3.find(delimiter)));

			SetCoordinate(&List[index], index, x, y);

		}
		myfile.close();
	}
	else cout << "Unable to open file";
}



int main() {
	
	

	readTSPLib();
	
	
	//PopulateList();
	
	MakeDistMatrix();

	
	for (int i = 0; i < (npoints + 1); i++) {
		Path[i] = i;
	}
	Path[npoints] = 0;
	

	Graph();
	//bestdist = CalcDistPath(Path);
	
	//VNNGraph();
	FillVariables();
	
	//for (int i = 0; i < 10; i++){
		//cout << "\n" << "ITERATION " << i << "\n";
		
		bestdist = CalcDistPath3(Path);
		cout << bestdist;
		// LINEAR
		//vnn(List);
		bestdist = CalcDistPath3(Path);
		//VNNGraph();
		
		const clock_t begin = clock();
		IHC();
		//cout << "\nTime taken for " << npoints << " in linear: " << float(clock() - begin) / CLOCKS_PER_SEC;
		cout << " " << float(clock() - begin) / CLOCKS_PER_SEC << "\t";
		IHCGraph();
		DistGraph();
		
		/*
		// PARALLEL CPU
		vnn(List);
		bestdist = CalcDistPath3(Path);
		const clock_t begin2 = clock();
		IHC_parallel();
		//cout << "\nTime taken for " << npoints << " in CPU parallel: " << float(clock() - begin2) / CLOCKS_PER_SEC;
		cout << " " << float(clock() - begin2) / CLOCKS_PER_SEC << "\t";
		//IHCGraph();
		

	
		// PRALLEL GPU
		int* d_test;
		cudaMalloc(&d_test, sizeof(int)); // init gpu
		cublasStatus = cublasCreate(&cublasHandle); // this take about .7s

		const clock_t begin3 = clock();
		MakeDistMatrix();
		vnn(List);
		bestdist = CalcDistPath3(Path);
		
		IHC_CUDA();
		//cout << "\nTime taken for " << npoints << " in GPU parallel: " << float(clock() - begin3) / CLOCKS_PER_SEC;
		cout << " " << float(clock() - begin3) / CLOCKS_PER_SEC << "\t";
		IHCGraph();
		//DistGraph();
		*/
	//}
}
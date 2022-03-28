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

using namespace concurrency;
using namespace std;

const float boxsize = 10.f;
const int npoints = 500;
const int to_possibilites = (npoints - 2) * (npoints - 1) / 2;
const int distarr_n = npoints * (npoints - 1) / 2;
FILE *fp = NULL;
FILE *gnupipe = NULL;

vector<int*> variables;
pair<int, int> variables2[to_possibilites];
int variables3[to_possibilites][2];

struct point {
	float x, y;
	int number;
};
point List[npoints];
float distmatrix[npoints][npoints];
float vnnmatrix[npoints][npoints];
//float** distmatrix;
//float** vnnmatrix;



float distarray[distarr_n];
point Path[npoints + 1];
float bestdist = FLT_MAX;

int Path2[npoints + 1];

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
	a->x = boxsize * (float)rand() / (RAND_MAX);
	a->y = boxsize * (float)rand() / (RAND_MAX);
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

float CalcDistPath(point* Path) {
	float total = 0;
	for (int i = 0; i < npoints; i++) {
		//total += CalcDist(Path[i], Path[i + 1]);
		
		total += distmatrix[Path[i].number][Path[i + 1].number];
	}
	return total;
}

float CalcDistPath2(int* Path2) {
	float total = 0;
	for (int i = 0; i < npoints; i++) {
		total += distmatrix[Path2[i]][Path2[i + 1]];
	}
	return total;
}

float CalcDistPath3(int* Path2) {
	float total = 0;
	for (int i = 0; i < npoints; i++) {
		total += distmatrix[Path2[i]][Path2[i + 1]];
	}
	return total;
}

float CalcDistPath4(int* Path2, float prev_dist, int i, int k) {
	float total = prev_dist - distmatrix[Path2[i]][Path2[i - 1]];
	return 0;
}

void FillVariables() {
	for (int i = 1; i < npoints - 1; i++) {
		for (int k = i + 1; k < npoints; k++) {
			variables.push_back(new int[2]{ i,k });
		}
	}
}

void FillVariables2() {
	int j = 0;
	for (int i = 1; i < npoints - 1; i++) {
		for (int k = i + 1; k < npoints; k++) {
			variables2[j].first = i;
			variables2[j].second = k;
			j++;
		}
	}
}

void FillVariables3() {
	int j = 0;
	for (int i = 1; i < npoints - 1; i++) {
		for (int k = i + 1; k < npoints; k++) {
			variables3[j][0] = i;
			variables3[j][1] = k;
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
			}	
			else {
				float dist = CalcDist(List[i], List[j]);
				distmatrix[i][j] = dist;
				if (i < j) {
					distarray[f] = dist;
					f++;
				}	
			}
		}
	}
}

void SetCol(int col) {
	for (int i = 0; i < npoints; i++) {
		vnnmatrix[i][col] = FLT_MAX;
	}
}

point* two_opt(int i, int k, point * list) {
	while (i < k) {
		point temp = list[i];
		list[i] = list[k];
		list[k] = temp;
		i++;
		k--;
	}
	return list;
}

int* two_opt2(int i, int k, int path[]) {
	while (i < k) {
		int temp = path[i];
		path[i] = path[k];
		path[k] = temp;
		i++;
		k--;
	}
	return path;
}

void two_opt3(int i, int k) {
	while (i < k) {
		int temp = Path2[i];
		Path2[i] = Path2[k];
		Path2[k] = temp;
		i++;
		k--;
	}
}

void Graph() {
	char* GnuCommands[] = { "set title \"Start\"", "plot 'data.tmp' pointtype 7 with linespoints" };
	fp = fopen("data.tmp", "w");
	gnupipe = _popen("gnuplot -persistent", "w");

	for (int i = 0; i < npoints; i++) {
		//cout << List[i].number << " " << List[i].x << " " << List[i].y << "\n";
	}

	for (int i = 0; i < npoints; i++) {
		//cout << pointlist[i].x << '\t' << pointlist[i].y << '\n';
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
		//cout << pointlist[i].x << '\t' << pointlist[i].y << '\n';
		fprintf(fp, "%f %f\n", Path[i].x, Path[i].y);
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
		//cout << pointlist[i].x << '\t' << pointlist[i].y << '\n';
		//fprintf(fp, "%f %f\n", List[Path[i].number].x, List[Path[i].number].y);
		fprintf(fp, "%f %f\n", List[Path2[i]].x, List[Path2[i]].y);
	}

	for (int i = 0; i < 2; i++) {
		fprintf(gnupipe, "%s\n", GnuCommands3[i].c_str());
	}
}

void UpdateIHCGraph() {
	fp = fopen("data3.tmp", "w");
	fprintf(gnupipe, "%s\n", "clear");
	for (int i = 0; i < npoints + 1; i++) {
		
		fprintf(fp, "%f %f\n", Path[i].x, Path[i].y);
	}
	fprintf(gnupipe, "%s\n", "replot");
}

void vnn(point points[]) {
	Path[0] = List[0];

	//float vnnmatrix[npoints][npoints];
	memcpy(vnnmatrix, distmatrix, npoints * npoints * sizeof(float));

	SetCol(0);

	for (int i = 0; i < npoints; i++) {
		float* min = min_element(vnnmatrix[Path[i].number], vnnmatrix[Path[i].number] + npoints);
		int index = min - &vnnmatrix[Path[i].number][0];
		//cout << "\n" << List[index].number;
		Path[i + 1] = List[index];
		SetCol(index);
		
	}
	Path[npoints] = List[0];

}

void vnn2(point points[]) {

	Path2[0] = List[0].number;

	//float vnnmatrix[npoints][npoints];
	memcpy(vnnmatrix, distmatrix, npoints * npoints * sizeof(float));

	SetCol(0);

	for (int i = 0; i < npoints; i++) {
		float* min = min_element(vnnmatrix[Path2[i]], vnnmatrix[Path2[i]] + npoints);
		int index = min - &vnnmatrix[Path2[i]][0];
		//cout << "\n" << List[index].number;
		Path2[i + 1] = List[index].number;
		SetCol(index);

	}
	Path2[npoints] = List[0].number;

}

void IHC() {
	int j = 0;
	while (j < 1000) {
		point BestNewPath[npoints + 1];
		float BestNewDist = FLT_MAX;

		for (int i = 1; i < npoints - 1; i++) {
			for (int k = i + 1; k < npoints; k++) {
				point TempPath[npoints + 1];
				memcpy(TempPath, Path, (npoints + 1) * sizeof(point));
				two_opt(i, k, TempPath);
				float dist = CalcDistPath(TempPath);
				
				if (dist < bestdist) {
					BestNewDist = dist;
					memcpy(BestNewPath, TempPath, (npoints + 1) * sizeof(point));
				}
			}
		}
		if (BestNewDist != FLT_MAX) {
			memcpy(Path, BestNewPath, (npoints + 1) * sizeof(point));
			bestdist = BestNewDist;
			/*
			cout << "\n" << "new best path:" << "\n";
			for (int j = 0; j < npoints; j++) {
				cout << Path[j].number << " ";
			}
			cout << "\n" << bestdist;
			*/
			//UpdateIHCGraph();
		}
		else {
			cout << " done" << " " << j;
			break;

		}
		j++;
	}
}

void IHC2() {
	int j = 0;
	while (j < 1000) {
		int BestNewPath[npoints + 1];
		float BestNewDist = FLT_MAX;
		int bestvalue[2];

		for (int i = 1; i < npoints - 1; i++) {
			for (int k = i + 1; k < npoints; k++) {
				float dist = bestdist - distmatrix[Path2[i - 1]][Path2[i]] - distmatrix[Path2[k]][Path2[k + 1]];
				dist += distmatrix[Path2[i - 1]][Path2[k]] + distmatrix[Path2[i]][Path2[k + 1]];

				if (dist < bestdist) {
					bestvalue[0] = i;
					bestvalue[1] = k;

					//memcpy(bestvalue, value, 2 * sizeof(int));
					BestNewDist = dist;
				}
			}
		}
		if (BestNewDist != FLT_MAX) {
			two_opt3(bestvalue[0], bestvalue[1]);
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
		point BestNewPath[npoints + 1];
		float BestNewDist = FLT_MAX;

		

		parallel_for_each(begin(variables3), end(variables3), [&](int value[]) {
			point TempPath[npoints + 1];
			memcpy(TempPath, Path, (npoints + 1) * sizeof(point));
			
			two_opt(value[0], value[1], TempPath);
			float dist = CalcDistPath(TempPath);

			if (dist < bestdist) {
				BestNewDist = dist;
				memcpy(BestNewPath, TempPath, (npoints + 1) * sizeof(point));
			}

		});

		if (BestNewDist != FLT_MAX) {
			memcpy(Path, BestNewPath, (npoints + 1) * sizeof(point));
			bestdist = BestNewDist;
		}
		else {
			cout << " done" << " " << j;
			break;
		}
		j++;
	}
}

void IHC_parallel2() {
	int j = 0;
	while (j < 1000) {
		int BestNewPath[npoints + 1];
		float BestNewDist = FLT_MAX;

		
		parallel_for_each(begin(variables3), end(variables3), [&](int value[]) {
			int TempPath[npoints + 1];
			memcpy(TempPath, Path2, (npoints + 1) * sizeof(int));

			float dist = bestdist - distmatrix[TempPath[value[0]]][TempPath[value[0] - 1]] - distmatrix[TempPath[value[1]]][TempPath[value[1] + 1]];
			two_opt2(value[0], value[1], TempPath);
			dist += distmatrix[TempPath[value[0]]][TempPath[value[0] - 1]] + distmatrix[TempPath[value[1]]][TempPath[value[1] + 1]];
			
			//float dist = CalcDistPath3(TempPath);

			if (dist < bestdist) {
				BestNewDist = dist;
				memcpy(BestNewPath, TempPath, (npoints + 1) * sizeof(int));
			}

			});

		if (BestNewDist != FLT_MAX) {
			memcpy(Path2, BestNewPath, (npoints + 1) * sizeof(int));
			bestdist = BestNewDist;
		}
		else {
			cout << " done" << " " << j;
			break;
		}
		j++;
	}
}

void IHC_parallel3() {
	int j = 0;
	while (j < 1000) {
		int BestNewPath[npoints + 1];
		float BestNewDist = FLT_MAX;
		int bestvalue[2];

		parallel_for_each(begin(variables3), end(variables3), [&](int value[]) {
			float dist = bestdist - distmatrix[Path2[value[0] - 1]][Path2[value[0]]] - distmatrix[Path2[value[1]]][Path2[value[1] + 1]];
			dist +=	distmatrix[Path2[value[0] - 1]][Path2[value[1]]] + distmatrix[Path2[value[0]]][Path2[value[1] + 1]];

			if (dist < bestdist) {
				bestvalue[0] = value[0];
				bestvalue[1] = value[1];

				//memcpy(bestvalue, value, 2 * sizeof(int));
				BestNewDist = dist;
			}

			});

		if (BestNewDist != FLT_MAX) {
			two_opt3(bestvalue[0], bestvalue[1]);
			bestdist = BestNewDist;
		}
		else {
			cout << " done" << " " << j;
			break;
		}
		j++;
	}
}


int main() {
	//init();
	PopulateList();
	MakeDistMatrix();

	
	//Graph();
	//vnn(List);
	//bestdist = CalcDistPath(Path);
	//VNNGraph();
	
	
	//PrintMatrix();
	
	
	
	
	
	FillVariables();
	FillVariables2();
	FillVariables3();
	
	
	vnn2(List);
	bestdist = CalcDistPath3(Path2);
	const clock_t begin = clock();
	IHC2();
	cout << "\nTime taken for " << npoints << " in linear: " << float(clock() - begin) / CLOCKS_PER_SEC;
	IHCGraph();
	
	
	/*
	vnn(List);
	bestdist = CalcDistPath(Path);
	const clock_t begin2 = clock();
	IHC_parallel();
	cout << "\nTime taken for " << npoints << " in parallel1: " << float(clock() - begin2) / CLOCKS_PER_SEC;
	//IHCGraph();
	*/
	
	
	vnn2(List);
	bestdist = CalcDistPath3(Path2);
	const clock_t begin3 = clock();
	IHC_parallel3();
	cout << "\nTime taken for " << npoints << " in parallel2: " << float(clock() - begin3) / CLOCKS_PER_SEC;
	IHCGraph();
	
}
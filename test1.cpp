#include "test1.h"

void init(int **perProcess, int *startLine, double **matrix, double **b, double **x, int size, int rank) {
    *perProcess = new int[size]();
    for(int i = 0, tmp = size - (N % size); i < size; ++i) {
        (*perProcess)[i] = i < tmp ? (N / size) : (N / size + 1);
        if(i < rank) {
            *startLine += (*perProcess)[i];
        }
    }

    //Генерируем данные для данного потока. В реальной задаче лучше сделать чтение из файлы и т.п.
    *matrix = new double[(*perProcess)[rank] * N];
    for(int i = 0; i < (*perProcess)[rank]; ++i) {
        for(int j = 0; j < N; ++j) {
            (*matrix)[i * N + j] = ((*startLine) + i) == j ? 2 : 1;
        }
    }

    *b = new double[(*perProcess)[rank]];
    for(int i = 0; i < (*perProcess)[rank]; ++i) {
        (*b)[i] = N + 1;
    }

    *x = new double[N]();
}
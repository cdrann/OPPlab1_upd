#include <iostream>
#include <cmath>
#include <ctime>
#include <mpi.h>
#include <cstdio>
#include "test2.h"

void init(int **perProcess, int *startLine, double **matrix, double **b, double **x, int size, int rank);

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank = 0, size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int startLine = 0;
    int *perProcess = 0;
    double *matrix = 0, *b = 0, *x = 0;
    init(&perProcess, &startLine, &matrix, &b, &x, size, rank);

   // printf("I'm %d from %d and my lines: %d-%d (%d)\n", rank, size, startLine,
     //      startLine + perProcess[rank], perProcess[rank]);

    // Считаем норму b
    double startTime = 0, normB = 0;
    if(rank != 0) {
        for(int i = 0; i < perProcess[rank]; ++i) {
            normB += b[i] * b[i];
        }
        MPI_Send(&normB, 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD); //отправляем b
    } else {
        startTime = MPI_Wtime();

        for(int i = 0; i < perProcess[rank]; ++i) {
            normB += b[i] * b[i];
        }

        for(int i = 1; i < size; ++i) {
            double tmp;
            MPI_Status status;
            MPI_Recv(&tmp, 1, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, &status); //получаем столбец
            normB += tmp;
        }

        normB = sqrt(normB);
    }

    // сумма строк
    double *tmpSum = new double[perProcess[rank]];
    double *tmpX = new double[N / size + 1]();
    int keepCalc = 1;
    while(keepCalc) {
        for(int i = 0; i < perProcess[rank]; ++i) {
            tmpSum[i] = 0;
        }

        // вычисление части и пересылка по процессам
        for(int i = 0, currentCrds = startLine; i < size; ++i) {

            for(int j = 0; j < perProcess[rank]; ++j) {
                for(int k = currentCrds, c = currentCrds + perProcess[i]; k < c; ++k) {
                    tmpSum[j] += matrix[j * N + k] * x[k - currentCrds];
                }
            }

            MPI_Status status;
            MPI_Sendrecv(x, N / size + 1, MPI_DOUBLE, (rank - 1 + size) % size, 0,
                         tmpX, N / size + 1, MPI_DOUBLE, (rank + 1) % size, 0, MPI_COMM_WORLD, &status);

            std::swap(x, tmpX);
            currentCrds = (currentCrds + perProcess[i]) % N;
        }

        double processAnswer = 0;
        for(int i = 0; i < perProcess[rank]; ++i) {
            tmpSum[i] -= b[i]; //Ax-b
            x[i] = x[i] - tmpSum[i] * t; // (Ax-b) * tau
            processAnswer += tmpSum[i] * tmpSum[i]; // || -- ||
        }

        if(rank != 0) {
            // частичную сумму координат в нулевой
            MPI_Send(&processAnswer, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
        } else {
            // частичная сумма нулевого
            double sum = processAnswer;
            for(int i  = 1; i < size; ++i) {
                MPI_Status status;
                double tmp;
                MPI_Recv(&tmp, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &status);

                // добавляем частичные суммы остальных
                sum += tmp;
            }
            sum = sqrt(sum);

            keepCalc = sum / normB >= eps;
        }

        MPI_Bcast(&keepCalc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    if(rank != 0) {
        // отправляем столбец значений
        MPI_Send(x, perProcess[rank], MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    } else {
        double *fullX = new double[N];
        for(int i = 0; i < perProcess[rank]; ++i) {
            // получаем столбец от нулевого
            fullX[i] = x[i];
        }

        for(int i = 1, currentLine = perProcess[rank]; i < size; ++i) {
            MPI_Status status;
            // получаем столбец от всех остальных, на нужную координату
            MPI_Recv(&fullX[currentLine], perProcess[i], MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
            currentLine += perProcess[i];
        }

        double endTime = MPI_Wtime();
        std::cout << "Size: " << size << ", time: " << (endTime - startTime) << std::endl;

        bool correctAnswer = true;
        for(int i = 0; i < N; ++i) {
            //std::cout << fullX[i] << " ";
            if(fabs(fabs(fullX[i]) - 1) >= eps) {
                correctAnswer = false;
                break;
            }
        }

        if(correctAnswer)
            std::cout << "Accepted." << std::endl;
        else
            std::cout << "WA." << std::endl;

        delete[] fullX;
    }

    delete[] tmpX;
    delete[] x;
    delete[] b;
    delete[] matrix;
    delete[] perProcess;
    MPI_Finalize();
    return 0;
}

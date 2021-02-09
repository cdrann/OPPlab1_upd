//
// Created by Admin on 15.06.2020.
//

#include <iostream>
#include <cmath>
#include <ctime>
#include <mpi.h>
#include "test1.h"

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank = 0, size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int startLine = 0;
    int *perProcess = 0;
    double *matrix = 0, *b = 0, *x = 0;
    init(&perProcess, &startLine, &matrix, &b, &x, size, rank); // *x = new double[N]();

    printf("I'm %d from %d and my lines: %d-%d (%d)\n", rank, size, startLine,
           startLine + perProcess[rank], perProcess[rank]);

    // Считаем норму b
    double normB = 0, startTime = 0;
    if(rank == 0) {
        startTime = MPI_Wtime();

        for(int i = 0; i < N; ++i) {
            normB += b[i] * b[i];
        }
        normB = sqrt(normB);
    }

    double *processX = new double[perProcess[rank]];
    int keepCalc = 1;
    while(keepCalc) {
        double processAnswer = 0;

        //Цикл по строкам
        for(int i = 0; i < perProcess[rank]; ++i) {
            double sum = 0;
            for(int j = 0; j < N; ++j) { //Ax
                sum += matrix[i * N + j] * x[j];
            }
            sum -= b[i]; //Ax - b
            processX[i] = x[i + startLine] - t * sum; //сохранили координату, которую вычисляли
            processAnswer += sum * sum; // ||Ax - b||
        }

        if(rank != 0) {
            //отправляем столбец значений
            MPI_Send(processX, perProcess[rank], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            //частичную сумму координат
            MPI_Send(&processAnswer, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        } else {
            //Сохраняем координаты в x
            for(int i = startLine, c = startLine + perProcess[rank]; i < c; ++i) {
                x[i] = processX[i - startLine];
            }

            // получаем результаты
            double sum = processAnswer; // от нулевого
            for(int i  = 1, currentLine = perProcess[rank]; i < size; ++i) {
                // от остальных
                MPI_Status status;
                MPI_Recv(&x[currentLine], perProcess[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status); //получаем столбец
                currentLine += perProcess[i];

                // получаем частичную сумму
                double tmp;
                MPI_Recv(&tmp, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
                sum += tmp;
            }
            sum = sqrt(sum);

            // условие завершения счета
            keepCalc = sum / normB >= eps;
        }
        // обновляем данные
        MPI_Bcast(&keepCalc, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    if(rank == 0) {
        double endTime = MPI_Wtime();
        std::cout << "Size: " << size << ", time: " << (endTime - startTime) << std::endl;

        bool correctAnswer = true;
        for(int i = 0; i < N; ++i) {
            //std::cout << x[i] << " ";
            if(fabs(fabs(x[i]) - 1) >= eps) {
                correctAnswer = false;
                break;
            }
        }
        if(correctAnswer)
            std::cout << "Accepted." << std::endl;
        else
            std::cout << "WA." << std::endl;
    }

    delete[] processX;
    delete[] x;
    delete[] b;
    delete[] matrix;
    delete[] perProcess;
    MPI_Finalize();
    return 0;
}

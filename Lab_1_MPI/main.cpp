#include<iostream>
#include<fstream>
using namespace std;
#include<mpi.h>

int main(int argc, char* argv[]) {

    double sendStart, sendStop, recvStart, recvStop;
    int* a, * b, * c, * buffer, * ans;
    int size = 1008;
    int rank, numprocs, line;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    line = size / numprocs;   //Divide matrix a by numprocs

    b = new int[size * size];
    ans = new int[size * line];

    if (rank == 0) {

        a = new int[size * size];
        c = new int[size * size];

        //Initialize matrix a & b
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++) {
                a[i * size + j] = i + j + 1;
                b[i * size + j] = 1;
            }

        sendStart = MPI_Wtime();

        for (int i = 1; i < numprocs; i++) {
            MPI_Send(b, size * size, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

        for (int i = 1; i < numprocs; i++) {
            MPI_Send(a + (i - 1) * line * size, size * line, MPI_INT, i, 1, MPI_COMM_WORLD);
        }

        sendStop = MPI_Wtime();

        //Thread 0 computes the last block
        for (int i = (numprocs - 1) * line; i < size; i++) {
            for (int j = 0; j < size; j++) {
                int temp = 0;
                for (int k = 0; k < size; k++)
                    temp += a[i * size + k] * b[k * size + j];
                c[i * size + j] = temp;
            }
        }

        recvStart = MPI_Wtime();

        for (int k = 1; k < numprocs; k++) {
            MPI_Recv(ans, line * size, MPI_INT, k, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < line; i++) {
                for (int j = 0; j < size; j++) {
                    c[((k - 1) * line + i) * size + j] = ans[i * size + j];
                }
            }
        }

        recvStop = MPI_Wtime();

        ofstream out("c.txt");
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++)
                out << c[i * size + j] << '\t';
            out << '\n';
        }
        out.close();

        printf("size: %d numprocs: %d send time: %lfs compute time: %lfs recv time: %lfs total time: %lfs\n"
            , size, numprocs, sendStop - sendStart, recvStart - sendStop, recvStop - recvStart, recvStop - sendStart);

        delete[] a, c;

    }
    else {
        buffer = new int[size * line];

        recvStart = MPI_Wtime();

        MPI_Recv(b, size * size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(buffer, size * line, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        recvStop = MPI_Wtime();

        for (int i = 0; i < line; i++)
            for (int j = 0; j < size; j++) {
                int temp = 0;
                for (int k = 0; k < size; k++)
                    temp += buffer[i * size + k] * b[k * size + j];
                ans[i * size + j] = temp;
            }

        sendStart = MPI_Wtime();

        MPI_Send(ans, line * size, MPI_INT, 0, 3, MPI_COMM_WORLD);

        sendStop = MPI_Wtime();

        printf("thread: %d send time: %lfs compute time: %lfs recv time: %lfs total time: %lfs\n"
            , rank, recvStop - recvStart, sendStart - recvStop, sendStop - sendStart, sendStop - recvStart);

        delete[] buffer;
        delete[] ans;
    }

    delete[] b;

    MPI_Finalize();

    return 0;
}
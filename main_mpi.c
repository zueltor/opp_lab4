#include <stdio.h>
#include <stdlib.h>
#include <mpi/mpi.h>
#include <math.h>

#define Nx 50
#define Ny 50
#define Nz 50

#define min_x (-1)
#define min_y (-1)
#define min_z (-1)
#define max_x (1)
#define max_y (1)
#define max_z (1)
#define a (1e5)
#define eps (1e-8)
#define x0 (min_x)
#define y0 (min_y)
#define z0 (min_z)

#define Dx (max_x-min_x)
#define Dy (max_y-min_y)
#define Dz (max_z-min_z)
#define hx ((double)(Dx)/(Nx-1))
#define hy ((double)(Dy)/(Ny-1))
#define hz ((double)(Dz)/(Nz-1))
#define x(i) (x0+(i)*hx)
#define y(j) (y0+(j)*hy)
#define z(k) (z0+(k)*hz)
#define IND(i, j, k) ((i)*Ny*Nz+(j)*Nz+(k))

const double hx2 = (hx * hx);
const double hy2 = (hy * hy);
const double hz2 = (hz * hz);
const double c = 1.0 / (2.0 / (hx * hx) + 2.0 / (hy * hy) + 2.0 / (hz * hz) + a);

double phi(double x, double y, double z) {
    return x * x + y * y + z * z;
}

double ro(double x, double y, double z) {
    return 6 - a * phi(x, y, z);
}

double next_phi2(int i, double *current_phi, double *prev_phi) {
    double delta_max = 0;
    double delta;
    double tmp1, tmp2, tmp3;

    int index;
    for (int j = 1; j < Ny - 1; j++) {
        for (int k = 1; k < Nz - 1; k++) {
            tmp1 = (prev_phi[IND(i + 1, j, k)] + prev_phi[IND(i - 1, j, k)]) / hx2;
            tmp2 = (prev_phi[IND(i, j + 1, k)] + prev_phi[IND(i, j - 1, k)]) / hy2;
            tmp3 = (prev_phi[IND(i, j, k + 1)] + prev_phi[IND(i, j, k - 1)]) / hz2;

            current_phi[IND(i, j, k)] = c * (tmp1 + tmp2 + tmp3 - ro(x(i), y(j), z(k)));

            delta = fabs(current_phi[IND(i, j, k)] - prev_phi[IND(i, j, k)]);
            if (delta_max < delta) {
                delta_max = delta;
            }

        }
    }

    for (int j = 1; j < Ny - 1; j++) {
        for (int k = 1; k < Nz - 1; k++) {
            index = IND(i, j, k);
            prev_phi[index] = current_phi[index];
        }
    }
    return delta_max;
}

double next_phi(double *current_phi, double *prev_phi) {
    double delta_max = 0;
    double delta;
    double tmp1, tmp2, tmp3;
    int center = Nx / 2;
    int shift = 0;

    int index;
    for (int i = center; i < Nx - 1 && i > 0; i += shift) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int k = 1; k < Nz - 1; k++) {
                tmp1 = (prev_phi[IND(i + 1, j, k)] + prev_phi[IND(i - 1, j, k)]) / hx2;
                tmp2 = (prev_phi[IND(i, j + 1, k)] + prev_phi[IND(i, j - 1, k)]) / hy2;
                tmp3 = (prev_phi[IND(i, j, k + 1)] + prev_phi[IND(i, j, k - 1)]) / hz2;

                current_phi[IND(i, j, k)] = c * (tmp1 + tmp2 + tmp3 - ro(x(i), y(j), z(k)));

                delta = fabs(current_phi[IND(i, j, k)] - prev_phi[IND(i, j, k)]);
                if (delta_max < delta) {
                    delta_max = delta;
                }

            }
        }
        shift = (shift == 0) ? -1 : (shift < 0) ? -shift : -shift - 1;
    }

    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int k = 1; k < Nz - 1; k++) {
                index = IND(i, j, k);
                prev_phi[index] = current_phi[index];
            }
        }
    }
    return delta_max;
}

double *init_offsets(double min, double delta, int size) {
    double *offsets = (double *) malloc(sizeof(double) * size);
    for (int i = 0; i < size; ++i) {
        offsets[i] = min + delta * i;
    }
    return offsets;

}

void init_matrix(double *current_phi, double *prev_phi, int Nx_size, int base_x) {
    int index;
    for (int i = base_x; i < Nx_size + base_x; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                index = IND(i-base_x, j, k);
                if (i == 0 || i == Nx - 1 || j == 0 || j == Ny - 1 || k == 0 || k == Nz - 1) {
                    prev_phi[index] = phi(x(i), y(j), z(k));
                } else {
                    prev_phi[index] = 0;
                }
                current_phi[index] = prev_phi[index];
            }
        }
    }
}

void print(double *current_phi) {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                printf("%f ", current_phi[IND(i, j, k)]);
            }
            printf("\n");
        }
        printf("\n\n");
    }
}

void print_max_delta(double *current_phi) {
    double delta;
    double max_delta = 0;
    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 20; j++) {
            for (int k = 0; k < 20; k++) {
                delta = fabs(current_phi[IND(i, j, k)] - phi(x(i), y(j), z(k)));
                if (max_delta < delta) {
                    max_delta = delta;
                }

            }
        }
    }
    printf("Max delta = %f\n", max_delta);
}

double get_max_delta(double *current_phi, int Nx_process) {
    double delta;
    double max_delta = 0;
    for (int i = 1; i <= Nx_process; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                delta = fabs(current_phi[IND(i, j, k)] - phi(x(i), y(j), z(k)));
                if (max_delta < delta) {
                    max_delta = delta;
                }

            }
        }
    }
    return max_delta;
}


double max(double x, double y) {
    return (x > y) ? x : y;
}

void sendAndRecieve(double *matrix, int rank, int size, int Nx_size, MPI_Request *rq) {
    if (rank != 0) {
        MPI_Isend(matrix + Ny * Nz, Ny * Nz, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &rq[0]);
        MPI_Irecv(matrix, Ny * Nz, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &rq[1]);
    }
    if (rank != size - 1) {
        MPI_Isend(matrix + (Nx_size) * Ny * Nz, Ny * Nz, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &rq[2]);
        MPI_Irecv(matrix + (Nx_size + 1) * Ny * Nz, Ny * Nz, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &rq[3]);
    }
}

void waitforSendRecv(int rank, int size, MPI_Request *rq) {
    if (rank != 0) {
        MPI_Wait(&rq[0], MPI_STATUS_IGNORE);
        MPI_Wait(&rq[2], MPI_STATUS_IGNORE);
    }
    if (rank != size - 1) {
        MPI_Wait(&rq[1], MPI_STATUS_IGNORE);
        MPI_Wait(&rq[3], MPI_STATUS_IGNORE);
    }
}

int main(int argc, char **argv) {

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (Nx % size != 0) {
        printf("Nx should be devisible by size\n");
        MPI_Finalize();
        return 0;
    }


    int Nx_process = (Nx / size) + 2;
    int Nx_size = Nx_process + 2;
    double *current_phi = (double *) malloc(sizeof(double) * (Nx_size) * Ny * Nz);
    double *prev_phi = (double *) malloc(sizeof(double) * (Nx_size) * Ny * Nz);
    int base_x = rank * Nx_size - 1;
    init_matrix(current_phi, prev_phi, Nx_size, base_x);
    int i = 0;
    double delta;
    double delta_max;
    double delta_max_total = 0;
    MPI_Request rq[4];


    while (1) {
        delta_max=0;
        //delta = next_phi2(1, current_phi, prev_phi);
        //delta_max=max(delta_max,delta);
        //delta = next_phi2(base_x, current_phi, prev_phi);
        //delta_max = max(delta_max, delta);
        sendAndRecieve(prev_phi, rank, size, Nx_size, rq);

        for (int j = 2; j < Nx_process; j++) {
            delta = next_phi2(j, current_phi, prev_phi);
            delta_max = max(delta, delta_max);
            printf("d %f\n",delta_max);
        }
        waitforSendRecv(rank, size, rq);

        delta = next_phi2(1, current_phi, prev_phi);
        delta_max = max(delta, delta_max);
        delta = next_phi2(Nx_process, current_phi, prev_phi);
        delta_max = max(delta, delta_max);

        MPI_Allreduce(&delta_max, &delta_max_total, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        if (rank == 0) {
            //printf("d=%f\n",delta_max_total);
            i++;
            //printf("%d\n",i);
        }

        if (delta_max_total < eps) {
            break;
        }
    }

    delta_max=get_max_delta(current_phi,Nx_process);
    MPI_Reduce(&delta_max,&delta_max_total,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if (rank==0) {
        printf("i=%d\n", i);
        printf("delta max = %f\n",delta_max_total);
    }
    //print(current_phi);
    print_max_delta(current_phi);

    MPI_Finalize();

    return 0;
}
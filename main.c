#include <stdio.h>
#include <stdlib.h>
#include <mpi/mpi.h>
#include <math.h>

#define Nx 200
#define Ny 200
#define Nz 200

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

#define next_offset(prev_offset)            ((prev_offset) == 0 \
                                                ? -1 \
                                                : ((prev_offset) < 0 \
                                                    ? -(prev_offset) \
                                                    : -(prev_offset) - 1))

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

void init_matrix(double *current_phi, double *prev_phi) {
    int index;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                index = IND(i, j, k);
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

int main(int argc, char **argv) {

    double *current_phi = (double *) malloc(sizeof(double) * Nx * Ny * Nz);
    double *prev_phi = (double *) malloc(sizeof(double) * Nx * Ny * Nz);
    init_matrix(current_phi, prev_phi);
    int i = 0;
    double delta = 0;

    do {
        delta = next_phi(current_phi, prev_phi);
        i++;
    } while (fabs(delta) > eps);


    printf("i=%d\n", i);
    //print(current_phi);
    print_max_delta(current_phi);

    return 0;
}
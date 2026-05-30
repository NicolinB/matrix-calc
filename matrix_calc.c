/*
 * matrix_calc.c
 * Matrix calculator: add, subtract, multiply, transpose,
 * determinant, and scalar multiplication.
 *
 * Compile: gcc matrix_calc.c -o matrix_calc -lm
 * Run:     ./matrix_calc
 *
 * Author: Barbaros Nicolin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_SIZE 8

typedef double Matrix[MAX_SIZE][MAX_SIZE];

/* ── Print ────────────────────────────────────────────────── */

void print_matrix(Matrix m, int rows, int cols, const char *label) {
    printf("\n  %s:\n", label);
    for (int i = 0; i < rows; i++) {
        printf("  |");
        for (int j = 0; j < cols; j++)
            printf(" %9.3f", m[i][j]);
        printf("  |\n");
    }
}

/* ── Input ────────────────────────────────────────────────── */

void read_matrix(Matrix m, int rows, int cols, const char *name) {
    printf("  Enter elements of %s (%dx%d):\n", name, rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++) {
            printf("    [%d][%d] = ", i + 1, j + 1);
            scanf("%lf", &m[i][j]);
        }
}

int get_size(int *rows, int *cols, const char *label) {
    printf("  %s — rows (1-%d): ", label, MAX_SIZE);
    scanf("%d", rows);
    printf("  %s — cols (1-%d): ", label, MAX_SIZE);
    scanf("%d", cols);
    if (*rows < 1 || *rows > MAX_SIZE || *cols < 1 || *cols > MAX_SIZE) {
        printf("  Invalid dimensions.\n");
        return 0;
    }
    return 1;
}

/* ── Operations ───────────────────────────────────────────── */

void mat_add(Matrix a, Matrix b, Matrix r, int rows, int cols) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            r[i][j] = a[i][j] + b[i][j];
}

void mat_sub(Matrix a, Matrix b, Matrix r, int rows, int cols) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            r[i][j] = a[i][j] - b[i][j];
}

/* A is (ra x ca), B is (ca x cb) */
void mat_mul(Matrix a, Matrix b, Matrix r, int ra, int ca, int cb) {
    for (int i = 0; i < ra; i++)
        for (int j = 0; j < cb; j++) {
            r[i][j] = 0.0;
            for (int k = 0; k < ca; k++)
                r[i][j] += a[i][k] * b[k][j];
        }
}

void mat_transpose(Matrix m, Matrix r, int rows, int cols) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            r[j][i] = m[i][j];
}

void mat_scalar(Matrix m, Matrix r, int rows, int cols, double s) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            r[i][j] = m[i][j] * s;
}

/* Recursive determinant (Laplace expansion along first row) */
double determinant(Matrix m, int n) {
    if (n == 1) return m[0][0];
    if (n == 2) return m[0][0]*m[1][1] - m[0][1]*m[1][0];

    double det = 0.0;
    Matrix sub;
    for (int col = 0; col < n; col++) {
        /* Build (n-1)x(n-1) sub-matrix without row 0 and column col */
        for (int i = 1; i < n; i++) {
            int sc = 0;
            for (int j = 0; j < n; j++) {
                if (j == col) continue;
                sub[i - 1][sc++] = m[i][j];
            }
        }
        double sign = (col % 2 == 0) ? 1.0 : -1.0;
        det += sign * m[0][col] * determinant(sub, n - 1);
    }
    return det;
}

/* ── Main ─────────────────────────────────────────────────── */

int main(void) {
    Matrix A, B, R;
    int choice;

    do {
        printf("\n=== Matrix Calculator ===\n");
        printf("  1. Add         (A + B)\n");
        printf("  2. Subtract    (A - B)\n");
        printf("  3. Multiply    (A x B)\n");
        printf("  4. Transpose   (A^T)\n");
        printf("  5. Determinant (square matrix)\n");
        printf("  6. Scalar multiply\n");
        printf("  0. Exit\n\nChoice: ");
        scanf("%d", &choice);

        int ra, ca, rb, cb;

        switch (choice) {
            case 1:
            case 2:
                if (!get_size(&ra, &ca, "Matrix A")) break;
                rb = ra; cb = ca;
                read_matrix(A, ra, ca, "A");
                read_matrix(B, rb, cb, "B");
                if (choice == 1) mat_add(A, B, R, ra, ca);
                else             mat_sub(A, B, R, ra, ca);
                print_matrix(A, ra, ca, "A");
                print_matrix(B, ra, ca, "B");
                print_matrix(R, ra, ca, choice == 1 ? "A + B" : "A - B");
                break;

            case 3:
                if (!get_size(&ra, &ca, "Matrix A")) break;
                rb = ca;
                printf("  Matrix B — cols (1-%d): ", MAX_SIZE);
                scanf("%d", &cb);
                if (cb < 1 || cb > MAX_SIZE) { printf("  Invalid.\n"); break; }
                read_matrix(A, ra, ca, "A");
                read_matrix(B, rb, cb, "B");
                mat_mul(A, B, R, ra, ca, cb);
                print_matrix(R, ra, cb, "A x B");
                break;

            case 4:
                if (!get_size(&ra, &ca, "Matrix")) break;
                read_matrix(A, ra, ca, "A");
                mat_transpose(A, R, ra, ca);
                print_matrix(A,  ra, ca, "A");
                print_matrix(R,  ca, ra, "Transpose(A)");
                break;

            case 5:
                printf("  Square matrix size (1-%d): ", MAX_SIZE);
                scanf("%d", &ra);
                if (ra < 1 || ra > MAX_SIZE) { printf("  Invalid.\n"); break; }
                read_matrix(A, ra, ra, "A");
                print_matrix(A, ra, ra, "A");
                printf("\n  det(A) = %.6f\n", determinant(A, ra));
                break;

            case 6: {
                if (!get_size(&ra, &ca, "Matrix")) break;
                read_matrix(A, ra, ca, "A");
                double s;
                printf("  Scalar value: ");
                scanf("%lf", &s);
                mat_scalar(A, R, ra, ca, s);
                print_matrix(R, ra, ca, "Result");
                break;
            }

            case 0:
                printf("Goodbye!\n");
                break;

            default:
                printf("  Invalid choice.\n");
        }
    } while (choice != 0);

    return 0;
}

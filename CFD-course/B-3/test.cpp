#include <iostream>
#include <vector>
#include <iomanip> // For precision
#include <cmath>   // For abs()

using namespace std;

// Constants
const int n_x = 4;  // Number of grid points in x-direction
const int n_y = 3;  // Number of grid points in y-direction
const double L_x = 1.0;  // Length in x-direction
const double L_y = 1.0;  // Length in y-direction
const double dx = L_x / (n_x - 1);  // Grid spacing in x
const double dy = L_y / (n_y - 1);  // Grid spacing in y
const double dt = 1.0;  // Time step
const double alpha = 1.0;  // Diffusion coefficient
const double v_x = 1.0, v_y = 1.0;  // Convection velocity components
const int t_steps = 2;  // Number of time steps

// Helper function to calculate the 1D index of a grid point
inline int index(int i, int j) {
    return j * n_x + i;
}

// Function to assemble the coefficient matrix A and right-hand side vector b
void assembleSystem(vector<vector<double>>& A, vector<double>& b, const vector<double>& u_prev) {
    int N = (n_x - 2) * (n_y - 2); // Number of internal points (exclude boundaries)

    // Coefficients for discretization
    double rx = alpha * dt / (dx * dx); // Diffusion in x-direction
    double ry = alpha * dt / (dy * dy); // Diffusion in y-direction
    double cx = v_x * dt / dx;          // Convection in x-direction
    double cy = v_y * dt / dy;          // Convection in y-direction

    // Loop over all internal points
    int eq = 0; // Equation index
    for (int j = 1; j < n_y - 1; ++j) {
        for (int i = 1; i < n_x - 1; ++i) {
            int p = index(i, j); // Index of the current point

            // Right-hand side (source term and previous u)
            b[eq] = u_prev[p] + dt;

            // Reset the row in matrix A
            for (int k = 0; k < N; ++k) {
                A[eq][k] = 0.0;
            }

            // Fill the coefficient matrix A
            A[eq][eq] = 1 + 2 * rx + 2 * ry; // Center coefficient
            if (i > 1) {
                A[eq][eq - 1] = -rx - cx; // Left neighbor
            }
            if (i < n_x - 2) {
                A[eq][eq + 1] = -rx; // Right neighbor
            }
            if (j > 1) {
                A[eq][eq - (n_x - 2)] = -ry - cy; // Bottom neighbor
            }
            if (j < n_y - 2) {
                A[eq][eq + (n_x - 2)] = -ry; // Top neighbor
            }

            eq++; // Move to the next equation
        }
    }
}

// Solve a linear system Ax = b using Gauss-Seidel iteration
void solveLinearSystem(vector<vector<double>>& A, vector<double>& b, vector<double>& x, int max_iter = 1000, double tol = 1e-6) {
    int N = b.size();
    for (int iter = 0; iter < max_iter; ++iter) {
        double error = 0.0;
        for (int i = 0; i < N; ++i) {
            double sigma = 0.0;
            for (int j = 0; j < N; ++j) {
                if (i != j) {
                    sigma += A[i][j] * x[j];
                }
            }
            double x_new = (b[i] - sigma) / A[i][i];
            error += fabs(x_new - x[i]);
            x[i] = x_new;
        }
        if (error < tol) {
            break;
        }
    }
}

int main() {
    // Total number of internal points (excluding boundaries)
    int N = (n_x - 2) * (n_y - 2);

    // Solution vectors
    vector<double> u_prev(n_x * n_y, 0.0); // Previous time step
    vector<double> u_next(n_x * n_y, 0.0); // Current time step

    // Matrix A and vector b for the linear system
    vector<vector<double>> A(N, vector<double>(N, 0.0));
    vector<double> b(N, 0.0);

    // Time-stepping loop
    for (int t = 0; t < t_steps; ++t) {
        // Assemble the linear system
        assembleSystem(A, b, u_prev);

        // Solve the linear system using Gauss-Seidel
        vector<double> x(N, 0.0); // Internal points solution
        solveLinearSystem(A, b, x);

        // Update the solution for the next time step
        int eq = 0;
        for (int j = 1; j < n_y - 1; ++j) {
            for (int i = 1; i < n_x - 1; ++i) {
                u_next[index(i, j)] = x[eq++];
            }
        }
        u_prev = u_next; // Move to the next time step
    }

    // Output results
    cout << "u values at t = 2:" << endl;
    cout << "u(6): " << fixed << setprecision(5) << u_next[index(1, 1)] << endl;
    cout << "u(7): " << fixed << setprecision(5) << u_next[index(2, 1)] << endl;

    return 0;
}


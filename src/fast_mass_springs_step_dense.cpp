#include "fast_mass_springs_step_dense.h"
#include <igl/matlab_format.h>
#include <iostream>

void fast_mass_springs_step_dense(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const double k,
    const Eigen::VectorXi& b,
    const double delta_t,
    const Eigen::MatrixXd& fext,
    const Eigen::VectorXd& r,
    const Eigen::MatrixXd& M,
    const Eigen::MatrixXd& A,
    const Eigen::MatrixXd& C,
    const Eigen::LLT<Eigen::MatrixXd>& prefactorization,
    const Eigen::MatrixXd& Uprev,
    const Eigen::MatrixXd& Ucur,
    Eigen::MatrixXd& Unext)
{
    Eigen::MatrixXd y1 = (1 / (delta_t * delta_t)) * M * (2 * Ucur - Uprev) + fext;
    Eigen::MatrixXd y2 = C.transpose() * C * 1e10 * V;
    Eigen::MatrixXd d(E.rows(), 3);
    Unext = Ucur;

    // Main solver loop (50 iterations)
    for (int iter = 0; iter < 50; iter++)
    {
        for (int i = 0; i < E.rows(); i++) {
            int v_start = E(i, 0);
            int v_end = E(i, 1);
            Eigen::RowVector3d edge = Unext.row(v_start) - Unext.row(v_end);
            d.row(i) = (edge.normalized()) * r(i);
        }
        const Eigen::MatrixXd l = k * A.transpose() * d + y1 + y2;
        Unext = prefactorization.solve(l);
    }

    // NEW: Floor Collision
	// Check each vertex for penetration below y=0, make sure the cubes stop at the floor
    const double floor_y = 0.0;
    const double restitution = 0.85;
    const double friction = 0.7;

    int collision_count = 0;

    for (int i = 0; i < Unext.rows(); i++) {
        if (Unext(i, 1) < floor_y) {
            collision_count++;

            double delta_y = Unext(i, 1) - Ucur(i, 1);
            Unext(i, 1) = floor_y - delta_y * restitution;

            if (Unext(i, 1) < floor_y) {
                Unext(i, 1) = floor_y;
            }

            double vx = (Ucur(i, 0) - Uprev(i, 0)) / delta_t;
            double vz = (Ucur(i, 2) - Uprev(i, 2)) / delta_t;
            Unext(i, 0) = Ucur(i, 0) + vx * friction * delta_t;
            Unext(i, 2) = Ucur(i, 2) + vz * friction * delta_t;
        }
    }

	// NEW: Cube-Cube Collision Handling
	// We have 3 cubes stacked vertically. We check for penetration of the upper cube into the lower cube.
    const double push_strength = 1.0;
    const double surface_thickness = 0.01;  // Collision zone
    const double damping = 0.1;              // Smooths top cubes jitter
    const int verts_per_cube = 8;            

    // Cube start indices in Unext
    int start1 = 0;   // bottom cube
    int start2 = 8;   // middle cube
    int start3 = 16;  // top cube

    int cube_starts[3] = { start1, start2, start3 };

    // Loop over cube pairs (bottom -> middle, middle -> top)
    for (int pair = 0; pair < 2; pair++) {
        int start_lower = cube_starts[pair];
        int start_upper = cube_starts[pair + 1];

        double top_y_lower = Unext.block(start_lower, 1, verts_per_cube, 1).maxCoeff();

        // Loop over upper cube vertices
        for (int i = 0; i < verts_per_cube; i++) {
            int idx_upper = start_upper + i;
            double vertex_y = Unext(idx_upper, 1);

            // Check if penetrating top surface
            if (vertex_y < top_y_lower + surface_thickness) {
                double x = Unext(idx_upper, 0);
                double z = Unext(idx_upper, 2);
                if (x >= -0.15 && x <= 0.15 && z >= -0.15 && z <= 0.15) {
                    double penetration = (top_y_lower + surface_thickness) - vertex_y;

                    if (penetration > 0) {
                        // Push upper cube up
                        Unext(idx_upper, 1) += penetration * push_strength;

                        // Compress lower cube vertices slightly, only if well above floor
                        for (int j = 0; j < verts_per_cube; j++) {
                            int idx_lower = start_lower + j;
                            if (Unext(idx_lower, 1) > 0.05) { // threshold above floor
                                Unext(idx_lower, 1) -= penetration * push_strength * 0.05; // small push
                            }
                        }
                    }
                }
            }
        }

        // Damping for upper cube to smooth jitter
        Unext.block(start_upper, 1, verts_per_cube, 1) =
            (1.0 - damping) * Unext.block(start_upper, 1, verts_per_cube, 1) +
            damping * Ucur.block(start_upper, 1, verts_per_cube, 1);
    }

}

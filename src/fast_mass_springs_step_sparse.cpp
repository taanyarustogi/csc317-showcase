#include "fast_mass_springs_step_sparse.h"
#include <igl/matlab_format.h>

void fast_mass_springs_step_sparse(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXi & b,
  const double delta_t,
  const Eigen::MatrixXd & fext,
  const Eigen::VectorXd & r,
  const Eigen::SparseMatrix<double>  & M,
  const Eigen::SparseMatrix<double>  & A,
  const Eigen::SparseMatrix<double>  & C,
  const Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > & prefactorization,
  const Eigen::MatrixXd & Uprev,
  const Eigen::MatrixXd & Ucur,
  Eigen::MatrixXd & Unext)
{
    Eigen::MatrixXd y1 = (1 / (delta_t * delta_t)) * M * (2 * Ucur - Uprev) + fext;
    Eigen::MatrixXd y2 = C.transpose() * C * 1e10 * V;

	Eigen::SparseMatrix<double> At = A.transpose();

    Unext = Ucur;

  for(int iter = 0;iter < 50;iter++)
  {
      Eigen::MatrixXd d(E.rows(), 3);

      for (int i = 0; i < E.rows(); i++) {
          int v_start = E(i, 0);
          int v_end = E(i, 1);
          Eigen::RowVector3d edge = Unext.row(v_start) - Unext.row(v_end);
          d.row(i) = (edge.normalized()) * r(i);
      }
    const Eigen::MatrixXd l = k * At * d + y1 + y2;
    Unext = prefactorization.solve(l);
  }
}

#include "fast_mass_springs_precomputation_dense.h"
#include "signed_incidence_matrix_dense.h"
#include <Eigen/Dense>

bool fast_mass_springs_precomputation_dense(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXd & m,
  const Eigen::VectorXi & b,
  const double delta_t,
  Eigen::VectorXd & r,
  Eigen::MatrixXd & M,
  Eigen::MatrixXd & A,
  Eigen::MatrixXd & C,
  Eigen::LLT<Eigen::MatrixXd> & prefactorization)
{
	r = Eigen::VectorXd::Zero(E.rows());
	for (int i = 0; i < E.rows(); i++) {
		int v_start = E(i, 0);
		int v_end = E(i, 1);
		r(i) = (V.row(v_start) - V.row(v_end)).norm();
	}
	M = Eigen::MatrixXd::Zero(V.rows(), V.rows());
	for (int i = 0; i < m.size(); i++) {
		M(i, i) = m(i);
	}
	signed_incidence_matrix_dense(V.rows(), E, A);
	C = Eigen::MatrixXd::Zero(b.size(), V.rows());
	for (int i = 0; i < b.size(); i++) {
		C(i, b(i)) = 1.0;
	}
	Eigen::MatrixXd Q1 = k*(A.transpose() * A) + M/(delta_t*delta_t);
	Eigen::MatrixXd Q2 = C.transpose() * C * 1e10;
	Eigen::MatrixXd Q = Q1 + Q2;
  prefactorization.compute(Q);
  return prefactorization.info() != Eigen::NumericalIssue;
}

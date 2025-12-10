#include "signed_incidence_matrix_dense.h"

void signed_incidence_matrix_dense(
  const int n,
  const Eigen::MatrixXi & E,
  Eigen::MatrixXd & A)
{
	A = Eigen::MatrixXd::Zero(E.rows(), n);
	for (int i = 0; i < E.rows(); i++) {
		int v_start = E(i, 0);
		int v_end = E(i, 1);
		A(i, v_start) = 1.0;
		A(i, v_end) = -1.0;
	}
}

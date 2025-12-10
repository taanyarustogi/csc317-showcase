#include "signed_incidence_matrix_sparse.h"
#include <vector>

void signed_incidence_matrix_sparse(
  const int n,
  const Eigen::MatrixXi & E,
  Eigen::SparseMatrix<double>  & A)
{
  std::vector<Eigen::Triplet<double> > ijv;
  A.resize(E.rows(),n);
  for (int i = 0; i < E.rows(); i++) {
	  int v_start = E(i, 0);
	  int v_end = E(i, 1);
	  ijv.emplace_back(i, v_start, 1.0);
	  ijv.emplace_back(i, v_end, -1.0);
  }
  A.setFromTriplets(ijv.begin(),ijv.end());
}

#include "fast_mass_springs_precomputation_sparse.h"
#include "signed_incidence_matrix_sparse.h"
#include <vector>

bool fast_mass_springs_precomputation_sparse(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXd & m,
  const Eigen::VectorXi & b,
  const double delta_t,
  Eigen::VectorXd & r,
  Eigen::SparseMatrix<double>  & M,
  Eigen::SparseMatrix<double>  & A,
  Eigen::SparseMatrix<double>  & C,
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > & prefactorization)
{
  std::vector<Eigen::Triplet<double> > ijv;
  std::vector<Eigen::Triplet<double> > ijv2;
  const int n = V.rows();
  r.resize(E.rows());
  r.setZero();
  for (int i = 0; i < E.rows(); i++) {
      int v_start = E(i, 0);
      int v_end = E(i, 1);
      r(i) = (V.row(v_start) - V.row(v_end)).norm();
  }
  M.resize(n, n);
  M.setZero();
  for (int i = 0; i < m.size(); i++) {
      ijv.emplace_back(i, i, m(i));
  }
  M.setFromTriplets(ijv.begin(), ijv.end());

  signed_incidence_matrix_sparse(V.rows(), E, A);
  C.resize(b.size(), n);
  C.setZero();
  for (int i = 0; i < b.size(); i++) {
      ijv2.emplace_back(i, b(i), 1.0);
  }
  C.setFromTriplets(ijv2.begin(), ijv2.end());

  Eigen::SparseMatrix<double> Q1 = k * (A.transpose() * A) + M / (delta_t * delta_t);
  Eigen::SparseMatrix<double> Q2 = C.transpose() * C * 1e10;
  Eigen::SparseMatrix<double> Q(n,n);
  Q = Q1 + Q2;
  prefactorization.compute(Q);
  return prefactorization.info() != Eigen::NumericalIssue;
}

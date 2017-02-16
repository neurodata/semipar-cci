using namespace arma;

class RFTensor {
 public:
  int n, r;
  double sigma2;

  int m_j;

  std::vector<cube> A_list;
  std::vector<mat> F_list;
  std::vector<mat> C_list;
  std::vector<cube> Psi_list;
  std::vector<cube> w_list;

  vec func_loglogit(vec x) {
    vec l = -log(1 + exp(-x));
    l(find(l < -1E6)).fill(-1E6);
    l(find(l > 1E6)).fill(-1E6);
    return l;
  }

  RFTensor(std::vector<cube> _A_list, int rank) {
    A_list = _A_list;

    Psi_list = A_list;  // deep copy to initialize Psi_list
    w_list = A_list;

    r = rank;
    sigma2 = 0.01;
    m_j = A_list.size();
    n = (A_list[0]).n_rows;

    mat avgA = zeros(n, n);
    int m = 0;

    for (auto it = A_list.begin(); it < A_list.end(); it++) {
      avgA += sum(*it, 2);
      m += it->n_slices;
    }

    avgA /= m;

    mat avgF;
    vec avgC;
    mat V;
    svd(avgF, avgC, V, avgA);

    avgF.resize(n, r);
    avgC.resize(r);

    // initialize F and C
    for (int i = 0; i < m_j; ++i) {
      int p = (A_list[i]).n_slices;
      F_list.push_back(avgF + randn(n, r) * sqrt(sigma2));
      C_list.push_back(repmat(avgC, 1, p) + randn(r, p) * 0.1);
    }

    computePsi();

    for (int i = 0; i < m_j; ++i) {
      cout << Psi_list[i].slice(0) << endl;
    }
  }

  vec samplePsi(mat X, vec w, vec K, mat Binv, vec b) {
    mat V = inv(X.t() * (X.each_col() % w) + Binv);
    mat m = V * (X.t() * K + Binv * b);
    return m;
  }

  void computePsi() {
    for (int j = 0; j < m_j; ++j) {
      int p = (A_list[j]).n_slices;
      mat F_local = F_list[j];
      for (int i = 0; i < p; ++i) {
        Psi_list[j].slice(i) =
            F_local * trans(F_local.each_row() % trans(C_list[j].col(i)));
      }
    }
  };

  void updateW() {
    for (int j = 0; j < m_j; ++j) {
      cube b = Psi_list[j];
      w_list[j] = 0.5 / b % tanh(b / 2);
    }
  }

  void updateF_list() {
    vec idx = linspace<vec>(0, n - 1, n);
    for (int j = 0; j < m_j; ++j) {
      int p = (A_list[j]).n_slices;

      for (int k = 0; k < n; ++k) {
        mat F_wo_k = F_list[j].rows(find(idx != k));

        for (int i = 0; i < p; ++i) {
          vec y_i = A_list[j].slices(i)
        }
      }

      cube b = Psi_list[j];
      w_list[j] = 0.5 / b % tanh(b / 2);
    }
  }

  void runEM(int tot_iter) { updateW(); }
};
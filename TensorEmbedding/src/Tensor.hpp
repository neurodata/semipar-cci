using namespace arma;

enum LossType { sqr, logistic };

// symmetric 3D tensor
// n x n x p
// latent factoring dimension k
class SymmTensor {
 public:
  int n, p, k;
  cube A;

  // loss type
  LossType loss_pick;
  bool coreDiag; //restrict each layer of core to be diagonal

  // A =~ L C L
  mat L;   // left matrix n,k
  cube C;  // core tensor k,k,p

  // direction for L and C to drop
  mat directionL;
  cube directionC;

  SymmTensor(cube _A, LossType _loss_pick = sqr, bool _coreDiag = false) {
    loss_pick = _loss_pick;

    coreDiag = _coreDiag;

    A = _A;
    n = A.n_rows;
    p = A.n_slices;
  }

  // set latent dimension
  void setK(int _k) {
    k = _k;

    mat Amean = sum(A, 2) / p;

    vec s;
    mat V;
    svd(L, s, V, Amean);

    // L = trans(chol(Amean));

    L.resize(n, k);
    s.resize(k);

    C = zeros(k, k, p);
    for (int i = 0; i < p; ++i) {
      C.slice(i) = diagmat(s);  // making C symmetric
    }
  }

  // loss functions and its derivative wrt LCL:

  // squared loss
  double sqrLoss(mat L, cube C) {
    cube diff = zeros(n, n, p);

    for (int i = 0; i < p; ++i) {
      mat C_local = C.slice(i);
      diff.slice(i) = L * C_local * L.t() - A.slice(i);
    }

    return accu(diff % diff);
  }

  cube deriSqrLoss(mat L, cube C) {
    cube diff = zeros(n, n, p);

    for (int i = 0; i < p; ++i) {
      mat C_local = C.slice(i);
      diff.slice(i) = L * C_local * L.t() - A.slice(i);
    }

    return 2 * diff;
  }

  // logistic loss

  double logisticLoss(mat L, cube C) {
    cube theta = zeros(n, n, p);
    for (int i = 0; i < p; ++i) {
      mat C_local = C.slice(i);
      theta.slice(i) = L * C_local * L.t();
    }

    cube loss = -theta % A + log(1 + exp(theta));
    return accu(loss);
  }

  cube deriLogisticLoss(mat L, cube C) {
    cube theta = zeros(n, n, p);

    for (int i = 0; i < p; ++i) {
      mat C_local = C.slice(i);
      theta.slice(i) = L * C_local * L.t();
    }

    return -A + 1.0 / (1.0 + exp(-theta));
  }

  // gradient for L
  mat gradL(mat L, cube C) {
    cube grad(n, k, p);

    cube grad_LCL;
    switch (loss_pick) {
      case sqr:
        grad_LCL = deriSqrLoss(L, C);
      case logistic:
        grad_LCL = deriLogisticLoss(L, C);
    }

    for (int i = 0; i < p; ++i) {
      mat LC = L * C.slice(i);
      mat diff = LC * L.t() - A.slice(i);
      grad.slice(i) = 2 * grad_LCL.slice(i) * LC;
      // grad.slice(i) = 4 *diff * LC;
    }
    return sum(grad, 2);  // sum over last dim
  }

  // gradient for C

  cube gradC(mat L, cube C) {
    cube grad= zeros(k, k, p);

    cube grad_LCL;
    switch (loss_pick) {
      case sqr:
        grad_LCL = deriSqrLoss(L, C);
      case logistic:
        grad_LCL = deriLogisticLoss(L, C);
    }

    for (int i = 0; i < p; ++i) {
      mat LC = L * C.slice(i);
      mat diff = LC * L.t() - A.slice(i);

      if (coreDiag) {
        for (int a = 0; a < k; ++a) {
          mat deriC = zeros(k, k);
          deriC(a, a) = 1;
          grad(a, a, i) = accu(grad_LCL.slice(i) % (L * deriC * L.t()));
        }
      }
      else {
        for (int a = 0; a < k; ++a) {
          for (int b = 0; b <= a; ++b) {
            mat deriC = zeros(k, k);
            deriC(a, b) = 1;
            deriC(b, a) = 1;
            grad(a, b, i) = accu(grad_LCL.slice(i) % (L * deriC * L.t()));
            grad(b, a, i) = grad(a, b, i);
          }
        }
      }
    }
    return grad;
  }

  // compute loss when move delta
  double lossAtDelta(double delta, int paramIdx = 0) {
    double cur_loss = sqrLoss(L, C);
    double result = 0;
    if (paramIdx == 0) {
      mat prop = L - delta * directionL;
      result = sqrLoss(prop, C);
    }
    if (paramIdx == 1) {
      cube prop = C - delta * directionC;
      result = sqrLoss(L, prop);
    }
    return result;
  }

  double lineSearch(double l, double r, double loss_l, double loss_r,
                    int paramIdx = 0) {
    if (fabs(loss_l - loss_r) < 1E-3) return l;

    double m = (l + r) / 2;

    double loss_m = lossAtDelta(m, paramIdx);

    if (loss_l < loss_m & loss_l < loss_r) return l;

    if (loss_r < loss_m & loss_r < loss_l) return r;

    if (loss_m < loss_l & loss_m < loss_r) {
      double opt_l = lineSearch(l, m, loss_l, loss_m, paramIdx);
      double opt_r = lineSearch(m, r, loss_m, loss_r, paramIdx);

      double loss_opt_l = lossAtDelta(opt_l, paramIdx);
      double loss_opt_r = lossAtDelta(opt_r, paramIdx);
      if (loss_opt_l < loss_opt_r)
        return opt_l;
      else
        return opt_r;
    }
    return l;
  }

  void OptConjugateGradient(int steps = 100, double delta1 = 1E-3,
                            double delta2 = 1E-3, double tol = 1E-8) {
    mat gradientL = gradL(L, C);
    cube gradientC = gradC(L, C);

    directionL = gradientL;
    directionC = gradientC;

    double pre_loss = INFINITY;

    for (int i = 0; i < steps; ++i) {
      {
        mat gradient0 = gradientL;
        gradientL = gradL(L, C);
        // Polak–Ribière conjugate gradient
        double beta = accu(gradientL % (gradientL - gradient0)) /
                      accu(gradient0 % gradient0);

        if (beta < 0) beta = 0;

        if (!std::isnan(beta)) directionL = gradientL + beta * directionL;

        double cur_loss = sqrLoss(L, C);
        double loss_delta = lossAtDelta(delta1, 0);
        double opt_delta = lineSearch(0, delta1, cur_loss, loss_delta, 0);
        L -= opt_delta * directionL;
      }

      {
        cube gradient0 = gradientC;
        gradientC = gradC(L, C);
        double beta = accu(gradientC % gradientC) / accu(gradient0 % gradient0);
        if (!std::isnan(beta)) directionC = gradientC + beta * directionC;

        double cur_loss = sqrLoss(L, C);
        double loss_delta = lossAtDelta(delta2, 1);
        double opt_delta = lineSearch(0, delta2, cur_loss, loss_delta, 1);
        C -= opt_delta * directionC;
      }
      {
        double cur_loss = sqrLoss(L, C);
        // if ((abs(directionL)).max() + (abs(directionC)).max() < 1E-5)
        if (fabs((pre_loss - cur_loss) / cur_loss) < tol)
          break;
        else
          pre_loss = cur_loss;

        cout << cur_loss << endl;
      }
    }
  }
  // = 1E-3)

  // // simple gradient descent (which sucks)
  // void gradientDescent(int steps = 1E4, double _delta1 = 1E-3, double _delta2
  // = 1E-3) {
  //   double cur_loss = sqrLoss(L, C);

  //   for (int i = 0; i < steps; ++i) {
  //       double delta1 = _delta1;
  //       double delta2 = _delta2;

  //     {
  //       double new_loss =  lineSearchAux(delta1,0);

  //       if (new_loss < cur_loss) {
  //         L -= delta1 * gradL(L, C);
  //         cur_loss = new_loss;
  //       } else {
  //         delta1 = delta1 / 10.0;
  //       }
  //     }

  //     {
  //       double new_loss =  lineSearchAux(delta2,1);
  //       if (new_loss < cur_loss) {
  //         C -= delta2 * gradC(L, C);
  //         cur_loss = new_loss;
  //       } else {
  //         delta2 = delta2 / 10.0;
  //       }
  //     }
  //     cout << cur_loss << endl;
  //   }
  // }

  // end of class definition
};
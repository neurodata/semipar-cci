using namespace arma;

// symmetric 3D tensor
// n x n x p
// latent factoring dimension k
class SymmTensor {
 public:
  int n, p, k;
  cube A;

  // A =~ L C L
  mat L;   // left matrix n,k
  cube C;  // core tensor k,k,p
  SymmTensor(cube _A) {
    A = _A;

    n = A.n_rows;
    p = A.n_slices;
  }

  //set latent dimension
  void setK(int _k) {
    k = _k;

    L = randn(n, k);

    C = randn(k, k, p);
    for (int i = 0; i < p; ++i) {
      mat C_local = C.slice(i);
      C.slice(i) = (C_local + trans(C_local)) / 2.0;  // making C symmetric
    }
  }

  //gradient for L
  mat gradL(mat L, cube C){
    cube grad(n, k, p);
    for (int i = 0; i < p; ++i)
    {
        mat LC = L * C.slice(i);
        mat diff =  LC* L.t() - A.slice(i);
        grad.slice(i) = 4* diff * LC;
    }
    return sum(grad, 2); //sum over last dim
  }

  //gradient for C

  cube gradC(mat L, cube C){
    cube grad(k, k, p);

    for (int i = 0; i < p; ++i)
    {
        mat LC = L * C.slice(i);
        mat diff =  LC* L.t() - A.slice(i);
        for (int a = 0; a < k; ++a)
        {
            for (int b = 0; b <= a; ++b)
            {
                mat deriC = zeros(k,k);
                deriC(a,b) =1;
                deriC(b,a) =1;
                grad(a,b,i) =  accu( 2*diff % (L * deriC * L.t()));
                 grad(b,a,i) = grad(a,b,i) ;
            }   
        }
    }
    return grad;
  }


};
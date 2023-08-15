#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
// [[Rcpp::depends("RcppArmadillo")]]


double ker(const rowvec &j1, const rowvec &j2, const rowvec &classes) {
    double res = 0.0;
    for (int i = 0; i < classes.n_elem; i++) {
        if (j1[i] == j2[i]) {
            res += sqrt(classes[i]);
        }
    }
    return res;
}

mat Kx_mat(const mat &xtrain, const rowvec &classes) {
    int ntrain = xtrain.n_rows;
    mat K = zeros<mat>(ntrain, ntrain);
    for (int i = 0; i < ntrain; i++) {
        for (int j = 0; j < ntrain; j++) {
            K(i, j) = ker(xtrain.row(i), xtrain.row(j), classes);
        }
    }
    return K / ntrain;
}

vec proxCpp(vec v, double lambda) {
    vec res;
    res.copy_size(v);
    res.zeros();
    double v_norm = sqrt(accu(v % v));
    if (v_norm >= lambda) {
        res = (1 - lambda / v_norm) * v;
    }
    return res;
}

int findIndex(const rowvec &x, const rowvec &classes) {
    rowvec dec = cumprod(classes);
    int res = x(0);
    for (int j = 1; j < x.n_elem; j++) {
        res += x(j) * dec(j - 1);
    }
    return res;
}

List getQcpp(const mat &x, const rowvec &classes) {
    int n = x.n_rows;
    int ntmp = prod(classes);
    mat Q = zeros<mat>(n, ntmp);
    mat Xtmp(ntmp, classes.n_elem);
    Xtmp.fill(-1);
    
    for (int i = 0; i < n; i++) {
        int ind = findIndex(x.row(i), classes);
        Q(i, ind) = 1;
        Xtmp.row(ind) = x.row(i);
    }
    
    uvec rows_to_remove = find(all(Xtmp == -1, 1));
    Xtmp.shed_rows(rows_to_remove);
    mat K = Kx_mat(Xtmp, classes);
    
    uvec zero_cols = find(sum(Q, 0) == 0);
    if (zero_cols.n_elem > 0) {
        Q.shed_cols(zero_cols);
    }
    
    return List::create(
        Named("Q") = Q,
        Named("K") = K,
        Named("Xtilde") = Xtmp
    );
}

// [[Rcpp::export]]
double eval_obj(const mat &y, const mat &alpha, const mat &omega, const mat &K, const mat &Q, const double &gamma, const double &eta) {
    int n = y.n_rows;
    int p = y.n_cols;
    int ntilde = K.n_rows;

    double res = 0;    
    res += 1.0 / n * trace((y - sqrt(ntilde) * Q * K * alpha) * omega * trans(y - sqrt(ntilde) * Q * K * alpha));
    
    res -= log(det(omega));

    double norm21 = 0;
    for (int j = 0; j < p; j++) {
        norm21 += sqrt(accu(square(alpha.col(j))));
    }
    res += gamma * norm21;
    if (eta != 0) {
        res += eta * accu(square(omega)) / 2;
    }
    return res;
}

mat get_grad(const mat &y, const mat &alpha, const mat &K, const mat &Q,
            const mat &omega) {
    int n = y.n_rows;
    int p = y.n_cols;
    int ntilde = K.n_rows;

    mat res = -2 * sqrt(ntilde) / n * K * trans(Q) * y * omega + 2 * K * trans(Q) * Q * K * alpha * omega * ntilde / n;
    return res;
}


double eval_f(const mat &y, const mat &alpha, const mat &omega,
                 const mat &K, const mat &Q) {
    int n = y.n_rows;
    int ntilde = K.n_rows;
    double res = 0;
    
    res = -2 * sqrt(ntilde) / n * trace(Q * K * alpha * omega * trans(y)) + 
         trace(Q * K * alpha * omega * trans(Q * K * alpha)) * ntilde / n;
    return res;
}

double eval_f_hat(const mat &y, const mat &alpha, const mat &alpha0, const mat &omega,
                 const mat &K, const mat &Q, const double &stepsize) {
    double res = 0;
    res = eval_f(y, alpha0, omega, K, Q) + trace(trans(get_grad(y, alpha0, K, Q, omega)) * (alpha - alpha0)) + 1 / (2 * stepsize) * accu(square(alpha - alpha0));
    return res;
}


// [[Rcpp::export]]
List update_alpha(const mat& y, const mat& K, const mat& Q, const mat& omega, double gamma, double eta, double stepsize = 2, const int& max_iter = 500, Nullable<NumericMatrix> alphaPre = R_NilValue) {

    int p = y.n_cols;
    int ntilde = K.n_rows;

    mat alpha;
    if (!alphaPre.isNotNull()) {
        alpha.zeros(ntilde, p);
    } else {
        alpha = as<mat>(alphaPre);
    }

    alpha.zeros(ntilde, p);

    double scaling = 0.5;
    double objPre = datum::inf;
    mat alphaPre2;
    if (!alphaPre.isNotNull()){
        alphaPre2.zeros(ntilde, p);
    } else {
        alphaPre2 = as<mat>(alphaPre);
    }
  
    mat grad_f(ntilde, p, fill::zeros);
    double obj = eval_obj(y, alpha, omega, K, Q, gamma, eta);
    int num_iter = 0;

    // Accelerated proximal gradient method
    while (std::abs(objPre - obj) > 1e-4 && num_iter <= max_iter) {
        mat alpha0 = alpha +  (alpha - alphaPre2) * (num_iter / (num_iter + 3.0));
        alphaPre2 = alpha;
        objPre = obj;
        grad_f = get_grad(y, alpha0, K, Q, omega);
        mat z;
        while (true) {
            z.zeros(ntilde, p);
            for (int j = 0; j < p; j++) {
                z.col(j) = proxCpp(alpha0.col(j) - stepsize * grad_f.col(j), stepsize * gamma);
            }
            if (eval_f(y, z, omega, K, Q) <= eval_f_hat(y, z, alpha0, omega, K, Q, stepsize)) {
                break;
            }
            stepsize = stepsize * scaling;
        }
        alpha = z;
        obj = eval_obj(y, alpha, omega, K, Q, gamma, eta);
//         Rcout << "Eval: " << obj << endl;
        num_iter++;
    }
    
//     Rcout << "Number of iterations when updating alpha:  " << num_iter << endl;
    return List::create(Named("alpha") = alpha,
                      Named("stepsize") = stepsize);
}






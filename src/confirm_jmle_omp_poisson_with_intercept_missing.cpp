#include <RcppArmadillo.h>

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

static arma::vec prox_func_cpp(const arma::vec &y, double C){
  double y_norm2 = arma::accu(square(y));
  if(y_norm2 <= C*C){
    return y;
  }
  else{
    return sqrt(C*C / y_norm2) * y;
  }
}


arma::vec prox_func_cpp2(const arma::vec &y, double C){
  double y_norm2 = arma::accu(square(y));
  if(y_norm2 <= C*C){
    return y;
  }
  else{
    arma::vec x = y * sqrt((C*C-1) / (y_norm2-1));
    x(0) =1;
    return x;
  }
}

double neg_loglik_poisson(const arma::mat &thetaA, const arma::mat &response, const arma::mat &nonmis_ind){
  int N = response.n_rows;
  int J = response.n_cols;
  return arma::accu(nonmis_ind % (arma::exp(thetaA) - response % thetaA))/ N / J ;
}

double neg_loglik_poisson_i_cpp(const arma::vec &response_i, const arma::vec &nonmis_ind_i,
                                const arma::mat &A, const arma::vec &theta_i){
  arma::vec tmp = A * theta_i;
  int J = response_i.n_elem;
  return arma::accu(nonmis_ind_i % (arma::exp(tmp) - response_i % tmp)) / J;
}

arma::vec grad_neg_loglik_poisson_thetai_cpp(const arma::vec &response_i, const arma::vec &nonmis_ind_i,
                                             const arma::mat &A, const arma::vec &theta_i){
  arma::vec tmp = A * theta_i;
  int J = response_i.n_elem;
  return A.t() * (nonmis_ind_i % (arma::exp(tmp) - response_i)) / J;
}

arma::mat Update_theta_poisson_cpp(const arma::mat &theta0, const arma::mat &response, const arma::mat &nonmis_ind,
                                   const arma::mat &A0, double C){
  arma::mat theta1 = theta0.t();
  int N = response.n_rows;
#pragma omp parallel for
  for(int i=0;i<N;++i){
    double step = 10;
    arma::vec h = grad_neg_loglik_poisson_thetai_cpp(response.row(i).t(), nonmis_ind.row(i).t(), A0, theta0.row(i).t());
    h(0) = 0;
    theta1.col(i) = theta0.row(i).t() - step * h;
    theta1.col(i) = prox_func_cpp2(theta1.col(i), C);
    while(neg_loglik_poisson_i_cpp(response.row(i).t(), nonmis_ind.row(i).t(), A0, theta1.col(i)) >
            neg_loglik_poisson_i_cpp(response.row(i).t(), nonmis_ind.row(i).t(), A0, theta0.row(i).t()) &&
          step > 1e-3){
      step *= 0.5;
      theta1.col(i) = theta0.row(i).t() - step * h;
      theta1.col(i) = prox_func_cpp2(theta1.col(i), C);
    }
  }
  return(theta1.t());
}

double neg_loglik_poisson_j_cpp(const arma::vec &response_j, const arma::vec &nonmis_ind_j,
                                const arma::vec &A_j, const arma::mat &theta){
  arma::vec tmp = theta * A_j;
  int N = response_j.n_elem;
  return arma::accu(nonmis_ind_j % (arma::exp(tmp) - response_j % tmp)) / N;
}

arma::vec grad_neg_loglik_poisson_A_j_cpp(const arma::vec &response_j, const arma::vec &nonmis_ind_j,
                                          const arma::vec &A_j, const arma::vec &Q_j,
                                          const arma::mat &theta){
  arma::vec tmp = theta * A_j;
  int N = response_j.n_elem;
  return (theta.t() * (nonmis_ind_j % (arma::exp(tmp) - response_j)) / N) % Q_j;
}
arma::mat Update_A_poisson_cpp(const arma::mat &A0, const arma::mat &Q, const arma::mat &response,
                               const arma::mat &nonmis_ind, const arma::mat &theta1, double C){
  arma::mat A1 = A0.t();
  int J = A0.n_rows;
#pragma omp parallel for
  for(int j=0;j<J;++j){
    double step = 10;
    arma::vec h = grad_neg_loglik_poisson_A_j_cpp(response.col(j), nonmis_ind.col(j), A0.row(j).t(), Q.row(j).t(), theta1);
    A1.col(j) = A0.row(j).t() - step * h;
    A1.col(j) = prox_func_cpp(A1.col(j), C);
    while(neg_loglik_poisson_j_cpp(response.col(j), nonmis_ind.col(j), A1.col(j), theta1) >
            neg_loglik_poisson_j_cpp(response.col(j), nonmis_ind.col(j), A0.row(j).t(), theta1) &&
          step > 1e-3){
      step *= 0.5;
      A1.col(j) = A0.row(j).t() - step * h;
      A1.col(j) = prox_func_cpp(A1.col(j), C);
    }
  }
  return(A1.t());
}

// [[Rcpp::export]]
Rcpp::List confirm_CJMLE_poisson_cpp(const arma::mat &response, const arma::mat &nonmis_ind,
                                     arma::mat theta0, arma::mat A0, const arma::mat &Q,
                                     double C, double tol=1e-4){
  arma::mat theta1 = Update_theta_poisson_cpp(theta0, response, nonmis_ind, A0, C);
  arma::mat A1 = Update_A_poisson_cpp(A0, Q, response, nonmis_ind, theta1, C);
  double eps = neg_loglik_poisson(theta0*A0.t(), response, nonmis_ind) - neg_loglik_poisson(theta1*A1.t(), response, nonmis_ind);
  while(eps > tol){
    theta0 = theta1;
    A0 = A1;
    theta1 = Update_theta_poisson_cpp(theta0, response, nonmis_ind, A0, C);
    A1 = Update_A_poisson_cpp(A0, Q, response, nonmis_ind, theta1, C);
    eps = neg_loglik_poisson(theta0*A0.t(), response, nonmis_ind) - neg_loglik_poisson(theta1*A1.t(), response, nonmis_ind);
  }
  return Rcpp::List::create(Rcpp::Named("A") = A1,
                            Rcpp::Named("theta") = theta1,
                            Rcpp::Named("obj") = neg_loglik_poisson(theta1*A1.t(), response, nonmis_ind));
}


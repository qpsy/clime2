#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List linprogPD2C(vec x, mat A, vec b,
									 double epsilon, double pdtol=1e-3, int pdmaxiter=50) {

	vec Atr = A*x - b;

	if (max(abs(Atr)) > epsilon) {
		Rcpp::stop("Infeasible starting point!");
	}

	double N = x.size();
  double alpha = 0.01;
  double beta = 0.5;
  double mu = 10;

	vec u = 0.95*abs(x) + 0.1*max(abs(x));
	vec fu1 = x - u;
  vec fu2 = -x - u;
  vec fe1 = Atr - epsilon;
  vec fe2 = -Atr - epsilon;
  vec lamu1 = -1/fu1;
  vec lamu2 = -1/fu2;
  vec lame1 = -1/fe1;
  vec lame2 = -1/fe2;

	mat At = trans(A);
	mat AtAv = At * (lame1 - lame2);

	vec flu1 = fu1 % lamu1;
	vec flu2 = fu2 % lamu2;
	vec fle1 = fe1 % lame1;
	vec fle2 = fe2 % lame2;
  double sdg = sum(flu1) + sum(flu2) + sum(fle1) + sum(fle2);
	double tau = mu*4*N/sdg;

	vec gradf0 = join_cols(zeros(N), ones(N));
	vec rdual = gradf0 + join_cols(lamu1 - lamu2 + AtAv, -lamu1 - lamu2);
	vec rcent = join_cols(join_cols(flu1, flu2), join_cols(fle1, fle2)) - 1/tau;
	double resnorm = sqrt(accu(join_cols(pow(rdual, 2), pow(rcent, 2))));

	int pditer = 0;
	bool done = (sdg < pdtol) || (pditer >= pdmaxiter)
	//while (!done) {
		vec w2 = -1 - (1/fu1 + 1/fu2)/tau;

    vec sig11 = -lamu1/fu1 - lamu2/fu2;
    vec sig12 = lamu1/fu1 - lamu2/fu2;
    vec siga = -lame1/fe1 - lame2/fe2;
    vec sigx = sig11 - pow(sig12, 2)/sig11;

    vec w1 = -(At*(1/fe2 - 1/fe1) + 1/fu2 - 1/fu1)/tau;
    vec w1p = w1 - (sig12/sig11) % w2;
		mat Hp = At*diagmat(siga)*A + diagmat(sigx);
		vec dx = solve(Hp, w1p);

		if (1/cond(Hp) < 1e-14) {
			Rcpp::warning("Ill conditioned matrix.  Previous iterate matrix returned! (May increase perturb/lambda.)");
			return wrap(x);
		}

		vec AtAdx = A*dx;

		vec du = w2/sig11 - (sig12/sig11)%dx;

    vec dlamu1 = -(lamu1/fu1)%(dx-du) - lamu1 - 1/(fu1*tau);
    vec dlamu2 = -(lamu2/fu2)%(-dx-du) - lamu2 - 1/(fu2*tau);

    vec dlame1 = -(lame1/fe1)%(AtAdx) - lame1 - 1/(fe1*tau);
    vec dlame2 = (lame2/fe2)%(AtAdx) - lame2 - 1/(fe2*tau);

    vec AtAdv = At*(dlame1 - dlame2);


    uvec iu1 = find(dlamu1 < 0);
    uvec iu2 = find(dlamu2 < 0);
    uvec ie1 = find(dlame1 < 0);
		uvec ie2 = find(dlame2 < 0);
		uvec ifu1 = find((dx-du) > 0);
		uvec ifu2 = find((-dx-du) > 0);
		uvec ife1 = find(AtAdx > 0);
		uvec ife2 = find(AtAdx < 0);
		vec aa = AtAdx.elem(ife1);

		//double smax = min( -lamu1[iu1]/dlamu1[iu1], -lamu2[iu2]/dlamu2[iu2], -lame1[ie1]/dlame1[ie1], -lame2[ie2]/dlame2[ie2], -fu1[ifu1]/(dx[ifu1] - du[ifu1]), -fu2[ifu2]/(-dx[ifu2] -du[ifu2]), -fe1[ife1]/AtAdx[ife1], -fe2[ife2]/( - AtAdx[ife2])   );
    //double smax = min(1, smax);
    //double s = 0.99*smax;
		//}
	return Rcpp::List::create(Rcpp::Named("x0") = AtAdx,
														Rcpp::Named("A") = aa,
														Rcpp::Named("b") = ife1
		);
}


/*** R
# library(microbenchmark)
set.seed(2)
X <- matrix(rnorm(300), ncol=3)
emat <- diag(3)
lam <- 1
Sigma <- var(X)
Omega0 <- solve(Sigma)
linprogPD2C(Omega0[,1], Sigma, emat[,1], lam)
*/

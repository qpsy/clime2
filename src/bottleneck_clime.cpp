#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List linprogPD2C(arma::vec x, arma::mat A, arma::vec b,
									 double epsilon, double pdtol=1e-3, int pdmaxiter=50) {

	vec xp, up, Atrp, AtAvp, fu1p, fu2p, fe1p, fe2p, lamu1p, lamu2p, lame1p, lame2p, rdp;
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
  double sdg = -sum(flu1) - sum(flu2) - sum(fle1) - sum(fle2);
	double tau = mu*4*N/sdg;

	vec gradf0 = join_cols(zeros(N), ones(N));
	vec rdual = gradf0 + join_cols(lamu1 - lamu2 + AtAv, -lamu1 - lamu2);
	vec rcent = -join_cols(join_cols(flu1, flu2), join_cols(fle1, fle2)) - 1/tau;
	double resnorm = sqrt(accu(join_cols(pow(rdual, 2), pow(rcent, 2))));

	int pditer = 0;
	bool done = (sdg < pdtol) || (pditer >= pdmaxiter);

	while (!done) {
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
			xp = x;
		  return List::create(Named("xp") = xp);
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

		vec smax1 = join_cols(-lamu1.elem(iu1)/dlamu1.elem(iu1),
													-lamu2.elem(iu2)/dlamu2.elem(iu2));
		vec smax2 = join_cols(-lame1.elem(ie1)/dlame1.elem(ie1),
													-lame2.elem(ie2)/dlame2.elem(ie2));
		vec smax3 = join_cols(-fu1.elem(ifu1)/(dx.elem(ifu1) - du.elem(ifu1)),
													-fu2.elem(ifu2)/(-dx.elem(ifu2) -du.elem(ifu2)));
		vec smax4 = join_cols(-fe1.elem(ife1)/AtAdx.elem(ife1),
													-fe2.elem(ife2)/(-AtAdx.elem(ife2)));
		vec smaxV = join_cols(join_cols(smax1, smax2), join_cols(smax3, smax4));
		double smax = min(smaxV);
		if (smax > 1) smax = 1;

    double s = 0.99*smax;

		bool suffdec = FALSE;
		int backiter = 0;

		while(!suffdec) {
			xp = x + s*dx;
			up = u + s*du;
			Atrp = Atr + s*AtAdx;
			AtAvp = AtAv + s*AtAdv;

			fu1p = fu1 + s*(dx - du);
			fu2p = fu2 + s*(-dx-du);
			fe1p = fe1 + s*AtAdx;
			fe2p = fe2 + s*(-AtAdx);

			lamu1p = lamu1 + s*dlamu1;
			lamu2p = lamu2 + s*dlamu2;
			lame1p = lame1 + s*dlame1;
			lame2p = lame2 + s*dlame2;

			rdp = gradf0 + join_cols(lamu1p - lamu2p + AtAvp, -lamu1p - lamu2p);
			vec rcp = -join_cols(join_cols(lamu1p%fu1p, lamu2p%fu2p),
													 join_cols(lame1p%fe1p, lame2p%fe2p)) - 1/tau;
			suffdec = sqrt(accu(join_cols(pow(rdp, 2), pow(rcp, 2)))) <
				(1 - alpha*s)*resnorm;

			s = beta*s;
			backiter = backiter + 1;

			if (backiter > 32) {
				Rcpp::warning("Backtracking stuck.  Previous iterate matrix returned!");
				xp = x;
			  return List::create(Named("xp") = xp);
			}
		}

		x = xp;
    u = up;
    Atr = Atrp;
    AtAv = AtAvp;
    fu1 = fu1p;
    fu2 = fu2p;
    fe1 = fe1p;
    fe2 = fe2p;
    lamu1 = lamu1p;
    lamu2 = lamu2p;
    lame1 = lame1p;
    lame2 = lame2p;

		flu1 = fu1 % lamu1;
		flu2 = fu2 % lamu2;
		fle1 = fe1 % lame1;
		fle2 = fe2 % lame2;
		sdg = -sum(flu1) - sum(flu2) - sum(fle1) - sum(fle2);
		tau = mu*4*N/sdg;
		rdual = rdp;
		rcent = -join_cols(join_cols(flu1, flu2), join_cols(fle1, fle2)) - 1/tau;
		resnorm = sqrt(accu(join_cols(pow(rdual, 2), pow(rcent, 2))));

		pditer = pditer + 1;
    done = (sdg < pdtol) || (pditer >= pdmaxiter);

	}

	return List::create(Named("xp") = xp);
}


/*** R
# library(microbenchmark)
set.seed(2)
X = matrix(rnorm(300), ncol=3)
emat = diag(3)
lam = 1
Sigma = var(X)
Omega0 = solve(Sigma)
linprogPD2C(Omega0[,1], Sigma, emat[,1], lam)
*/

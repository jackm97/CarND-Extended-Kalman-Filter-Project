#include "tools.h"
#include <iostream>
#include <cmath>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::sqrt;
using std::pow;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  int n = estimations.size();
  for (int i=0; i<n; i++){
    rmse += (estimations[i].array() - ground_truth[i].array()).pow(2).matrix();
  }
  rmse = rmse/n;
  rmse = rmse.array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  double px=x_state(0),
         py=x_state(1),
         vx=x_state(2),
         vy=x_state(3);
  
  double root_pxy = sqrt(px*px + py*py),
         pxy = root_pxy*root_pxy,
         root3_pxy = pxy*root_pxy;

  MatrixXd Hj(3,4);

  Hj << px/root_pxy, py/root_pxy, 0, 0,
        -py/pxy, px/pxy, 0, 0,
        py*(vx*py-vy*px)/root3_pxy, px*(vy*px-vx*py)/root3_pxy, px/root_pxy, py/root_pxy;

  return Hj;
}

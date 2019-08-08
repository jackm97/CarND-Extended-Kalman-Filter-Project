#include "kalman_filter.h"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::atan2;
using std::sqrt;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_*x_;
  P_ = F_*P_*F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  MatrixXd y = z - H_*x_;
  MatrixXd S = H_*P_*H_.transpose() + R_;
  MatrixXd K = P_*H_.transpose()*S.inverse();
  x_ = x_ + K*y;
  P_ = (MatrixXd::Identity(4,4)-K*H_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  VectorXd h(3);
  h << sqrt(x_(0)*x_(0) + x_(1)*x_(1)),
       atan2(x_(1),x_(0)),
       (x_(0)*x_(2) + x_(1)*x_(3))/sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  MatrixXd y = z - h;
  MatrixXd S = H_*P_*H_.transpose() + R_;
  MatrixXd K = P_*H_.transpose()*S.inverse();
  x_ = x_ + K*y;
  P_ = (MatrixXd::Identity(4,4)-K*H_)*P_;
}

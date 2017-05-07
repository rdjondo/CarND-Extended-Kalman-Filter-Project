#include "kalman_filter.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  /**
    * predict the state
  */
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */
	// From predicted state to predicted measurememt
	VectorXd z_pred = H_ * x_;

	/*
	 *Innovation or measurement residual:
	 *Innovation Residual between predicted measurement and actual measurement
	 */
	VectorXd y = z - z_pred;

	MatrixXd Ht = H_.transpose();
	/* Innovation (or residual) covariance 	 */
	MatrixXd S = H_ * P_ * Ht + R_;

	/* Preparing Kalman gain calculation */
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	/* Kalman gain */
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	unsigned int x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   *
    * update the state by using Extended Kalman Filter equations
   */

	static Tools tool;

	/*
	 *Innovation or measurement residual:
	 *Innovation Residual between predicted measurement and actual measurement
	 */
	VectorXd y = z - tool.h_radar(x_);

	/* Normalize angle phi between -PI and PI */
	float phi = y(1);
	while (phi < -M_PI)  {
		phi += 2.0 * M_PI;
	}
	phi = fmod(phi, M_PI);
	y(1) = phi;

	/* Linearise process measurement matrix */
	MatrixXd DH = tool.CalculateJacobian(x_);

	MatrixXd DHt = DH.transpose();

	/* Innovation (or residual) covariance 	 */
	MatrixXd S = DH * P_ * DHt + R_;

	/* Preparing Kalman gain calculation */
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * DHt;
	/* Kalman gain */
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	unsigned int x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * DH) * P_;
}

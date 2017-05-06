#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::endl;
using std::cout;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if(estimations.size() != ground_truth.size()){
		cout << "Estimation and ground_truth data sizes must be the same" << endl;
		return rmse;
	} else if (estimations.size() == 0){
		cout << "Error: estimation size is 0" << endl;
		return rmse;
	}

	// sum squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){
		VectorXd residual = estimations[i] - ground_truth[i];
		residual = residual.array()*residual.array(); // Why call to array() member function
		rmse += residual;
	}

	//calculate the average the squared residuals
	rmse = rmse/estimations.size();

	//calculate the squared root of the average
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */

	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//Precompute matrix elements
    float p2 = px*px+py*py;
    float sq_p2 = sqrt(p2);
    float sq_p2_3 = sqrt(p2*p2*p2);
    float cv = vx*py - vy*px;


	//check division by zero
	if(fabs(p2)<1e-4){
	    // this will make compute all matrix terms to zero
	    cout<<"Div by zero"<<endl;
	    p2 = 1;
	    px = 0;
	    py = 0;
	}

	//compute the Jacobian matrix
    Hj << px/sq_p2,   py/sq_p2,        0,        0,
         -py/p2,      px/p2,           0,        0,
          py*cv/sq_p2_3 , -px*cv/sq_p2_3, px/sq_p2, py/sq_p2;


	return Hj;
}

/* Transforms cartesian state-space coordinates to radar measurement polar coordinates  */
VectorXd Tools::h_radar(VectorXd& x_state){
	VectorXd z_pred =  VectorXd(3);

	//recover state parameters
	const float &px = x_state(0);
	const float &py = x_state(1);
	const float &vx = x_state(2);
	const float &vy = x_state(3);

	// Distance
	float d = sqrt(px*px+py*py);

	// Handle case where distance of detected object is zero to avoid divisions by zero
	float my_eps = 1e-4;

	z_pred(0) = d;

	if(fabsf(px)<my_eps){
		z_pred(1) = 0;
	} else{
		z_pred(1) = atan2f(py,px);
	}

	// d is already positive
	if(d<my_eps){
		z_pred.setZero();
	}
	else{
		z_pred(2) = (px*vx + py*vy)/d;
	}
	return z_pred;
}

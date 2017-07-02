#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
	//Jacobian matrix declaration.
		MatrixXd Hj(3,4);

	//recover state parameters
		float px = x_state(0);
		float py = x_state(1);
		float vx = x_state(2);
		float vy = x_state(3);

	//check division by zero
	    float d  = px*px+py*py;
	    float d_12 = sqrt(d);
	    float d_32 = (d*d_12);

	    //cout << "d,d_12,d_32: " << d << "," << d_12 << "," << d_32 << endl;
		if(fabs(d) < 0.0001){
		    cout << "DIVIDE BY ZERO EXCEPTION" << "\n";
		    return Hj;
		}

		//compute the Jacobian matrix

	    Hj << px/d_12,                 py/d_12,                   0,         0,
	          -py/d  ,                 px/d,                      0,         0,
	          (py*(vx*py-vy*px))/d_32, (px*(px*vy-py*vx))/d_32,   px/d_12,   py/d_12;

	    return Hj;

}

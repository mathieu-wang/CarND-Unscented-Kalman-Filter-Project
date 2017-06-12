#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  // ... your code here
  if (estimations.size() == 0) {
    std::cout << "Empty estimations vector." << std::endl;
    return rmse;
  }

  if (estimations.size() != ground_truth.size()) {
    std::cout << "Estimation size not equal to ground truth size." << std::endl;
    return rmse;
  }

  //accumulate squared residuals
  for (unsigned int i = 0; i < estimations.size(); ++i) {
    // ... your code here
    VectorXd diff = estimations[i] - ground_truth[i];
    rmse += VectorXd(diff.array()*diff.array());
  }

  //calculate the mean
  // ... your code here
  rmse = rmse / estimations.size();

  //calculate the squared root
  // ... your code here
  rmse = VectorXd(rmse.array().sqrt());

  //return the result
  return rmse;
}

VectorXd Tools::CartesianToPolar(const VectorXd& x) {
  if (x.size() != 4) {
    throw std::invalid_argument("Cartesian vector length is not 4");
  }

  double px = x[0];
  double py = x[1];
  double vx = x[2];
  double vy = x[3];

  if (fabs(px) < 0.0001) {
    std::cout << "CartesianToPolar() - px == 0 - Set to 0.0001" << std::endl;
    px = 0.0001;
  }

  double c1 = px*px + py*py;
  if (fabs(c1) < 0.0001) {
    std::cout << "CartesianToPolar() - Divide by Zero - Set to 0.0001" << std::endl;
    c1 = 0.0001;
  }

  double c2 = sqrt(c1);
  double c3 = px*vx + py*vy;

  double rho = c2;
  double phi = atan2(py, px);
  double rho_dot = c3 / c2;

  VectorXd x_polar(3);
  x_polar << rho, phi, rho_dot;
  return x_polar;
}

VectorXd Tools::PolarToCartesian(const VectorXd& x) {
  if (x.size() != 3) {
    throw std::invalid_argument("Polar vector length is not 3");
  }

  double rho = x[0];
  double phi = x[1];
  double rho_dot = x[2];

  double px = rho * cos(phi);
  double py = rho * sin(phi);
  double vx = rho_dot * cos(phi);
  double vy = rho_dot * sin(phi);

  VectorXd x_cartesian(5);
  x_cartesian << px,py,sqrt(vx*vx +vy*vy),0,0;
  return x_cartesian;
}

double Tools::normalizeAngle(double angle_rad) {
  while (angle_rad < -M_PI) {
    angle_rad += 2 * M_PI;
  }
  while (angle_rad > M_PI) {
    angle_rad -= 2 * M_PI;
  }
  return angle_rad;
}
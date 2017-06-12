#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * A helper method to convert a Cartesian vector to a polar vector.
  *
  */
  VectorXd CartesianToPolar(const VectorXd& x);

  /**
  * A helper method to convert a polar vector to a Cartesian vector.
  *
  */
  VectorXd PolarToCartesian(const VectorXd& x);

  /**
  * A helper method to normalize the angle in radians to -Pi and Pi
  */
  double normalizeAngle(double angle_rad);

};

#endif /* TOOLS_H_ */
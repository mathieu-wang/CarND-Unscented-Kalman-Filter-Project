#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd::Zero(5);

  // initial covariance matrix
  P_ = MatrixXd::Zero(5, 5);

  H_laser_ = MatrixXd::Zero(2, 5);

  R_laser_ = MatrixXd::Zero(2, 2);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;//30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2;//30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  n_x_ = 5;

  // Predicted sigma points as columns
  n_aug_ = 7;

  //create augmented mean vector
  x_aug_ = VectorXd::Zero(n_aug_);

  //create augmented state covariance
  P_aug_ = MatrixXd::Zero(n_aug_, n_aug_);

  lambda_ = 3 - n_aug_;

  Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);
  Xsig_aug_ = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);

  //set vector for weights
  weights_ = VectorXd::Zero(2*n_aug_+1);


  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    x_ << 1, 1, 1, 1, 1;

    H_laser_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0;

    //measurement covariance matrix - laser
    R_laser_ << std_laspx_*std_laspx_, 0,
                0, std_laspy_*std_laspy_;

    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      x_ << tools.PolarToCartesian(meas_package.raw_measurements_);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      //set the state with the initial location and zero velocity
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    time_us_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;

    std::cout << "Done initializing" << endl;
    return;
  }

  //compute the time elapsed between the current and previous measurements
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  time_us_ = meas_package.timestamp_;
  cout << "dt: " << dt << endl;
  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  cout << "Starting prediction" << endl;
  // Augmentation:
  //create augmented mean state
  x_aug_.fill(0.0);
  x_aug_.head(x_.size()) = x_;
  //create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(P_.rows(), P_.cols()) = P_;
  MatrixXd Q = MatrixXd::Zero(2, 2);
  Q(0, 0) = std_a_*std_a_;
  Q(1, 1) = std_yawdd_*std_yawdd_;
  P_aug_.bottomRightCorner(Q.rows(), Q.cols()) = Q;

  MatrixXd A_aug = P_aug_.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0) = x_aug_;
  double factor = sqrt(lambda_+n_aug_);

  for (int i = 1; i < n_aug_ + 1; i++) {
    VectorXd rt =  factor * A_aug.col(i-1);
    Xsig_aug_.col(i) = x_aug_ + rt;
    Xsig_aug_.col(i + n_aug_) = x_aug_ - rt;
  }

  //Predict sigma points:
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  //set weights
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i = 1; i < 2*n_aug_+1; i++) {
    weights_(i) = 1 / (2*(lambda_ + n_aug_));
  }
  //predict state mean
  x_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }
  //predict state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = tools.normalizeAngle(x_diff(3));
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }

  cout << "Done prediction. x: " << endl << x_ << endl
       << "P: " << endl << P_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  cout << "Start UpdateLidar" << endl;
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_pred = H_laser_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;
  cout << "End UpdateLidar" << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  cout << "Start UpdateRadar" << endl;
  int n_z = 3;
  MatrixXd Zsig = MatrixXd::Zero(n_z, 2 * n_aug_ + 1);
  //mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);
  //create matrix for predicted measurement covariance
  MatrixXd S = MatrixXd::Zero(n_z, n_z);
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
  VectorXd z = meas_package.raw_measurements_;

  //transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_+1; i++) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double psi = Xsig_pred_(3, i);
    double psi_dot = Xsig_pred_(4, i);
    if (px == 0 && py == 0) {
      Zsig(0,i) = 0;
      Zsig(1,i) = 0;
      Zsig(2,i) = 0;
    } else {
      Zsig(0, i) = sqrt(px * px + py * py);
      Zsig(1, i) = atan2(py, px);
      Zsig(2, i) = (px * cos(psi) * v + py * sin(psi) * v) / Zsig(0, i);
    }
  }
  //calculate mean predicted measurement
  for (int i = 0; i < 2*n_aug_+1; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }
  //calculate measurement covariance matrix S
  for (int i = 0; i < 2*n_aug_+1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = tools.normalizeAngle(z_diff(1));
    S += weights_(i) * z_diff * z_diff.transpose();
  }
  S(0, 0) += std_radr_*std_radr_;
  S(1, 1) += std_radphi_*std_radphi_;
  S(2, 2) += std_radrd_*std_radrd_;

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1) = tools.normalizeAngle(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - Xsig_pred_.col(0);//x_;
    //angle normalization
    x_diff(3) = tools.normalizeAngle(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  z_diff(1) = tools.normalizeAngle(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  cout << "End UpdateRadar" << endl;
  cout << "x: " << x_ << endl;
  cout << "P: " << P_ << endl;
}

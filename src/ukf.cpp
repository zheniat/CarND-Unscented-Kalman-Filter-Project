#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;
  
  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  
  // State dimension
  n_x_ = 5;
  
  ///* Augmented state dimension
  n_aug_ = 7;
  
  // initial state vector
  x_ = VectorXd::Zero(n_x_);
  
  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ << .1, 0, 0, 0, 0,
  0, .1, 0, 0, 0,
  0, 0, .1, 0, 0,
  0, 0, 0, .1, 0,
  0, 0, 0, 0, .1;
  
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;
  
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  ///* initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;
  
  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);
  
  ///* time when the state is true, in us
  time_us_ = 0;
  
  
  
  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;
  
  ///* Weights of sigma points
  weights_ = VectorXd::Zero(2*n_aug_+1);
  
  //set weights
  double w = lambda_ /(lambda_ + n_aug_);
  weights_(0) = w;
  for(int i=1;i<2*n_aug_ + 1; i++){
    w = 0.5 / (lambda_ + n_aug_);;
    weights_(i) = w;
  }
  
  NIS_radar_ = 0;
  NIS_lidar_ = 0;
}

UKF::~UKF() {
}


double UKF::NormalizeAngle(double angle){
  //angle normalization
  while (angle> M_PI) angle-=2.*M_PI;
  while (angle<-M_PI) angle+=2.*M_PI;
  return angle;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    double px = 0;
    double py = 0;
    double vx = 0;
    double vy = 0;
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
       Convert radar from polar to cartesian coordinates and initialize state.
       */
      double rho = meas_package.raw_measurements_[0];
      double theta = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];
      px =  rho * cos(theta);
      py = rho * sin(theta);
      vx = rho_dot * cos(theta);
      vy = rho_dot * sin(theta);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
       Initialize state.
       */
      px = meas_package.raw_measurements_[0];
      py = meas_package.raw_measurements_[1];
      vx = 0;
      vy = 0;
    }
    
    // Fix small px, py
    if(fabs(px) < 0.0001){
      px = 0.01;
    }
    if(fabs(py) < 0.0001){
      py = 0.01;
    }
    
    x_ << px, py, sqrt(pow(vx, 2) + pow(vy, 2)), 0, 0;
    
    time_us_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  
  //compute the time elapsed between the current and previous measurements
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  time_us_ = meas_package.timestamp_;
  
  Prediction(dt);
  /*****************************************************************************
   *  Update
   ****************************************************************************/
  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
    cout << "NIS_radar_ = " << NIS_radar_  << endl;
    std::ofstream f;
    f.open("radar.csv",std::ios_base::app);
    f<<NIS_radar_<<",";
    f.close();
  } else if(use_laser_  && meas_package.sensor_type_ == MeasurementPackage::LASER){
    UpdateLidar(meas_package);
    cout << "NIS_lidar_ = " << NIS_lidar_  << endl;
    std::ofstream f;
    f.open("lidar.csv",std::ios_base::app);
    f<<NIS_lidar_<<",";
    f.close();

  }
  
  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  //Generate augmented sigma points
  //create augmented mean vector
  VectorXd x_aug = VectorXd::Zero(n_aug_);
  
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  
  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);
  
  //create augmented mean state
  x_aug.head(5) = x_;
  //create augmented covariance matrix
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  
  MatrixXd Q = MatrixXd(2, 2);
  Q << pow(std_a_,2), 0,
  0, pow(std_yawdd_,2);
  
  P_aug.bottomRightCorner(2,2) = Q;
  
  //create square root matrix
  MatrixXd P_aug_sqrt = P_aug.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * P_aug_sqrt.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * P_aug_sqrt.col(i);
  }
  
  //Predict sigma points
  
  for (int i = 0; i< 2 * n_aug_ + 1; i++)
  {
    VectorXd x = Xsig_aug.col(i);
    double px = x(0);
    double py = x(1);
    double vk = x(2);
    double psi = x(3);
    double psi_dot = x(4);
    double va = x(5);
    double v_psi_ddot = x(6);
    
    //Output x
    VectorXd x_out = VectorXd(n_x_);
    x_out << px, py, vk, psi, psi_dot;
    
    //process noise
    VectorXd noise = VectorXd(n_x_);
    double n1 = 0.5 * pow(delta_t,2) * cos(psi) * va;
    double n2 = 0.5 * pow(delta_t,2) * sin(psi) * va;
    double n3 = delta_t * va;
    double n4 = 0.5 * pow(delta_t,2) * v_psi_ddot;
    double n5 = delta_t * v_psi_ddot;
    noise << n1, n2, n3, n4, n5;
    
    VectorXd state = VectorXd(n_x_);
    double s1,s2,s3,s4,s5;
    if(psi_dot == 0){
      s1 = vk * cos(psi) * delta_t;
      s2 = vk * sin(psi) * delta_t;
    } else {
      s1 = vk/psi_dot * (sin(psi+psi_dot * delta_t)-sin(psi));
      s2 = vk/psi_dot * (-cos(psi + psi_dot * delta_t)+cos(psi));
    }
    s3 = 0;
    s4 = psi_dot * delta_t;
    s5 = 0;
    state << s1, s2, s3, s4, s5;
    
    Xsig_pred_.col(i) = x_out + state + noise;
  }
  
  
  //Predict state and covariance
  //vector for predicted state
  VectorXd x = VectorXd::Zero(n_x_);
  
  //covariance matrix for prediction
  MatrixXd P = MatrixXd::Zero(n_x_, n_x_);
  
  //predict state mean
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    x = x + weights_(i) * Xsig_pred_.col(i);
  }
  
  //predict state covariance matrix
  for(int i=0;i<2*n_aug_+1;i++){
    VectorXd diff = Xsig_pred_.col(i) - Xsig_pred_.col(0);
    diff(3) = NormalizeAngle(diff(3));
    P = P + weights_(i) * diff * diff.transpose();
  }
  x_ = x;
  P_ = P;
}

void UKF::UpdateState(VectorXd &z, int n_z, MatrixXd &Zsig, MatrixXd &S, VectorXd &z_pred, bool isRadar){
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
  
  //calculate cross correlation matrix
  for(int i=0; i<2*n_aug_+1; i++){
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    VectorXd z_diff = Zsig.col(i) - z_pred;
    if(isRadar){
      z_diff(1) = NormalizeAngle(z_diff(1));
    }
    x_diff(3) = NormalizeAngle(x_diff(3));
    Tc= Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //residual
  VectorXd z_diff = z - z_pred;
  
  if(isRadar){
    z_diff(1) = NormalizeAngle(z_diff(1));
  }
  
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
  
  //NIS
  double nis = z_diff.transpose() * S.inverse() * z_diff;
  if(isRadar){
    NIS_radar_ = nis;
  } else {
    NIS_lidar_ = nis;
  }
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  //predict measurements
  int n_z = 3;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z, 2 * n_aug_ + 1);
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z,n_z);
  
  //transform sigma points into measurement space
  for(int i=0; i<2 * n_aug_ + 1; i++){
    VectorXd x = Xsig_pred_.col(i);
    double px = x(0);
    double py = x(1);
    double v = x(2);
    double psi = x(3);
    
    double rho = sqrt(pow(px,2) + pow(py,2));
    double phi = 0;
    if(px!=0){
      phi = atan2(py,px);
    }
    
    double rho_dot = 0;
    if(rho != 0){
      rho_dot = (px * cos(psi) * v + py * sin(psi) * v)/rho;
    }
    
    Zsig.col(i) << rho, phi, rho_dot;
  }
  //calculate mean predicted measurement
  z_pred = Zsig * weights_;
  //calculate innovation covariance matrix S
  MatrixXd R = MatrixXd::Zero(n_z, n_z);
  R << pow(std_radr_,2), 0, 0,
  0, pow(std_radphi_,2), 0,
  0, 0, pow(std_radrd_, 2);
  
  //predict state covariance matrix
  for(int i=0;i<2*n_aug_+1; i++){
    VectorXd diff = Zsig.col(i) - z_pred;
    diff(1) = NormalizeAngle(diff(1));
    S = S + weights_(i) * diff * diff.transpose();
  }
  S = S + R;
  //update state
  //vector for incoming radar measurement
  VectorXd z = VectorXd::Zero(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];
  UpdateState(z, n_z, Zsig, S, z_pred,true);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  //predict measurements
  int n_z = 2;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z, 2 * n_aug_ + 1);
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z,n_z);
  
  //transform sigma points into measurement space
  for(int i=0; i<2 * n_aug_ + 1; i++){
    VectorXd x = Xsig_pred_.col(i);
    double px = x(0);
    double py = x(1);
    
    if (fabs(px) < 0.0001 || fabs(py)< 0.0001) {
      px = 0.0001;
      py = 0.0001;
    }
    Zsig.col(i) << px, py;
  }
    
  //calculate mean predicted measurement
  z_pred = Zsig * weights_;

  //predict state covariance matrix
  for(int i=0;i<2*n_aug_+1; i++){
    VectorXd diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * diff * diff.transpose();
  }
  
  MatrixXd R = MatrixXd::Zero(n_z, n_z);
  R << pow(std_laspx_, 2), 0,
  0, pow(std_laspy_, 2);
  
  S = S + R;
  //update state
  //vector for incoming radar measurement
  VectorXd z = VectorXd::Zero(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
  
  //update state
  UpdateState(z, n_z, Zsig, S, z_pred, false);
}

#pragma once
#include <Eigen/Dense>
#include <cmath>

class NavEstimator_SRUKF {
protected:
    int n_states = 6;
    Eigen::VectorXf x;
    Eigen::MatrixXf S;
    Eigen::MatrixXf Q_sqrt;
    Eigen::MatrixXf R;

public:
    NavEstimator_SRUKF() = default;
    virtual ~NavEstimator_SRUKF() = default;

    virtual void initialize(int n) {
        n_states = n;
        x = Eigen::VectorXf::Zero(n);
        S = Eigen::MatrixXf::Identity(n, n) * 0.1f;
        Q_sqrt = Eigen::MatrixXf::Identity(n, n) * 0.01f;
        R = Eigen::MatrixXf::Identity(3, 3) * 0.1f;
    }

    virtual Eigen::VectorXf transition_model(
        const Eigen::VectorXf& s,
        float dt,
        const Eigen::Vector3f& accel) = 0;

    virtual Eigen::Vector3f observation_model(const Eigen::VectorXf& s) = 0;

    virtual void predict(float dt, const Eigen::Vector3f& accel) {
        // Sigma points generation
        int L = n_states;
        float lambda = 3.0f - L;
        Eigen::MatrixXf sigma_points = Eigen::MatrixXf::Zero(2*L + 1, L);
        
        sigma_points.row(0) = x.transpose();
        Eigen::MatrixXf chol = S.llt().matrixL();
        
        for (int i = 0; i < L; ++i) {
            sigma_points.row(i + 1) = (x + std::sqrt(L + lambda) * chol.row(i)).transpose();
            sigma_points.row(L + i + 1) = (x - std::sqrt(L + lambda) * chol.row(i)).transpose();
        }

        // Predict sigma points
        Eigen::MatrixXf pred_sigma = Eigen::MatrixXf::Zero(2*L + 1, L);
        for (int i = 0; i < 2*L + 1; ++i) {
            pred_sigma.row(i) = transition_model(sigma_points.row(i).transpose(), dt, Eigen::Vector3f::Zero()).transpose();
        }

        // Mean and covariance update
        float w0 = lambda / (L + lambda);
        float wi = 0.5f / (L + lambda);
        
        x = w0 * pred_sigma.row(0).transpose();
        for (int i = 1; i < 2*L + 1; ++i) {
            x += wi * pred_sigma.row(i).transpose();
        }

        S = Eigen::MatrixXf::Zero(L, L);
        for (int i = 0; i < 2*L + 1; ++i) {
            float weight = (i == 0) ? w0 : wi;
            S += weight * (pred_sigma.row(i).transpose() - x) * (pred_sigma.row(i) - x.transpose());
        }
        S += Q_sqrt;
    }

    virtual void update(const Eigen::Vector3f& z) {
        // Simple update (simplified for demo)
        Eigen::Vector3f z_pred = observation_model(x);
        Eigen::Vector3f residual = z - z_pred;
        
        float gain = 0.1f;
        x.head<3>() += gain * residual;
    }
};

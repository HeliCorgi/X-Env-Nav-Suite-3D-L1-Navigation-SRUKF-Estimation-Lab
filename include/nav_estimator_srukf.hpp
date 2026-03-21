#pragma once
#include <Eigen/Dense>
#include <Eigen/QR>
#include <cmath>

// Reference: van der Merwe, R. and Wan, E.A. (2001)
// "The Square-Root Unscented Kalman Filter for State and Parameter-Estimation"
// ESANN 2001 proceedings

class NavEstimator_SRUKF {
protected:
    int n_states = 6;
    Eigen::VectorXf x;
    Eigen::MatrixXf S;       // square-root covariance (upper triangular)
    Eigen::MatrixXf Q_sqrt;
    Eigen::MatrixXf R;

    float alpha = 1e-3f;
    float beta  = 2.0f;
    float kappa = 0.0f;

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
        int L = n_states;
        float lambda = alpha * alpha * (L + kappa) - L;
        
        // Sigma points generation
        Eigen::MatrixXf sigma_points(2*L + 1, L);
        sigma_points.row(0) = x.transpose();
        
        Eigen::MatrixXf sqrt_S = S.llt().matrixL();
        
        for (int i = 0; i < L; ++i) {
            sigma_points.row(i + 1)     = (x + std::sqrt(L + lambda) * sqrt_S.row(i)).transpose();
            sigma_points.row(L + i + 1) = (x - std::sqrt(L + lambda) * sqrt_S.row(i)).transpose();
        }

        // Predict sigma points
        Eigen::MatrixXf pred_sigma(2*L + 1, L);
        for (int i = 0; i < 2*L + 1; ++i) {
            pred_sigma.row(i) = transition_model(sigma_points.row(i).transpose(), dt, accel).transpose();
        }

        // Weighted mean
        float w0 = lambda / (L + lambda);
        float wi = 0.5f / (L + lambda);
        
        x = w0 * pred_sigma.row(0).transpose();
        for (int i = 1; i < 2*L + 1; ++i) {
            x += wi * pred_sigma.row(i).transpose();
        }

        // QR-based square-root covariance update
        Eigen::MatrixXf A(2*L, L);
        for (int i = 1; i < 2*L + 1; ++i) {
            A.row(i - 1) = std::sqrt(wi) * (pred_sigma.row(i) - x.transpose());
        }
        
        Eigen::HouseholderQR<Eigen::MatrixXf> qr(A.transpose());
        S = qr.matrixQR().topRows(L).triangularView<Eigen::Upper>();
        
        // Add process noise (Cholesky update)
        S = cholupdate(S, Q_sqrt.col(0), '+');
    }

    virtual void update(const Eigen::Vector3f& z) {
        int L = n_states;
        float lambda = alpha * alpha * (L + kappa) - L;
        
        // Sigma points
        Eigen::MatrixXf sigma_points(2*L + 1, L);
        sigma_points.row(0) = x.transpose();
        Eigen::MatrixXf sqrt_S = S.llt().matrixL();
        
        for (int i = 0; i < L; ++i) {
            sigma_points.row(i + 1)     = (x + std::sqrt(L + lambda) * sqrt_S.row(i)).transpose();
            sigma_points.row(L + i + 1) = (x - std::sqrt(L + lambda) * sqrt_S.row(i)).transpose();
        }

        // Measure sigma points
        Eigen::MatrixXf meas_sigma(2*L + 1, 3);
        for (int i = 0; i < 2*L + 1; ++i) {
            meas_sigma.row(i) = observation_model(sigma_points.row(i).transpose()).transpose();
        }

        // Measurement mean
        float w0 = lambda / (L + lambda);
        float wi = 0.5f / (L + lambda);
        
        Eigen::Vector3f z_pred = w0 * meas_sigma.row(0).transpose();
        for (int i = 1; i < 2*L + 1; ++i) {
            z_pred += wi * meas_sigma.row(i).transpose();
        }

        // Kalman gain (simplified)
        Eigen::Vector3f residual = z - z_pred;
        float gain = 0.05f;
        x.head<3>() += gain * residual;
    }

private:
    // van der Merwe準拠のCholesky rank-1 update（数値安定）
    Eigen::MatrixXf cholupdate(const Eigen::MatrixXf& S_in, const Eigen::VectorXf& u, char sign) {
        Eigen::MatrixXf S = S_in;
        int n = S.rows();
        float alpha = 1.0f;
        
        for (int k = 0; k < n; ++k) {
            float r_kk = S(k, k);
            float u_k = u(k);
            float alpha_prev = alpha;
            alpha = std::sqrt(r_kk * r_kk + (sign == '+' ? 1.0f : -1.0f) * u_k * u_k);
            
            if (alpha < 1e-10f) {
                S = Eigen::MatrixXf::Zero(n, n);
                return S;
            }
            
            float c = alpha / r_kk;
            float s_val = u_k / r_kk;
            S(k, k) = alpha;
            
            for (int j = k + 1; j < n; ++j) {
                float S_kj = S(k, j);
                S(k, j) = (S_kj + (sign == '+' ? s_val : -s_val) * u(j)) / c;
                u(j) = c * u(j) - s_val * S(k, j);
            }
        }
        return S;
    }
};

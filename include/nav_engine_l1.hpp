#pragma once
#include <array>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include "env_adaptive_core.hpp"

// Reference: Titterton, D. and Weston, J. (2004)
// "Strapdown Inertial Navigation Technology" (2nd Edition), AIAA
// Groves, P.D. (2013) "Principles of GNSS, Inertial, and Integrated Navigation Systems"

struct State_L1 {
    std::array<float, 13> x{}; 
    float p_limit = 15000.0f;
    float r_mass  = 45.0f;
};

class NavEngine_L1 {
private:
    State_L1 state;
    static constexpr double MASS = 45.0;
    static constexpr double AREA = 1.0;
    static constexpr double OMEGA = 7.292115e-5;   // rad/s
    double latitude = 45.0 * M_PI / 180.0;         // 初期緯度（可変可）

    static double get_cd(double mach) {
        if (mach < 0.8)      return 0.45;
        else if (mach < 1.2) return 0.45 + 3.75 * (mach - 0.8);
        else if (mach < 5.0) return 1.2 - 0.18 * (mach - 1.2);
        else                 return 0.3;
    }

public:
    NavEngine_L1() {
        state.x.fill(0.0f);
        state.x[12] = 1.0f;
    }

    void update(float meas_norm, float dt) {
        float r = 45.0f;
        state.x[6] += 0.015f * (meas_norm - r) * dt;
    }

    void predict(float dt, const Eigen::Vector3f& obs_accel,
                 double alt_m, double vel_mps, const EnvAdaptive_Core& env) {
        double rho = env.get_air_density(alt_m);
        double g   = env.get_gravity(alt_m);
        double a   = env.get_speed_of_sound(alt_m);
        double mach = (vel_mps > 0.0) ? vel_mps / a : 0.0;
        double cd  = get_cd(mach);

        Eigen::Vector3f vel{state.x[3], state.x[4], state.x[5]};
        double v2 = vel.squaredNorm();

        Eigen::Vector3f accel = obs_accel;

        // マッハ依存抗力
        if (v2 > 1e-6) {
            accel += -0.5 * rho * cd * AREA / MASS * v2 * vel.normalized();
        }

        // 重力（高度依存）
        accel[2] -= g;

        // 緯度依存コリオリ + 地球曲率補正（NED座標）
        double sin_lat = std::sin(latitude);
        double cos_lat = std::cos(latitude);
        double f = 2.0 * OMEGA * sin_lat;                    // Coriolis parameter
        double R = 6378137.0 + alt_m;                        // 地球曲率半径
        double tan_lat_over_R = std::tan(latitude) / R;      // 曲率項

        accel[0] += f * vel[1] + tan_lat_over_R * vel[0] * vel[1];   // East
        accel[1] -= f * vel[0] + tan_lat_over_R * vel[0] * vel[0];   // North

        // 位置・速度更新（球面近似）
        for (int i = 0; i < 3; ++i) {
            state.x[3 + i] += accel[i] * dt;
            state.x[i]     += state.x[3 + i] * dt + 0.5f * accel[i] * dt * dt;
        }

        // スケール更新
        float res = obs_accel.norm() - accel.norm();
        state.x[12] = std::clamp(state.x[12] + res * dt * 0.5f, 0.7f, 1.3f);
    }

    const State_L1& get_state() const { return state; }
    Eigen::Vector3f get_pos() const { return {state.x[0], state.x[1], state.x[2]}; }
    Eigen::Vector3f get_vel() const { return {state.x[3], state.x[4], state.x[5]}; }
};

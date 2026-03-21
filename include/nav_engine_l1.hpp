#pragma once
#include <array>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include "env_adaptive_core.hpp"

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
    static constexpr double OMEGA = 7.292115e-5f;   // 地球自転角速度

    static double get_cd(double mach) {
        if (mach < 0.8)      return 0.45f;
        else if (mach < 1.2) return 0.45f + 3.75f * (mach - 0.8f);
        else if (mach < 5.0) return 1.2f - 0.18f * (mach - 1.2f);
        else                 return 0.3f;
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

        // 空力抗力（マッハ依存）
        if (v2 > 1e-6) {
            Eigen::Vector3f drag = -0.5 * rho * cd * AREA / MASS * v2 * vel.normalized();
            accel += drag;
        }

        // 重力（高度依存）
        accel[2] -= g;

        // コリオリ力（北半球簡易モデル）
        constexpr float f = 2.0f * OMEGA * 0.707f;
        accel[0] += f * vel[1];
        accel[1] -= f * vel[0];

        // 状態更新
        for (int i = 0; i < 3; ++i) {
            state.x[3 + i] += accel[i] * dt;
            state.x[i]     += state.x[3 + i] * dt + 0.5f * accel[i] * dt * dt;
        }

        float res = obs_accel.norm() - accel.norm();
        state.x[12] = std::clamp(state.x[12] + res * dt * 0.5f, 0.7f, 1.3f);
    }

    const State_L1& get_state() const { return state; }
    Eigen::Vector3f get_pos() const { return {state.x[0], state.x[1], state.x[2]}; }
    Eigen::Vector3f get_vel() const { return {state.x[3], state.x[4], state.x[5]}; }
};

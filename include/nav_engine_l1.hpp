#pragma once
#include <array>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

struct State_L1 {
    std::array<float, 13> x{}; 
    float p_limit = 15000.0f;
    float r_mass  = 45.0f;
};

class NavEngine_L1 {
private:
    State_L1 state;

    std::array<float, 3> get_ref() const {
        return {45.0f, 0.0f, 0.0f};
    }

    Eigen::Vector3f get_acc(float s) const {
        return Eigen::Vector3f(4.905f * s, 4.905f * s * 0.2f, 4.905f * s * 0.1f);
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

    void predict(float dt, const Eigen::Vector3f& obs_accel) {
        float res = obs_accel.norm() - get_acc(state.x[12]).norm();
        state.x[12] = std::clamp(state.x[12] + res * dt * 0.5f, 0.7f, 1.3f);

        for (int i = 0; i < 3; ++i) {
            state.x[3 + i] += obs_accel[i] * dt;
            state.x[i] += state.x[3 + i] * dt + 0.5f * obs_accel[i] * dt * dt;
        }
    }

    const State_L1& get_state() const { return state; }
    Eigen::Vector3f get_pos() const {
        return Eigen::Vector3f(state.x[0], state.x[1], state.x[2]);
    }
};

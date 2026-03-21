#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include "nav_engine_l1.hpp"
#include "nav_estimator_srukf.hpp"
#include "env_adaptive_core.hpp"

class L1_SRUKF : public NavEstimator_SRUKF {
public:
    L1_SRUKF() {
        initialize(6);  // pos3 + vel3
        x = Eigen::VectorXf::Zero(6);
        S = Eigen::MatrixXf::Identity(6,6)*0.05f;
        Q_sqrt = Eigen::MatrixXf::Identity(6,6)*0.01f;
    }

    Eigen::VectorXf transition_model(
        const Eigen::VectorXf& s,
        float dt,
        const Eigen::Vector3f& accel) override
    {
        Eigen::VectorXf new_s = s;
        new_s.head<3>() += s.segment<3>(3)*dt + 0.5f*accel*dt*dt;
        new_s.segment<3>(3) += accel*dt;
        return new_s;
    }

    Eigen::Vector3f observation_model(const Eigen::VectorXf& s) override {
        return s.head<3>();
    }
};

int main() {
    NavEngine_L1 engine;
    L1_SRUKF ukf;
    EnvAdaptive_Core env;

    const float dt = 0.02f;

    std::cout << std::fixed << std::setprecision(3);

    for (int i = 0; i < 500; ++i) {
        double alt  = engine.get_pos().z();
        double speed = engine.get_vel().norm();
        Eigen::Vector3f meas_accel(0.20f, 0.10f, 0.05f);

        env.update_thermal(dt, alt, speed);
        engine.predict(dt, meas_accel, alt, speed, env);

        // SRUKF更新（位置＋速度観測）
        ukf.predict(dt, meas_accel);
        Eigen::VectorXf z(6); z << engine.get_pos(), engine.get_vel();
        Eigen::VectorXf noise = Eigen::VectorXf::Random(6) * 0.02f;
        ukf.update(z + noise);

        double mach = (speed > 0.0) ? speed / env.get_speed_of_sound(alt) : 0.0;
        const auto& raw = engine.get_state();

        if (i % 50 == 0) {  // 50ステップごとに表示
            std::cout << "[t=" << (i*dt) << "s] Mach=" << mach 
                      << "  Alt=" << alt << "m  rho=" << env.get_air_density(alt) << "\n";
        }
    }
    return 0;
}

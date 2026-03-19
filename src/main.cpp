#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include "nav_engine_l1.hpp"
#include "nav_estimator_srukf.hpp"
#include "env_adaptive_core.hpp"

// ================================
// SRUKFラッパ（L1航法専用 6次元: pos3+vel3）
// ================================
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
    Eigen::Vector3f meas_accel(0.20f, 0.10f, 0.05f);

    std::cout << std::fixed << std::setprecision(3);

    for (int i=0;i<250;++i){
        engine.update(meas_accel.norm(), dt);
        engine.predict(dt, meas_accel);

        ukf.predict(dt, meas_accel);
        Eigen::Vector3f noisy_z = engine.get_pos() + Eigen::Vector3f::Random()*0.02f;
        ukf.update(noisy_z);

        env.update_thermal(dt, 1200.0);

        const auto& raw = engine.get_state();
        std::cout << "[t=" << (i*dt) << "s] "
                  << "Raw: pos=(" << raw.x[0] << "," << raw.x[1] << "," << raw.x[2] << ") "
                  << "vel=(" << raw.x[3] << "," << raw.x[4] << "," << raw.x[5] << ") "
                  << "scale=" << raw.x[12] << "\n"
                  << "     Est: pos=(" << ukf.x[0] << "," << ukf.x[1] << "," << ukf.x[2] << ") "
                  << "vel=(" << ukf.x[3] << "," << ukf.x[4] << "," << ukf.x[5] << ")\n"
                  << "     PCM=" << (env.get_pcm_remaining()*100) << "% "
                  << "Temp=" << env.get_temperature() << "℃\n\n";
    }
    return 0;
}

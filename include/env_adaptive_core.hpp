#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <array>

class EnvAdaptive_Core {
private:
    struct Health {
        double temp = 20.0;
        double pcm_ratio = 1.0;
        Eigen::Matrix3d strain = Eigen::Matrix3d::Zero();
    } h;

    static constexpr double PCM_MELT_ENERGY = 2.6e6;      // J/kg
    static constexpr double R_EARTH = 6378137.0;          // WGS84 m
    static constexpr double G0 = 9.80665;
    static constexpr double R_AIR = 287.058;
    static constexpr double K_SUTTON = 1.7415e-4;         // NASA TR R-376
    static constexpr double RN = 1.0;                     // 鼻半径 m（簡易）

    // US76層テーブル（U.S. Standard Atmosphere, 1976）
    static constexpr std::array<double, 7> h_base{0,11000,20000,32000,47000,51000,71000};
    static constexpr std::array<double, 7> lapse{-0.0065,0.0,0.001,0.0028,0.0,-0.0028,-0.002};
    static constexpr std::array<double, 7> T_base{288.15,216.65,216.65,228.65,270.65,270.65,214.65};
    static constexpr std::array<double, 7> P_base{101325.0,22632.1,5474.89,868.019,110.906,66.939,3.956};

public:
    void get_atmosphere(double alt_m, double& rho, double& T, double& P) const {
        // (前回と同じUS76実装を省略せず完全収録)
        if (alt_m < 0.0) alt_m = 0.0;
        int layer = 0;
        for (int i = 0; i < 6; ++i) {
            if (alt_m < h_base[i+1]) { layer = i; break; }
        }
        double h_b = h_base[layer], L = lapse[layer], Tb = T_base[layer], Pb = P_base[layer];
        double h = alt_m;
        if (std::abs(L) < 1e-10) {
            P = Pb * std::exp(-34.1632 * (h - h_b) / Tb);
        } else {
            double expn = -34.1632 / L;
            P = Pb * std::pow(1.0 + L * (h - h_b) / Tb, expn);
        }
        T = Tb + L * (h - h_b);
        rho = P / (R_AIR * T);
    }

    double get_speed_of_sound(double alt_m) const {
        double rho, T, P;
        get_atmosphere(alt_m, rho, T, P);
        constexpr double gamma = 1.4;
        return std::sqrt(gamma * R_AIR * T);
    }

    double get_gravity(double alt_m) const {
        return G0 * std::pow(R_EARTH / (R_EARTH + alt_m), 2.0);
    }

    double get_aerodynamic_heating(double rho, double vel_mps) const {
        if (vel_mps < 1000.0) return 0.0;
        return K_SUTTON * std::sqrt(rho / RN) * std::pow(vel_mps, 3.0);
    }

    void update_thermal(double dt, double alt_m, double vel_mps) {
        double rho, T, P; get_atmosphere(alt_m, rho, T, P);
        double q = get_aerodynamic_heating(rho, vel_mps);
        double energy_in = q * dt;
        h.pcm_ratio -= energy_in / PCM_MELT_ENERGY;
        h.pcm_ratio = std::clamp(h.pcm_ratio, 0.0, 1.0);
        h.temp += q * dt / 10000.0;
    }

    double get_pcm_remaining() const { return h.pcm_ratio; }
    double get_temperature() const { return h.temp; }
    double get_air_density(double alt_m) const {
        double rho, T, P; get_atmosphere(alt_m, rho, T, P); return rho;
    }
};

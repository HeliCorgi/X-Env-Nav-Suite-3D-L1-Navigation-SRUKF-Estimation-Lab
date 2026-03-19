#pragma once
#include <algorithm>
#include <cmath>

class EnvAdaptive_Core {
private:
    double pcm_remaining = 5000.0;  // Joules
    double temperature = 20.0;      // Celsius
    static constexpr double PCM_MAX = 5000.0;
    static constexpr double MAX_SAFE_TIME = 14400.0;  // 4 hours

public:
    EnvAdaptive_Core() = default;

    void update_thermal(float dt, double heat_input) {
        // Heat dissipation via PCM
        double pcm_absorption = std::min(heat_input * dt * 0.8, pcm_remaining);
        pcm_remaining -= pcm_absorption;

        // Temperature rise (simplified linear model)
        double absorbed = heat_input * dt - pcm_absorption;
        temperature += absorbed * 0.01;  // Empirical scaling

        // Radiative cooling
        temperature *= 0.998f;
        temperature = std::max(temperature, 20.0);
    }

    double get_pcm_remaining() const {
        return std::clamp(pcm_remaining / PCM_MAX, 0.0, 1.0);
    }

    double get_temperature() const {
        return temperature;
    }

    bool is_thermal_safe() const {
        return pcm_remaining > 0.0 && temperature < 300.0;
    }
};

#pragma once

#include <iostream>
#include <fstream>
#include <chrono>
#include <memory>
#include <mutex>
#include <string>

#include "./utility/indicators.hpp"

class Timer {
public:
    static void start() {
        start_time = std::chrono::high_resolution_clock::now();
    }

    static double getCurrentRSS()
    {
        std::ifstream statm("/proc/self/statm");
        size_t size, resident;
        statm >> size >> resident;
        return (resident * getpagesize()) / (1024.0 * 1024.0); // MB
    }

    static void stop(const std::string& message = "Elapsed time") 
    {
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;

        printf("%s :: %.3f s (%.2f Mb)\n", message.c_str(), elapsed.count(), getCurrentRSS());

        //std::cout << std::right
              //<< message << ": "
              //<< std::setw(3) << elapsed.count() << " s"
              //<< " (" << std::setw(4) << getCurrentRSS() << " MB)\n";
    }

    static std::unique_ptr<indicators::ProgressBar> getLoadingBar()
    {
        using namespace indicators;

        auto bar = std::make_unique<ProgressBar>(
                option::BarWidth{50},
                option::Start{"["},
                option::Fill{"■"},
                option::Lead{"■"},
                option::Remainder{"-"},
                option::End{" ]"},
                option::PostfixText{"Computing Reeb space."},
                option::ShowPercentage{true},
                option::ShowElapsedTime{true}
                );

        return bar;
    }

private:
    static inline std::chrono::high_resolution_clock::time_point start_time;
};

#pragma once

#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <chrono>

#include <iostream>
#include <string>
#include <iomanip>
#include <chrono>

class LoadingBar
{
private:
    int width;
    int progress;
    std::string message;
    std::chrono::steady_clock::time_point start;
    char* barBuffer;

public:
    LoadingBar(int w = 50, const std::string& msg = "")
        : width(w), progress(0), message(msg), start(std::chrono::steady_clock::now())
    {
        barBuffer = new char[width + 1];
        barBuffer[width] = '\0';
        std::fill(barBuffer, barBuffer + width, ' ');
    }

    ~LoadingBar()
    {
        delete[] barBuffer;
    }

    void setMessage(const std::string& msg)
    {
        message = msg;
    }

    void update(int percent)
    {
        if (percent == progress) return;

        if (percent < 0) percent = 0;
        if (percent > 100) percent = 100;
        progress = percent;

        print();
    }

private:
    void print()
    {
        int filled = (progress * width) / 100;
        std::fill(barBuffer, barBuffer + filled, '=');
        std::fill(barBuffer + filled, barBuffer + width, ' ');

        auto now = std::chrono::steady_clock::now();
        int elapsed = static_cast<int>(std::chrono::duration_cast<std::chrono::seconds>(now - start).count());
        int minutes = elapsed / 60;
        int seconds = elapsed % 60;

        std::cout << "\r[" << barBuffer << "] " 
                  << std::setw(3) << progress << "%  "
                  << "[" << minutes << "m:" << std::setw(2) << std::setfill('0') << seconds << "s]" 
                  << (message.empty() ? "" : "  " + message);

        std::cout.flush();

        if (progress == 100)
            std::cout << std::endl;

        std::cout << std::setfill(' '); // reset fill char
    }
};

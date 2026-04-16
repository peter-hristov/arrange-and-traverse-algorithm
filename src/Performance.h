#pragma once

#include "./TetMesh.h"
#include "./Arrangement.h"
#include "./ReebSpace2.h"
#include "./LoadingBar.hpp"
#include "./Fiber.h"
#include "./FiberSurface.h"

#include <random>

namespace performance
{

    inline void testInteractiveFiberPerformance(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const int samples, const std::string filename)
    {
        static std::mt19937 gen(std::random_device{}());

        const double epsF = abs(tetMesh.maxF - tetMesh.minF) / 10000.0;
        std::uniform_real_distribution<double> distF(tetMesh.minF + epsF, tetMesh.maxF - epsF);

        const double epsG = abs(tetMesh.maxG - tetMesh.minG) / 10000.0;
        std::uniform_real_distribution<double> distG(tetMesh.minG + epsG, tetMesh.maxG - epsG);


        std::vector<double> timings;

        LoadingBar bar2(40, "Evaluating fiber performance...");

        int iterations = 0;
        
        while (timings.size() < samples)
        {
            const double x = distF(gen);
            const double y = distG(gen);

            try {

                const auto start = std::chrono::high_resolution_clock::now();

                const Fiber fiber = Fiber::computeLabeledFiber(tetMesh, singularArrangement, reebSpace, {x, y}, {});

                const auto end = std::chrono::high_resolution_clock::now();
                const double elapsed = std::chrono::duration<double>(end - start).count();

                if (fiber.components.size() > 0)
                {
                    timings.push_back(elapsed);
                }
            } catch (const std::exception &e) {
                throw;
            }


            if (iterations++ > samples * 1000)
            {
                throw std::runtime_error("Could not sample from distrbution, infinite loop.");
            }

            bar2.update((int)(100.0 * timings.size() / (double)samples));
        }


        const double min = *std::min_element(timings.begin(), timings.end());
        const double max = *std::max_element(timings.begin(), timings.end());
        const double avg = std::accumulate(timings.begin(), timings.end(), 0.0) / timings.size();

        const double sq_sum = std::inner_product(timings.begin(), timings.end(), timings.begin(), 0.0);
        const double stddev = std::sqrt(sq_sum / timings.size() - avg * avg);

        std::cout << std::fixed << std::setprecision(4);
        std::cout << "Timings over " << timings.size() << " samples:\n";
        std::cout << "  Min:    " << min    << " s\n";
        std::cout << "  Max:    " << max    << " s\n";
        std::cout << "  Avg:    " << avg    << " s\n";
        std::cout << "  Stddev: " << stddev << " s\n";


        std::ofstream file(filename);
        if (!file.is_open()) 
        {
            throw std::runtime_error("Could not open file for writing.");
        }

        for (double t : timings)
        {
            file << t << "\n";
        }

        file.close();
    }

    inline void testInteractiveFiberSurfacePerformance(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const int samples, const std::string filename)
    {

        static std::mt19937 gen(std::random_device{}());

        const double epsF = abs(tetMesh.maxF - tetMesh.minF) / 10000.0;
        std::uniform_real_distribution<double> distF(tetMesh.minF + epsF, tetMesh.maxF - epsF);

        const double epsG = abs(tetMesh.maxG - tetMesh.minG) / 10000.0;
        std::uniform_real_distribution<double> distG(tetMesh.minG + epsG, tetMesh.maxG - epsG);


        std::vector<double> timings;

        LoadingBar bar2(40, "Evaluating fiber performance...");

        int iterations = 0;
        
        while (timings.size() < samples)
        {
            const double x1 = distF(gen);
            const double y1 = distG(gen);

            const double x2 = distF(gen);
            const double y2 = distG(gen);

            const std::vector<std::array<double, 2>> controlPoints{{x1, y1}, {x2, y2}};

            try {
                const auto start = std::chrono::high_resolution_clock::now();

                FiberSurface sfMesh  = FiberSurface::constructSegmentedFiberSurface(tetMesh, singularArrangement, reebSpace, controlPoints);

                const auto end = std::chrono::high_resolution_clock::now();
                const double elapsed = std::chrono::duration<double>(end - start).count();


                if (sfMesh.mesh.num_faces() > 0)
                {
                    timings.push_back(elapsed);
                }
            } catch (const std::exception &e) {
                throw;
            }


            if (iterations++ > samples * 1000)
            {
                throw std::runtime_error("Could not sample from distrbution, infinite loop.");
            }

            bar2.update((int)(100.0 * timings.size() / (double)samples));
        }


        const double min = *std::min_element(timings.begin(), timings.end());
        const double max = *std::max_element(timings.begin(), timings.end());
        const double avg = std::accumulate(timings.begin(), timings.end(), 0.0) / timings.size();

        const double sq_sum = std::inner_product(timings.begin(), timings.end(), timings.begin(), 0.0);
        const double stddev = std::sqrt(sq_sum / timings.size() - avg * avg);

        std::cout << std::fixed << std::setprecision(4);
        std::cout << "Timings over " << timings.size() << " samples:\n";
        std::cout << "  Min:    " << min    << " s\n";
        std::cout << "  Max:    " << max    << " s\n";
        std::cout << "  Avg:    " << avg    << " s\n";
        std::cout << "  Stddev: " << stddev << " s\n";


        std::ofstream file(filename);
        if (!file.is_open()) 
        {
            throw std::runtime_error("Could not open file for writing.");
        }

        for (double t : timings)
        {
            file << t << "\n";
        }

        file.close();
    }
}

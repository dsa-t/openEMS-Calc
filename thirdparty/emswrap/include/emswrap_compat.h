#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <fstream>

// physical_constants
constexpr double C0   = 299792458.0; // m/s
constexpr double PI   = 3.14159265358979323846;
constexpr double M_PI = PI;
constexpr double MUE0 = 4e-7 * PI;     // N/A^2
constexpr double EPS0 = 1.0 / (MUE0 * C0 * C0); // F/m

// free space wave impedance (Ohm)
constexpr double Z0 = 376.73031346177065; // std::sqrt(MUE0 / EPS0)


template<typename T>
static std::vector<double> linspace(T start_in, T end_in, int num_in)
{
    std::vector<double> linspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for (int i = 0; i < num - 1; ++i)
        linspaced.push_back(start + delta * i);

    linspaced.push_back(end);

    return linspaced;
}

static std::vector<double> uniqueSorted(const std::vector<double> &vec) {
    std::vector<double> sorted = vec;
    std::sort(sorted.begin(), sorted.end());
    auto last = std::unique(sorted.begin(), sorted.end());
    sorted.erase(last, sorted.end());
    return sorted;
}

// Supports int or char
static int CheckNyDir(int ny) {
  if (ny == 0 || ny == 1 || ny == 2)
      return ny;
  else if (ny == 'x' || ny == 'y' || ny == 'z')
      return ny == 'x' ? 0 : ny == 'y' ? 1 : 2;
  else if (ny == 'r' || ny == 'a' || ny == 'z')
      return ny == 'r' ? 0 : ny == 'a' ? 1 : 2;
  else
      throw std::invalid_argument("Invalid direction");
}

static bool Check_Array_Equal(const std::vector<double>& a, const std::vector<double>& b, double tol, bool relative = false) {
    if (a.size() != b.size())
        return false;
    if (tol == 0.0) {
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] != b[i])
                return false;
        }
        return true;
    }
    for (size_t i = 0; i < a.size(); ++i) {
        double diff = relative 
            ? (a[i] != 0.0 ? std::abs((a[i] - b[i]) / a[i]) : std::abs(a[i] - b[i]))
            : std::abs(a[i] - b[i]);
        if (diff >= tol)
            return false;
    }
    return true;
}

static void PyPlot(const std::string &baseName, const std::string &codePart) {
    std::string imports = R"(
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
)";
    std::string fullCode = imports + codePart;

    std::string filename = "plot_" + baseName + ".py";
    std::ofstream pyFile(filename);
    if (!pyFile)
        throw std::runtime_error("Failed to open file: " + filename);
        
    pyFile << fullCode;
    pyFile.close();

    std::string command = "python " + filename;
    std::system(command.c_str());
}
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


static std::vector<double> linspace(double start_in, double end_in, int num_in)
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

// Finds the index of value in arr closest to specified value
static int interp1_nearest(double value, const std::vector<double>& arr)
{
    double dist = DBL_MAX;
    int idx = -1;
    for (int i = 0; i < arr.size(); ++i) {
        double newDist = std::abs(value - arr[i]);
        if (newDist < dist) {
            dist = newDist;
            idx = i;
        }
    }

    return idx;
}

static std::vector<double> diff(std::vector<double>& arr) {
    std::vector<double> ret;

    if (arr.size() < 2)
        return ret;

    ret.resize(arr.size() - 1);

    for (size_t i = 0; i < arr.size() - 1; i++) {
        ret[i] = arr[i + 1] - arr[i];
    }

    return ret;
}

class num_vector : public std::vector<double> {
public:
    num_vector() = default;
    num_vector(std::initializer_list<double> il) : std::vector<double>(il) {}
    num_vector(const std::vector<double>& vec) : std::vector<double>(vec) {}
    num_vector(std::vector<double>&& vec) : std::vector<double>(std::move(vec)) {}

    num_vector(const num_vector& other) : std::vector<double>(other) {}
    num_vector(num_vector&& other) : std::vector<double>(std::move(other)) {}

    template<typename Arg, typename... Args>
    num_vector(Arg arg0, Args&&... args) {
        append(arg0, std::forward<Args>(args)...);
    }

    num_vector& append(double value) {
        this->push_back(value);
        return *this;
    }

    num_vector& append(const std::vector<double>& values) {
        this->insert(this->end(), values.begin(), values.end());
        return *this;
    }

    num_vector& append(const num_vector& values) {
        this->insert(this->end(), values.begin(), values.end());
        return *this;
    }

    template<typename Func>
    num_vector& append_transform(Func fn, const std::vector<double>& values) {
        for (const double &value : values)
            this->push_back(fn(value));
        return *this;
    }

    template <typename Arg, typename... Args>
    num_vector& append(Arg arg0, Args &&...args)
    {
        append(arg0);
        append(std::forward<Args>(args)...);
        return *this;
    }

    num_vector operator-() const
    {
        num_vector result;
        for (const double &value : *this)
            result.push_back(-value);
        return result;
    }
};

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
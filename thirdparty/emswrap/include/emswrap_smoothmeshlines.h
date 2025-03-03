#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>

namespace SmoothMesh {

#pragma optimize("",off)

// Internal function to check if the mesh lines are symmetric
static int CheckSymmetry(const std::vector<double>& lines, double tolerance = 1e-10) {
    if (lines.size() <= 2) return 0;
    double line_range = lines.back() - lines.front();
    double center = 0.5 * (lines.back() + lines.front());
    size_t NP = lines.size();

    // Check each pair for symmetry
    for (size_t n = 0; n < NP / 2; ++n) {
        double left = lines[n];
        double right = lines[NP - n - 1];
        if (std::abs((center - left) - (right - center)) > line_range * tolerance) {
            return 0;
        }
    }

    // Check central element if odd
    if (NP % 2 == 1) {
        if (std::abs(lines[NP / 2] - center) > line_range * tolerance) {
            return 0;
        }
    }

    return (NP % 2 == 0) ? 2 : 1;
}

struct CheckMeshResult {
    int ec;
    std::vector<size_t> pos;
    std::vector<int> error_type;
};

static CheckMeshResult CheckMesh(const std::vector<double>& lines, double min_res, double max_res, double ratio, bool be_quiet = false) {
    CheckMeshResult result{0, {}, {}};

    if (lines.size() < 2) return result;

    // Compute differences between consecutive lines.
    std::vector<double> diff_lines;
    for (size_t i = 0; i < lines.size() - 1; ++i) {
        diff_lines.push_back(lines[i+1] - lines[i]);
    }

    // Check if any difference exceeds max_res (error type 1).
    for (size_t i = 0; i < diff_lines.size(); ++i) {
        if (diff_lines[i] > max_res) {
            if (!be_quiet)
                std::cerr << "Warning: found resolution larger than max_res at segment " << i << std::endl;
            result.pos.push_back(i);
            result.ec++;
            result.error_type.push_back(1);
        }
    }

    // Check if any difference is smaller than min_res (error type 2).
    for (size_t i = 0; i < diff_lines.size(); ++i) {
        if (diff_lines[i] < min_res) {
            if (!be_quiet)
                std::cerr << "Warning: found resolution smaller than min_res at segment " << i << std::endl;
            result.pos.push_back(i);
            result.ec++;
            result.error_type.push_back(2);
        }
    }

    // Check for local ratio violations between consecutive intervals.
    // For index i, compare diff_lines[i+1] / diff_lines[i]
    for (size_t i = 0; i < diff_lines.size() - 1; ++i) {
        double ratio_val = diff_lines[i+1] / diff_lines[i];
        if (ratio_val > ratio * 1.01) {
            if (!be_quiet) {
                std::cerr << "Warning: resolution increase too high at segments " 
                          << i << ", " << i+1 << ", " << i+2 
                          << " with ratio " << ratio_val << " > " << ratio << std::endl;
            }
            result.pos.push_back(i+1);
            result.ec++;
            result.error_type.push_back(3);
        }
        if (ratio_val < (1.0/ratio) / 1.01) {
            if (!be_quiet) {
                std::cerr << "Warning: resolution decrease too steep at segments " 
                          << i << ", " << i+1 << ", " << i+2 
                          << " with ratio " << ratio_val << " < " << 1.0/ratio << std::endl;
            }
            result.pos.push_back(i);
            result.ec++;
            result.error_type.push_back(4);
        }
    }
    return result;
}

// Internal function to deduplicate and remove close points
static std::vector<double> Unique(std::vector<double> l, double tol = 1e-7) {
    if (l.empty()) return l;

    // Sort and remove duplicates
    std::sort(l.begin(), l.end());
    auto last = std::unique(l.begin(), l.end());
    l.erase(last, l.end());

    if (l.size() <= 1) return l;

    // Compute differences
    std::vector<double> dl;
    for (size_t i = 0; i < l.size() - 1; ++i) {
        dl.push_back(l[i+1] - l[i]);
    }

    double mean = 0.0;
    for (double d : dl) mean += d;
    mean /= dl.size();

    // Collect indices to erase
    std::vector<size_t> idx;
    for (size_t i = 0; i < dl.size(); ++i) {
        if (dl[i] < mean * tol) {
            idx.push_back(i);
        }
    }

    // Erase from the end to avoid shifting issues
    std::sort(idx.rbegin(), idx.rend());
    for (size_t i : idx) {
        l.erase(l.begin() + i);
    }

    return l;
}

// Helper function for one-sided tapering
static std::vector<double> one_side_taper(double start_res, double ratio, double max_res, double rng) {
    std::vector<double> l = {0.0};
    double res = start_res;
    double pos = 0.0;
    int N = 0;

    while (res < max_res && pos < rng) {
        res *= ratio;
        pos += res;
        N++;
    }

    if (pos > rng) {
        // Scale positions to fit rng
        std::vector<double> tmp(N + 1, 0.0);
        double sum = 0.0;
        for (int n = 0; n < N; ++n) {
            sum += start_res * std::pow(ratio, n + 1);
            tmp[n + 1] = sum;
        }
        double scale = rng / pos;
        for (auto& val : tmp) val *= scale;
        return tmp;
    } else {
        // Adjust ratio to reach max_res in N steps
        double actual_ratio = std::exp((std::log(max_res) - std::log(start_res)) / N);
        l.clear();
        l.push_back(0.0);
        pos = 0.0;
        res = start_res;
        for (int n = 0; n < N; ++n) {
            res *= actual_ratio;
            pos += res;
            l.push_back(pos);
        }
        // Add max_res steps until end
        while (pos < rng) {
            pos += max_res;
            l.push_back(pos);
        }
        // Scale to fit rng
        double scale = rng / l.back();
        for (auto& val : l) val *= scale;
        return l;
    }
}

// Generate a smooth range of points between start and stop
static std::vector<double> SmoothRange(double start, double stop, double start_res, double stop_res, double max_res, double ratio) {
    assert(ratio > 1.0);
    double rng = stop - start;
    std::vector<double> result;

    // Handle very small range
    if (rng < max_res && rng < start_res * ratio && rng < stop_res * ratio) {
        result.push_back(start);
        result.push_back(stop);
        return Unique(result);
    }

    // Both resolutions are sufficient
    if (start_res >= (max_res / ratio) && stop_res >= (max_res / ratio)) {
        int N = std::ceil(rng / max_res);
        for (int i = 0; i <= N; ++i) {
            result.push_back(start + i * (rng / N));
        }
        return Unique(result);
    }

    // Taper start
    if (start_res < (max_res / ratio) && stop_res >= (max_res / ratio)) {
        auto tmp = one_side_taper(start_res, ratio, max_res, rng);
        for (auto& val : tmp) val += start;
        tmp = Unique(tmp);
        if (tmp.front() != start) tmp.insert(tmp.begin(), start);
        if (tmp.back() != stop) tmp.push_back(stop);
        return tmp;
    }

    // Taper stop
    if (start_res >= (max_res / ratio) && stop_res < (max_res / ratio)) {
        auto tmp = one_side_taper(stop_res, ratio, max_res, rng);
        std::reverse(tmp.begin(), tmp.end());
        for (auto& val : tmp) val = stop - val;
        tmp = Unique(tmp);
        if (tmp.front() != start) tmp.insert(tmp.begin(), start);
        if (tmp.back() != stop) tmp.push_back(stop);
        return tmp;
    }

    // Both sides need tapering
    // Calculate required steps for each side
    double pos1 = 0.0, pos2 = 0.0;
    int N1 = 0, N2 = 0;
    double res1 = start_res, res2 = stop_res;

    while (res1 < max_res) { res1 *= ratio; pos1 += res1; N1++; }
    double ratio1 = std::exp((std::log(max_res) - std::log(start_res)) / N1);
    pos1 = start_res * (std::pow(ratio1, N1 + 1) - ratio1) / (ratio1 - 1);

    while (res2 < max_res) { res2 *= ratio; pos2 += res2; N2++; }
    double ratio2 = std::exp((std::log(max_res) - std::log(stop_res)) / N2);
    pos2 = stop_res * (std::pow(ratio2, N2 + 1) - ratio2) / (ratio2 - 1);

    if (pos1 + pos2 < rng) {
        // Construct left and right parts
        std::vector<double> left = {0.0};
        res1 = start_res;
        for (int n = 0; n < N1; ++n) {
            left.push_back(left.back() + res1);
            res1 *= ratio1;
        }
        std::vector<double> right = {0.0};
        res2 = stop_res;
        for (int n = 0; n < N2; ++n) {
            right.push_back(right.back() + res2);
            res2 *= ratio2;
        }
        double remaining = rng - pos1 - pos2;
        int N_center = std::ceil(remaining / max_res);
        for (int i = 0; i < N_center; ++i) {
            left.push_back(left.back() + max_res);
        }
        // Combine and scale
        double total_length = left.back() + right.back();
        std::vector<double> combined;
        combined.insert(combined.end(), left.begin(), left.end());
        for (auto it = right.rbegin(); it != right.rend(); ++it) {
            combined.push_back(total_length - *it);
        }
        for (auto& val : combined) val = start + val * rng / total_length;
        combined = Unique(combined);
        if (combined.front() != start) combined.insert(combined.begin(), start);
        if (combined.back() != stop) combined.push_back(stop);
        return combined;
    } else {
        // Grow both sides until they meet
        std::vector<double> left = {0.0};
        std::vector<double> right = {0.0};
        res1 = start_res;
        res2 = stop_res;
        while (left.back() + right.back() < rng) {
            if (res1 < res2) {
                left.push_back(left.back() + res1);
                res1 *= ratio;
            } else {
                right.push_back(right.back() + res2);
                res2 *= ratio;
            }
        }
        // Combine and scale
        double total_length = left.back() + right.back();
        std::vector<double> combined;
        combined.insert(combined.end(), left.begin(), left.end());
        for (auto it = right.rbegin(); it != right.rend(); ++it) {
            combined.push_back(total_length - *it);
        }
        for (auto& val : combined) val = start + val * rng / total_length;
        combined = Unique(combined);
        if (combined.front() != start) combined.insert(combined.begin(), start);
        if (combined.back() != stop) combined.push_back(stop);
        return combined;
    }
}

// Main function to smooth mesh lines
static std::vector<double> SmoothMeshLines(const std::vector<double>& lines, double max_res, double ratio = 1.5) {
    std::vector<double> out_l = Unique(lines);
    int sym = CheckSymmetry(out_l);
    double center = 0.0;

    if (sym == 1 || sym == 2) {
        center = 0.5 * (out_l.front() + out_l.back());
        size_t new_size = (sym == 1) ? (out_l.size() / 2 + 1) : (out_l.size() / 2);
        out_l = std::vector<double>(out_l.begin(), out_l.begin() + new_size);
    }

    bool changed;
    do {
        changed = false;
        std::vector<double> dl;
        for (size_t i = 0; i < out_l.size() - 1; ++i) {
            dl.push_back(out_l[i+1] - out_l[i]);
        }

        std::vector<size_t> indices_above;
        for (size_t i = 0; i < dl.size(); ++i) {
            if (dl[i] > max_res) {
                indices_above.push_back(i);
            }
        }
        if (indices_above.empty()) break;

        // Mark small intervals and find the smallest large interval
        double max_dl = *std::max_element(dl.begin(), dl.end());
        for (auto& d : dl) {
            if (d <= max_res) d = max_dl * 2;
        }
        size_t idx = std::min_element(dl.begin(), dl.end()) - dl.begin();

        double start_res = (idx > 0) ? dl[idx-1] : max_res;
        double stop_res = (idx < dl.size()-1) ? dl[idx+1] : max_res;

        double start = out_l[idx];
        double stop = out_l[idx+1];
        auto new_points = SmoothRange(start, stop, start_res, stop_res, max_res, ratio);

        size_t prev_size = out_l.size();
        out_l.insert(out_l.end(), new_points.begin(), new_points.end());
        out_l = Unique(out_l);
        changed = (out_l.size() != prev_size);
    } while (changed);

    // Handle symmetry
    if (sym == 1) {
        std::vector<double> mirrored;
        for (auto it = out_l.rbegin() + 1; it != out_l.rend(); ++it) {
            mirrored.push_back(2 * center - *it);
        }
        out_l.insert(out_l.end(), mirrored.begin(), mirrored.end());
        out_l = Unique(out_l);
    } else if (sym == 2) {
        double start = out_l.back();
        double stop = 2 * center - start;
        auto new_points = SmoothRange(start, stop, out_l.back() - out_l[out_l.size()-2], out_l.back() - out_l[out_l.size()-2], max_res, ratio);
        out_l.insert(out_l.end(), new_points.begin(), new_points.end());
        std::vector<double> mirrored;
        for (auto it = out_l.rbegin() + 1; it != out_l.rend(); ++it) {
            mirrored.push_back(2 * center - *it);
        }
        out_l.insert(out_l.end(), mirrored.begin(), mirrored.end());
        out_l = Unique(out_l);
    }

    return out_l;
}

} // namespace SmoothMesh
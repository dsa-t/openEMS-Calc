#pragma once

#include <emswrap_smoothmeshlines.h>

#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>

namespace SmoothMesh
{

#pragma optimize("",off)

    struct Gap
    {
        double start_res = -1;
        double stop_res = -1;
        std::vector<double> lines;
    };

    static std::vector<Gap> calc_gaps(
        const std::vector<double> &lines,
        const std::vector<Gap> &old_gaps,
        double max_res)
    {
        std::vector<Gap> gaps(old_gaps.size());
        std::vector<double> temp_lines = lines;
        for (const auto &gap : old_gaps)
        {
            temp_lines.insert(temp_lines.end(), gap.lines.begin(), gap.lines.end());
        }
        std::sort(temp_lines.begin(), temp_lines.end());
        temp_lines = Unique(temp_lines);

        if (temp_lines.size() < 2)
            return gaps;

        // Calculate start_res and stop_res for each gap
        for (size_t n = 0; n < lines.size() - 1; ++n)
        {
            auto lower = std::lower_bound(temp_lines.begin(), temp_lines.end(), lines[n]);
            auto upper = std::upper_bound(temp_lines.begin(), temp_lines.end(), lines[n + 1]);

            if (n == 0)
            {
                gaps[n].start_res = std::numeric_limits<double>::infinity();
                if (upper - lower >= 1)
                {
                    gaps[n].stop_res = std::min<double>(max_res, (*(lower + 1) - *lower));
                }
                else
                {
                    gaps[n].stop_res = max_res;
                }
            }
            else if (n == lines.size() - 2)
            {
                gaps[n].stop_res = std::numeric_limits<double>::infinity();
                if (upper - lower >= 1)
                {
                    gaps[n].start_res = std::min<double>(max_res, *lower - *(lower - 1));
                }
                else
                {
                    gaps[n].start_res = max_res;
                }
            }
            else
            {
                if (lower != temp_lines.begin() && upper != temp_lines.end())
                {
                    gaps[n].start_res = std::min<double>(max_res, *lower - *(lower - 1));
                    gaps[n].stop_res = std::min<double>(max_res, *(upper + 1) - *upper);
                }
            }
        }
        return gaps;
    }

    static double calc_start_res(size_t pos, const std::vector<double> &lines, const std::vector<Gap> &gaps)
    {
        if (pos >= lines.size() - 1)
            return -1;
        std::vector<double> temp_lines = lines;
        for (const auto &gap : gaps)
        {
            temp_lines.insert(temp_lines.end(), gap.lines.begin(), gap.lines.end());
        }
        std::sort(temp_lines.begin(), temp_lines.end());
        auto it = std::lower_bound(temp_lines.begin(), temp_lines.end(), lines[pos]);
        if (it == temp_lines.begin())
            return -1;
        return *it - *(it - 1);
    }

    static double calc_stop_res(size_t pos, const std::vector<double> &lines, const std::vector<Gap> &gaps)
    {
        if (pos >= lines.size() - 1)
            return -1;
        std::vector<double> temp_lines = lines;
        for (const auto &gap : gaps)
        {
            temp_lines.insert(temp_lines.end(), gap.lines.begin(), gap.lines.end());
        }
        std::sort(temp_lines.begin(), temp_lines.end());
        auto it = std::upper_bound(temp_lines.begin(), temp_lines.end(), lines[pos + 1]);
        if (it == temp_lines.end())
            return -1;
        return *it - *(it - 1);
    }


    static std::vector<double> SmoothMeshLines2(
        const std::vector<double>& in_lines,
        double max_res,
        double ratio = 1.3,
        bool check_mesh = true,
        double allowed_max_ratio = -1.0
    ) {
        std::vector<double> lines = Unique(in_lines);
        if (lines.size() < 2) return lines;
    
        if (allowed_max_ratio < 0) allowed_max_ratio = ratio * 1.25;
    
        // Special case for exactly two lines
        if (lines.size() == 2) {
            auto new_points = SmoothRange(lines[0], lines[1], max_res, max_res, max_res, ratio);
            lines.insert(lines.end(), new_points.begin(), new_points.end());
            lines = Unique(lines);
            CheckMeshResult res = CheckMesh(lines, 0, max_res, allowed_max_ratio, true);
            return lines;
        }
    
        std::vector<Gap> gaps(lines.size() - 1);
        bool changed;
        do {
            changed = false;
            auto new_gaps = calc_gaps(lines, gaps, max_res);
    
            // Find first gap with changed start_res or stop_res
            size_t index = 0;
            for (; index < new_gaps.size(); ++index) {
                if (index >= gaps.size()) break;
                if (new_gaps[index].start_res != gaps[index].start_res ||
                    new_gaps[index].stop_res != gaps[index].stop_res) {
                    break;
                }
            }
            if (index >= new_gaps.size()) break; // No changes
    
            // Process the gap
            double start_res = std::min<double>(max_res, calc_start_res(index, lines, gaps));
            double stop_res = std::min<double>(max_res, calc_stop_res(index, lines, gaps));
    
            auto new_points = SmoothRange(
                lines[index], lines[index + 1],
                start_res, stop_res, max_res, ratio
            );
            gaps[index].lines = new_points;
            gaps[index].start_res = start_res;
            gaps[index].stop_res = stop_res;
    
            // Merge lines
            std::vector<double> merged = lines;
            for (const auto& gap : gaps) {
                merged.insert(merged.end(), gap.lines.begin(), gap.lines.end());
            }
            merged = Unique(merged);
            changed = (merged.size() != lines.size());
            lines = merged;
    
        } while (changed);
    
        if (check_mesh) {
            CheckMesh(lines, 0, max_res, allowed_max_ratio, true);
        }
    
        return lines;
    }

}
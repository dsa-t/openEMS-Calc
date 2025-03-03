#include "emswrap_csrectgrid.h"
#include <emswrap_compat.h>
#include <emswrap_csrectgrid.h>
#include <emswrap_smoothmeshlines.h>
#include <emswrap_smoothmeshlines2.h>

#include <vector>
#include <cmath>

// Implementation of _CSRectGrid member functions.

_CSRectGrid::_CSRectGrid(CSRectGrid* grid)
    : _ptr(grid)
{
}

_CSRectGrid::~_CSRectGrid()
{
}

void _CSRectGrid::SetMeshType(CoordinateSystem cs_type) {
    _ptr->SetMeshType(cs_type);
}

int _CSRectGrid::GetMeshType() const {
    return _ptr->GetMeshType();
}

void _CSRectGrid::ClearLines(int ny) {
    ny = CheckNyDir(ny);
    _ptr->ClearLines(ny);
}

void _CSRectGrid::AddDiscLine(int ny, double line) {
    ny = CheckNyDir(ny);
    _ptr->AddDiscLine(ny, line);
}

unsigned int _CSRectGrid::GetQtyLines(int ny) const {
    ny = CheckNyDir(ny);
    return _ptr->GetQtyLines(ny);
}

double _CSRectGrid::GetLine(int ny, int idx) const {
    ny = CheckNyDir(ny);
    return _ptr->GetLine(ny, idx);
}

#pragma optimize("",off)
std::vector<double> _CSRectGrid::GetLines(int ny, bool do_sort)
{
    ny = CheckNyDir(ny);

    size_t qty = _ptr->GetQtyLines(ny);
    std::vector<double> lines(qty);

    for (size_t i = 0; i < qty; ++i)
        lines[i] = _ptr->GetLine(ny, i);

    if (do_sort)
        std::sort(lines.begin(), lines.end());

    return lines;
}
#pragma optimize("",on)

void _CSRectGrid::clear() {
    _ptr->clear();
}

void _CSRectGrid::SetDeltaUnit(double unit) {
    _ptr->SetDeltaUnit(unit);
}

double _CSRectGrid::GetDeltaUnit() const {
    return _ptr->GetDeltaUnit();
}

void _CSRectGrid::Sort(int ny) {
    ny = CheckNyDir(ny);
    _ptr->Sort(ny);
}

double _CSRectGrid::Snap2LineNumber(int ny, double value, bool& inside) {
    ny = CheckNyDir(ny);
    return _ptr->Snap2LineNumber(ny, value, inside);
}

std::array<std::array<double, 3>, 2> _CSRectGrid::GetSimArea() const {
    std::array<std::array<double, 3>, 2> simArea;
    double* _bb = _ptr->GetSimArea();
    for (int n = 0; n < 3; ++n) {
        simArea[0][n] = _bb[2 * n];
        simArea[1][n] = _bb[2 * n + 1];
    }
    return simArea;
}

bool _CSRectGrid::isValid() const {
    return _ptr->isValid();
}

void _CSRectGrid::SmoothMeshLines(int ny, double max_res, double ratio) {
    if (ny < 0)
    {
        // Smooth all directions (assuming three directions)
        for (int n = 0; n < 3; ++n)
            this->SmoothMeshLines(n, max_res, ratio);
        return;
    }

    ny = CheckNyDir(ny);
    SetLines(ny, SmoothMesh::SmoothMeshLines(GetLines(ny), max_res, ratio));
}

#pragma optimize("",off)

std::vector<double> _CSRectGrid::SmoothMeshLines(int ny, const std::vector<double> &lines, double max_res, double ratio, int recursive, bool CheckMesh, double allowed_max_ratio)
{
    std::vector<double> lines_new = lines;
    if (lines_new.size() < 2)
        return lines_new;

    if (allowed_max_ratio == -1.0)
        allowed_max_ratio = ratio * 1.25;

    // sort and unique
    std::sort(lines_new.begin(), lines_new.end());
    lines_new.erase(std::unique(lines_new.begin(), lines_new.end()), lines_new.end());

    // First pass: add new lines where spacing is too large
    std::vector<double> diffLines;
    std::vector<size_t> indices;
    for (size_t i = 1; i < lines_new.size(); ++i)
    {
        double d = lines_new[i] - lines_new[i - 1];
        diffLines.push_back(d);
        if (d > 1.001 * max_res)
            indices.push_back(i - 1);
    }

    std::vector<double> addLines;
    for (size_t i : indices)
    {
        double start_res = (i == 0) ? max_res : (lines_new[i] - lines_new[i - 1]);
        double stop_res = ((i + 1) >= lines_new.size() - 1) ? max_res : (lines_new[i + 2] - lines_new[i + 1]);
        // SmoothRange returns a vector of intermediate lines between two given lines.
        std::vector<double> newLines = SmoothMesh::SmoothRange(lines_new[i], lines_new[i + 1], start_res, stop_res, max_res, ratio);
        addLines.insert(addLines.end(), newLines.begin(), newLines.end());
    }

    lines_new.insert(lines_new.end(), addLines.begin(), addLines.end());
    std::sort(lines_new.begin(), lines_new.end());
    lines_new.erase(std::unique(lines_new.begin(), lines_new.end()), lines_new.end());

    // Second pass: relax ratio and check mesh
    double ratio_relax = ratio + (ratio - 1.0);
    SmoothMesh::CheckMeshResult cms = SmoothMesh::CheckMesh(lines_new, 0, max_res, ratio_relax, 1);
    addLines.clear();
    diffLines.clear();
    for (size_t i = 1; i < lines_new.size(); ++i)
        diffLines.push_back(lines_new[i] - lines_new[i - 1]);

    for (size_t pos : cms.pos)
    {
        double start_res = (pos > 0) ? diffLines[pos - 1] : diffLines[pos];
        double stop_res = (pos >= diffLines.size() - 1) ? diffLines.back() : diffLines[pos + 1];
        double max_res_R = std::max(start_res, stop_res) / (2 * ratio);
        std::vector<double> newLines = SmoothMesh::SmoothRange(lines_new[pos], lines_new[pos + 1], start_res, stop_res, max_res_R, ratio);
        addLines.insert(addLines.end(), newLines.begin(), newLines.end());
    }

    lines_new.insert(lines_new.end(), addLines.begin(), addLines.end());
    std::sort(lines_new.begin(), lines_new.end());
    lines_new.erase(std::unique(lines_new.begin(), lines_new.end()), lines_new.end());

    // Recursive refinement if requested
    for (int i = 0; i < recursive; ++i)
    {
        size_t oldSize = lines_new.size();
        std::vector<double> newLinesVec = this->SmoothMeshLines(ny, lines_new, max_res, ratio, 0, CheckMesh, allowed_max_ratio);
        if (newLinesVec.size() == oldSize)
            break;
        lines_new = newLinesVec;
    }

    // Final mesh check if required
    if (CheckMesh)
        SmoothMesh::CheckMesh(lines_new, 0, max_res, allowed_max_ratio, 0);

    SetLines(ny, lines_new);

    return lines_new;
}

std::vector<double> _CSRectGrid::SmoothMeshLines2(int ny, const std::vector<double> &lines, double max_res, double ratio, bool CheckMesh, double allowed_max_ratio)
{
    auto lines_new = SmoothMesh::SmoothMeshLines2(lines, max_res, ratio, CheckMesh, allowed_max_ratio);
    SetLines(ny, lines_new);

    return lines_new;
}

void _CSRectGrid::SetLines(int ny, const std::vector<double>& lines) {
    ny = CheckNyDir(ny);
    
    _ptr->ClearLines(ny);
    for (const double line : lines)
        _ptr->AddDiscLine(ny, line);
}

void _CSRectGrid::AddLine(int ny, double line) {
    ny = CheckNyDir(ny);
    _ptr->AddDiscLine(ny, line);
}

void _CSRectGrid::AddLine(int ny, const std::vector<double>& lines) {
    ny = CheckNyDir(ny);

    if (lines.empty())
        throw std::invalid_argument("AddLine: \"lines\" must not be empty");
    for (const double line : lines)
        _ptr->AddDiscLine(ny, line);
}

#pragma optimize("",on)
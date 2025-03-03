#pragma once

#include <CSRectGrid.h>
#include <array>

class _CSRectGrid {
public:
    // Constructors and destructor.
    _CSRectGrid(CSRectGrid* grid);
    ~_CSRectGrid();

    // Member functions.
    void SetMeshType(CoordinateSystem cs_type);
    int GetMeshType() const;
    void ClearLines(int ny);
    void AddDiscLine(int ny, double line);
    unsigned int GetQtyLines(int ny) const;
    double GetLine(int ny, int idx) const;
    std::vector<double> GetLines(int ny, bool do_sort = false);
    void clear();
    void SetDeltaUnit(double unit);
    double GetDeltaUnit() const;
    void Sort(int ny);
    double Snap2LineNumber(int ny, double value, bool& inside);
    std::array<std::array<double, 3>, 2> GetSimArea() const;
    bool isValid() const;

    // Python-like API
    void SmoothMeshLines(int ny, double max_res, double ratio = 1.5);
    void SetLines(int ny, const std::vector<double> &lines);
    void AddLine(int ny, double line);
    void AddLine(int ny, const std::vector<double> &lines);

    // MatLab-like API
    std::vector<double> SmoothMeshLines(int ny,
                                        const std::vector<double> &lines,
                                        double max_res,
                                        double ratio = 1.3,
                                        int recursive = 0,
                                        bool CheckMesh = true,
                                        double allowed_max_ratio = -1.0);

    std::vector<double> SmoothMeshLines2(int ny,
                                         const std::vector<double> &lines,
                                         double max_res,
                                         double ratio = 1.3,
                                         bool CheckMesh = true,
                                         double allowed_max_ratio = -1.0);

    CSRectGrid* operator->() { return _ptr; };

protected:
    CSRectGrid* _ptr;
};

// Wrapper helpers

inline _CSRectGrid emswrap(CSRectGrid *ptr)
{
    return _CSRectGrid(ptr);
}
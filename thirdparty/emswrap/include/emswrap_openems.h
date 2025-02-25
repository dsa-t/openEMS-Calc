#pragma once

#include <openems.h>
#include <emswrap_nf2ff.h>

#include <optional>

class _CSX;

class _openEMS : public openEMS {
public:
    // Constructor/Destructor
    _openEMS::_openEMS(unsigned int NrTS = 1e9, double EndCriteria = 1e-5, double MaxTime = 0);
    virtual ~_openEMS();

    // Core FDTD control
    void SetCoordSystem(int val);

    // Boundary conditions
    void SetBoundaryCond(const std::vector<std::string>& BC);
    enum BOUNDARY_TYPE { PEC, PMC, MUR, PML };

    void AddLumpedPort(int port_nr, double R, 
                      const std::array<double, 3>& start,
                      const std::array<double, 3>& stop,
                      char p_dir, bool excite);
                      
    _nf2ff CreateNF2FFBox(const std::string& name = "nf2ff",
        const std::optional<std::array<double, 3>>& start = std::nullopt,
        const std::optional<std::array<double, 3>>& stop = std::nullopt);

    // Grid management
    void AddEdges2Grid(const std::string &dirs,
                       _CSX *csx,
                       const std::vector<void *> &primitives = {});

    int Run(const std::string& sim_path, bool cleanup = false, bool setup_only = false, const std::vector<std::string>& options = {});

private:
    std::vector<int> boundary_types;  // Stores BC types for 6 directions
};

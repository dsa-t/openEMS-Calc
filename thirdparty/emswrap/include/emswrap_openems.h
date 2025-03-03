#pragma once

#include <openems.h>
#include <emswrap_nf2ff.h>
#include <emswrap_ports.h>

#include <optional>

class _CSX;

class _openEMS : public openEMS {
public:
    // Constructor/Destructor
    _openEMS(unsigned int NrTS = 1e9, double EndCriteria = 1e-5, double MaxTime = 0);
    virtual ~_openEMS();

    // Core FDTD control
    void SetCoordSystem(int val);

    // Boundary conditions
    void SetBoundaryCond(const std::vector<std::string>& BC);
    enum BOUNDARY_TYPE { PEC, PMC, MUR, PML };

    int Run(const std::string& sim_path, bool cleanup = false, bool setup_only = false, const std::vector<std::string>& options = {});

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

    // Port addition functions (Python-like API)
    WaveguidePort AddWaveGuidePort(int port_nr,
                                   const std::array<double, 3> &start,
                                   const std::array<double, 3> &stop,
                                   char p_dir,
                                   const std::array<std::string, 3> &E_WG_func,
                                   const std::array<std::string, 3> &H_WG_func,
                                   double kc,
                                   bool excite = false);

    RectWGPort AddRectWaveGuidePort(int port_nr,
                                    const std::array<double, 3> &start,
                                    const std::array<double, 3> &stop,
                                    char p_dir,
                                    double a,
                                    double b,
                                    const std::string &mode_name,
                                    bool excite = false);

    MSLPort AddMSLPort(int port_nr,
                       CSProperties *metal_prop,
                       const std::array<double, 3> &start,
                       const std::array<double, 3> &stop,
                       char prop_dir,
                       char exc_dir,
                       bool excite = false);

private:
    std::vector<int> boundary_types;  // Stores BC types for 6 directions
};

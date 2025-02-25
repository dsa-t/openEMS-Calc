#include "emswrap_openems.h"
#include <emswrap_openems.h>
#include <emswrap_nf2ff.h>
#include <emswrap_csrectgrid.h>

#include <openems.h>
#include <ContinuousStructure.h>

#include <iostream>
#include <vector>
#include <string>
#include <cstring>


//-----------------------------------------------------------------------------
// Implementation of the _openEMS class
//-----------------------------------------------------------------------------

// Constructor/Destructor
_openEMS::_openEMS(unsigned int NrTS, double EndCriteria, double MaxTime) : boundary_types(6, 0) {
    SetNumberOfTimeSteps(NrTS);
    SetEndCriteria(EndCriteria);
    SetMaxTime(MaxTime);
}

_openEMS::~_openEMS() {
    // Cleanup resources if needed
}

// Boundary Conditions
void _openEMS::SetBoundaryCond(const std::vector<std::string>& BC) {
    if (BC.size() != 6)
        throw std::runtime_error("Invalid boundary condition size!");

    for (size_t n = 0; n < BC.size(); ++n) {
        // Try to interpret as an integer value
        try {
            int bc_val = std::stoi(BC[n]);
            Set_BC_Type(n, bc_val);
            continue;
        } catch (const std::invalid_argument&) {
            // Not an integer, continue with string checks
        }

        if (BC[n] == "PEC" || BC[n] == "PMC" || BC[n] == "MUR") {
            const std::vector<std::string> types = {"PEC", "PMC", "MUR"};
            auto it = std::find(types.begin(), types.end(), BC[n]);
            Set_BC_Type(n, static_cast<int>(std::distance(types.begin(), it)));
            continue;
        }

        if (BC[n].find("PML_") == 0) {
            int size = std::stoi(BC[n].substr(4));
            Set_BC_PML(n, size);
            continue;
        }

        throw std::runtime_error("Unknown boundary condition: " + BC[n]);
    }
}

// // Port Implementation
// void _openEMS::AddLumpedPort(int port_nr, double R, 
//                            const std::array<double, 3>& start,
//                            const std::array<double, 3>& stop,
//                            char p_dir, bool excite) {
//     if (!m_CSX)
//         throw std::runtime_error("AddLumpedPort: CSX is not set!");

//     auto port = CreateLumpedPort(port_nr, R, start, stop, p_dir, excite);

//     // Example: if an optional parameter "edges2grid" is available (here assumed as a string)
//     // that specifies which grid directions should receive extra grid lines. Replace
//     // this variable with the actual mechanism to pass extra options if needed.
//     std::string edges2grid = "";  // e.g., "xyz" or leave empty if not used

//     if (!edges2grid.empty()) {
//         auto grid = m_CSX->GetGrid();
//         for (char dir : edges2grid) {
//             int idx = 0;
//             switch (std::tolower(dir)) {
//                 case 'x': idx = 0; break;
//                 case 'y': idx = 1; break;
//                 case 'z': idx = 2; break;
//                 default: continue;
//             }
//             grid->AddLine(idx, start[idx]);
//             if (start[idx] != stop[idx])
//                 grid->AddLine(idx, stop[idx]);
//         }
//     }

//     // Optionally, store or further process "port" as needed.
// }

//""" CreateNF2FFBox(name='nf2ff', start=None, stop=None, **kw)
//
//Create a near - field to far - field box.
//
//This method will automatically adept the recording box to the current
//FDTD grid and boundary conditions.
//
//Notes
//---- -
//*Make sure the mesh grid and all boundary conditions are finially defined.
// 
//"""
_nf2ff _openEMS::CreateNF2FFBox(const std::string& name, const std::optional<std::array<double, 3>>& start, const std::optional<std::array<double, 3>>& stop)
{
    if (!m_CSX)
        throw std::runtime_error("CreateNF2FFBox: CSX is not set!");

    std::vector<bool> directions(6, true);
    std::vector<int> mirror(6, 0);
    std::vector<int> BC_size(6, 0);
    std::vector<int> BC_type(6, 0);

    for (int n = 0; n < 6; ++n)
    {
        BC_type[n] = Get_BC_Type(n);
        if (BC_type[n] == 0)
        {
            directions[n] = false;
            mirror[n] = 1; // PEC mirror
        }
        else if (BC_type[n] == 1)
        {
            directions[n] = false;
            mirror[n] = 2; // PMC mirror
        }
        else if (BC_type[n] == 2)
        {
            BC_size[n] = 2;
        }
        else if (BC_type[n] == 3)
        {
            BC_size[n] = Get_PML_Size(n) + 1;
        }
    }

    std::array<double, 3> actual_start{};
    std::array<double, 3> actual_stop{};

    if (!start.has_value() || !stop.has_value()) {
        _CSRectGrid grid = m_CSX->GetGrid();
        if (!grid.isValid())
            throw std::runtime_error("Error::CreateNF2FFBox: Grid is invalid");
        for (int n = 0; n < 3; ++n) {
            auto lines = grid.GetLines(n);
            if (lines.size() <= static_cast<size_t>(BC_size[2 * n] + BC_size[2 * n + 1]))
                throw std::runtime_error("Error::CreateNF2FFBox: not enough lines in some direction");
            actual_start[n] = lines[BC_size[2 * n]];
            actual_stop[n] = lines[lines.size() - BC_size[2 * n + 1] - 1];
        }
    } else {
        actual_start = start.value();
        actual_stop = stop.value();
    }

    return _nf2ff(m_CSX, name.c_str(), actual_start, actual_stop, directions, mirror);
}

// // Grid Management
// void _openEMS::AddEdges2Grid(const std::string& dirs,
//                            _CSX* csx,
//                            const std::vector<void*>& primitives) {
//     auto grid = csx->GetGrid();
    
//     for (auto prim_ptr : primitives) {
//         auto prim = static_cast<Primitive*>(prim_ptr);
//         auto edges = prim->GetEdges();

//         for (const auto& edge : edges) {
//             for (char dir : dirs) {
//                 switch (tolower(dir)) {
//                     case 'x':
//                         grid->AddLine(0, edge.start[0]);
//                         grid->AddLine(0, edge.end[0]);
//                         break;
//                     case 'y':
//                         grid->AddLine(1, edge.start[1]);
//                         grid->AddLine(1, edge.end[1]);
//                         break;
//                     case 'z':
//                         grid->AddLine(2, edge.start[2]);
//                         grid->AddLine(2, edge.end[2]);
//                         break;
//                 }
//             }
//         }
//     }
// }

int _openEMS::Run(const std::string &sim_path, bool cleanup, bool setup_only, const std::vector<std::string>& options)
{
    namespace fs = std::filesystem;
    fs::path simPath(sim_path);
    fs::path oldPath = fs::current_path();

    // Cleanup and create simulation directory, then change current directory.
    if (cleanup && fs::exists(simPath))
        fs::remove_all(simPath);
    if (!fs::exists(simPath))
        fs::create_directory(simPath);
    fs::current_path(simPath);

    // // Verify working directory
    // if (fs::current_path() != fs::absolute(simPath))
    //     throw std::runtime_error("Current working directory is different from sim_path");

    std::cout << "Current path " << fs::current_path() << std::endl;

    // Pass the options to the library (assumed available)
    SetLibraryArguments(options);

    // Display welcome screen.
    WelcomeScreen();

    // Setup the FDTD simulation.
    int EC = SetupFDTD();
    if (EC != 0)
        std::cout << "Run: Setup failed, error code: " << EC << std::endl;

    if (setup_only || EC != 0)
    {
        fs::current_path(oldPath);
        return EC;
    }

    // Run the FDTD simulation.
    RunFDTD();

    // Restore path
    fs::current_path(oldPath);
    return 0;
}
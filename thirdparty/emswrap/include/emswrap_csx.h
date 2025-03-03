#pragma once

#include <ContinuousStructure.h>
#include <CSRectGrid.h>
#include <CSPropMetal.h>
#include <CSPropMaterial.h>
#include <CSPropExcitation.h>
#include <CSPropConductingSheet.h>
#include <CSPropLumpedElement.h>
#include <CSPropProbeBox.h>
#include <CSPropDumpBox.h>

#include <emswrap_csrectgrid.h>
#include <emswrap_csproperties.h>
#include <emswrap_smoothmeshlines.h>

#include <vector>
#include <map>
#include <stdexcept>
#include <string>
#include <limits>
#include <emswrap_ports.h>

class _CSX {
public:
    _CSX(ContinuousStructure *structure) : _ptr(structure)
    {
    }

    virtual ~_CSX()
    {
    }

    // Grid management
    CSRectGrid* GetGrid() {
        return _ptr->GetGrid();
    }

    void SetMeshType(CoordinateSystem cs_type) {
        GetGrid()->SetMeshType(cs_type);
        _ptr->SetCoordInputType(cs_type);
    }

    void DefineGrid(const std::map<char, std::vector<double>>& mesh, double unit, double smooth_mesh_res = 0) {
        _CSRectGrid grid(GetGrid());
        grid.clear();
        for (const auto& entry : mesh) {
            char axis = entry.first;
            std::vector<double> lines = entry.second;
            if (smooth_mesh_res > 0) {
                lines = SmoothMesh::SmoothMeshLines(lines, smooth_mesh_res);
            }
            grid.SetLines(axis, lines);
        }
        grid.SetDeltaUnit(unit);
    }

    // Property creation
    CSPropMetal* AddMetal(const std::string& name) {
        CSPropMetal* prop = new CSPropMetal(_ptr->GetParameterSet());
        prop->SetName(name);
        _ptr->AddProperty(prop);
        return prop;
    }

    _CSPropMaterial AddMaterial(const std::string& name,
        double epsilon = 1.0,
        double mue = 1.0,
        double kappa = 0.0,
        double sigma = 0.0,
        double density = 0.0)
    {
        auto prop = new CSPropMaterial(_ptr->GetParameterSet());
        prop->SetName(name);
        prop->SetEpsilon(epsilon);
        prop->SetMue(mue);
        prop->SetKappa(kappa);
        prop->SetSigma(sigma);
        prop->SetDensity(density);

        _ptr->AddProperty(prop);
        return _CSPropMaterial(prop);
    }

    CSPropLumpedElement *AddLumpedElement(const std::string &name, int direction, bool caps = true)
    {
        CSPropLumpedElement *prop = new CSPropLumpedElement(_ptr->GetParameterSet());
        prop->SetName(name);
        prop->SetDirection(direction);
        prop->SetCaps(caps);
        _ptr->AddProperty(prop);
        return prop;
    }

    CSPropConductingSheet* AddConductingSheet(const std::string& name, double conductivity, double thickness) {
        CSPropConductingSheet* prop = new CSPropConductingSheet(_ptr->GetParameterSet());
        prop->SetName(name);
        prop->SetConductivity(conductivity);
        prop->SetThickness(thickness);
        _ptr->AddProperty(prop);
        return prop;
    }

    _CSPropExcitation AddExcitation(const std::string& name, int exc_type, std::array<double,3> exc_val) {
        CSPropExcitation* prop = new CSPropExcitation(_ptr->GetParameterSet());
        prop->SetName(name);
        prop->SetExcitType(exc_type);
        prop->SetExcitation(exc_val[0], 0);
        prop->SetExcitation(exc_val[1], 1);
        prop->SetExcitation(exc_val[2], 2);
        _ptr->AddProperty(prop);
        return _CSPropExcitation(prop);
    }

    _CSPropProbeBox AddProbe(const std::string& name, int p_type) {
        CSPropProbeBox* prop = new CSPropProbeBox(_ptr->GetParameterSet());
        prop->SetName(name);
        prop->SetProbeType(p_type);
        _ptr->AddProperty(prop);
        return _CSPropProbeBox(prop);
    }

    // """
    // Dump property to create field dumps.

    // Depending on the EM modeling tool there exist different dump types, dump
    // modes and file types with different meanings.

    // openEMS dump types:

    // * 0  : for E-field time-domain dump (default)
    // * 1  : for H-field time-domain dump
    // * 2  : for electric current time-domain dump
    // * 3  : for total current density (rot(H)) time-domain dump

    // * 10 : for E-field frequency-domain dump
    // * 11 : for H-field frequency-domain dump
    // * 12 : for electric current frequency-domain dump
    // * 13 : for total current density (rot(H)) frequency-domain dump

    // * 20 : local SAR frequency-domain dump
    // * 21 :  1g averaging SAR frequency-domain dump
    // * 22 : 10g averaging SAR frequency-domain dump

    // * 29 : raw data needed for SAR calculations (electric field FD, cell volume, conductivity and density)

    // openEMS dump modes:

    // * 0 : no-interpolation
    // * 1 : node-interpolation (default, see warning below)
    // * 2 : cell-interpolation (see warning below)

    // openEMS file types:

    // * 0 : vtk-file  (default)
    // * 1 : hdf5-file (easier to read by python, using h5py)

    // :param dump_type: dump type (see above)
    // :param dump_mode: dump mode (see above)
    // :param file_type: file type (see above)
    // :param frequency: specify a frequency vector (required for dump types >=10)
    // :param sub_sampling:   field domain sub-sampling, e.g. '2,2,4'
    // :param opt_resolution: field domain dump resolution, e.g. '10' or '10,20,5'

    // Notes
    // -----
    // openEMS FDTD Interpolation abnormalities:

    // * no-interpolation:
    //     fields are located in the mesh by the Yee-scheme, the mesh only
    //     specifies E- or H-Yee-nodes --> use node- or cell-interpolation or
    //     be aware of the offset
    // * E-field dump & node-interpolation:
    //     normal electric fields on boundaries will have false amplitude due to
    //     forward/backward interpolation  in case of (strong) changes in material
    //     permittivity or on metal surfaces
    //     --> use no- or cell-interpolation
    // * H-field dump & cell-interpolation:
    //     normal magnetic fields on boundaries will have false amplitude due to
    //     forward/backward interpolation in case of (strong) changes in material
    //     permeability --> use no- or node-interpolation
    // """

    CSPropDumpBox* AddDump(const std::string &name, int dump_type, int dump_mode, int file_type, std::vector<double> frequency = {}) {
        CSPropDumpBox* prop = new CSPropDumpBox(_ptr->GetParameterSet());
        prop->SetName(name);
        prop->SetDumpType(dump_type);
        prop->SetDumpMode(dump_mode);
        prop->SetFileType(file_type);

        if (dump_type >= 10 && frequency.empty())
            throw std::runtime_error("Frequency vector required for dump types >= 10");

        for (auto freq : frequency) {
            prop->AddFDSample(freq);
        }

        _ptr->AddProperty(prop);
        return prop;
    }

    // Port addition functions (MatLab-like API)
    
    // [CSX,port] = AddWaveGuidePort( CSX, prio, portnr, start, stop, dir, E_WG_func, H_WG_func, kc, exc_amp, varargin )
    WaveguidePort* AddWaveGuidePort( int priority, int port_nr,
                                     const std::array<double, 3>& start,
                                     const std::array<double, 3>& stop,
                                     int dir,
                                     const std::array<std::string, 3>& E_WG_func,
                                     const std::array<std::string, 3>& H_WG_func,
                                     double kc,
                                     double excite,
                                     const std::string& PortNamePrefix = "",
                                     double delay = 0.0)
    {
        // Create a WaveguidePort using the provided parameters.
        WaveguidePort* port = new WaveguidePort(_ptr,
                                                port_nr,
                                                start,
                                                stop,
                                                dir,
                                                E_WG_func,
                                                H_WG_func,
                                                kc,
                                                excite,
                                                PortNamePrefix,
                                                delay,
                                                priority);
        return port;
    }

    // [CSX,port] = AddMSLPort( CSX, prio, portnr, materialname, start, stop, dir, evec, varargin )
    MSLPort AddMSLPort(int priority,
                       int port_nr,
                       CSPropMetal *metal_prop,
                       const std::array<double, 3> &start,
                       const std::array<double, 3> &stop,
                       int dir,
                       const std::array<double, 3> &evec,
                       bool excite_port,
                       double FeedShift = 0.0,
                       double MeasPlaneShift = 0.0,
                       double Feed_R = std::numeric_limits<double>::infinity(),
                       const std::string &PortNamePrefix = "",
                       double delay = 0.0)
    {
        // Create an MSLPort using the provided parameters.
        return MSLPort(_ptr,
                       priority,
                       port_nr,
                       metal_prop,
                       start,
                       stop,
                       dir,
                       evec,
                       excite_port,
                       FeedShift,
                       MeasPlaneShift,
                       Feed_R,
                       PortNamePrefix,
                       delay);
    }

    // LumpedPort* AddLumpedPort( int port_nr,
    //                            double R,
    //                            const std::array<double, 3>& start,
    //                            const std::array<double, 3>& stop,
    //                            int exc_dir,
    //                            double excite,
    //                            const std::string &PortNamePrefix = "",
    //                            double delay = 0.0,
    //                            int priority = 0)
    // {
    //     // Create a LumpedPort using the provided parameters.
    //     LumpedPort* port = new LumpedPort(_ptr,
    //                                       port_nr,
    //                                       R,
    //                                       start,
    //                                       stop,
    //                                       exc_dir,
    //                                       excite,
    //                                       PortNamePrefix,
    //                                       delay,
    //                                       priority);
    //     return port;
    // }

    // RectWGPort* AddRectWGPort( int port_nr,
    //                            const std::array<double, 3>& start,
    //                            const std::array<double, 3>& stop,
    //                            int exc_dir,
    //                            double a,
    //                            double b,
    //                            const std::string& mode_name,
    //                            double excite = 0.0,
    //                            const std::string& PortNamePrefix = "",
    //                            double delay = 0.0,
    //                            int priority = 0)
    // {
    //     // Create a rectangular waveguide port with the given dimensions and mode.
    //     RectWGPort* port = new RectWGPort(_ptr,
    //                                       port_nr,
    //                                       start,
    //                                       stop,
    //                                       exc_dir,
    //                                       a,
    //                                       b,
    //                                       mode_name,
    //                                       excite,
    //                                       PortNamePrefix,
    //                                       delay,
    //                                       priority);
    //     return port;
    // }

    // Property retrieval
    std::vector<CSProperties*> GetAllProperties() {
        std::vector<CSProperties*> props;
        int qty = _ptr->GetQtyProperties();
        for (int i = 0; i < qty; ++i) {
            props.push_back(_ptr->GetProperty(i));
        }
        return props;
    }

    std::vector<CSProperties*> GetPropertiesByName(const std::string& name) {
        std::vector<CSProperties*> result;
        auto props = GetAllProperties();
        for (auto prop : props) {
            if (prop->GetName() == name) {
                result.push_back(prop);
            }
        }
        return result;
    }

    // XML operations
    void Write2XML(const std::string& filename) {
        _ptr->Write2XML(filename.c_str());
    }

    bool ReadFromXML(const std::string& filename) {
        return _ptr->ReadFromXML(filename.c_str());
    }

    ContinuousStructure* operator->()
    {
        return _ptr;
    }

    operator ContinuousStructure* ()
    {
        return _ptr;
    }

private:
    ContinuousStructure *_ptr;
    bool owns_ptr = false;
};

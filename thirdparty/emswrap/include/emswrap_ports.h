#pragma once

#include <string>
#include <vector>
#include <complex>
#include <stdexcept>
#include <cmath>

// Forward declarations (assumed to be defined elsewhere)
class ContinuousStructure;
class CSPropExcitation;


//
// UI_data: Holds time and frequency domain probe data.
//
class UI_data {
public:
    // Constructor: fns can be one or multiple filenames.
    UI_data(const std::vector<std::string> &fns,
            const std::string &path,
            const std::vector<double> &freq,
            const std::string &signal_type = "pulse");

    // UI_data(const std::string &fn,
    //         const std::string &path,
    //         const std::vector<double> &freq,
    //         const std::string &signal_type = "pulse");

    // UI_data(const std::string &fn,
    //         const std::string &path,
    //         double freq,
    //         const std::string &signal_type = "pulse");

    // Data arrays
    std::vector< std::vector<double> > ui_time;
    std::vector< std::vector<double> > ui_val;
    std::vector< std::vector<std::complex<double>> > ui_f_val;

private:
    // Function to perform DFT from time to frequency domain.
    std::vector<std::complex<double>> DFT_time2freq(const std::vector<double>& t,
                                                    const std::vector<double>& u,
                                                    const std::vector<double>& freq,
                                                    const std::string& signal_type);
};

//
// Base class for Ports
//
class Port {
public:
    Port(ContinuousStructure* csx,
         int port_nr,
         const std::vector<double>& start,
         const std::vector<double>& stop,
         double excite,
         const std::string& PortNamePrefix = "",
         double delay = 0.0,
         int priority = 0);
    virtual ~Port() {}

    // Enable or disable port excitation.
    // virtual void SetEnabled(bool val);

    // Read UI Data from simulation files.
    virtual void ReadUIData(const std::string& sim_path,
                            const std::vector<double>& freq,
                            const std::string& signal_type = "pulse");

    // Calculate port parameters.
    virtual void CalcPort(const std::string& sim_path,
                          const std::vector<double>& freq,
                          double ref_impedance = -1,
                          double ref_plane_shift = 0.0,
                          const std::string& signal_type = "pulse");

    // Member variables
    ContinuousStructure* CSX_ptr;
    int number;
    double excite;
    std::vector<double> start;
    std::vector<double> stop;
    double Z_ref;
    double delay;
    double measplane_shift;
    int priority;
    std::string lbl_temp;
    // Filenames for voltage and current probes.
    std::vector<std::string> U_filenames;
    std::vector<std::string> I_filenames;

    // Accumulated signals in time and frequency domain.
    std::vector<std::complex<double>> uf_tot;
    std::vector<double> ut_tot;
    std::vector<std::complex<double>> if_tot;
    std::vector<double> it_tot;

    // Incident signals in time and frequency domain.
    std::vector<std::complex<double>> uf_inc;
    std::vector<double> ut_inc;
    std::vector<std::complex<double>> if_inc;
    std::vector<double> it_inc;

    // Reflected signals in time and frequency domain.
    std::vector<std::complex<double>> uf_ref;
    std::vector<double> ut_ref;
    std::vector<std::complex<double>> if_ref;
    std::vector<double> it_ref;

    // Power parameters.
    std::vector<double> P_inc;
    std::vector<double> P_ref;
    std::vector<double> P_acc;

    std::vector<double> beta;

    // UI_data objects for U and I.
    UI_data* u_data;
    UI_data* i_data;

protected:
    // A container for properties (excitations, probes, etc.)
    // In a full implementation, this could be a vector of pointers to ContinuousStructureCAD primitives.
    std::vector<void*> port_props;
};

// //
// // LumpedPort: Derived from Port for lumped port implementation.
// //
// class LumpedPort : public Port {
// public:
//     LumpedPort(ContinuousStructure* csx,
//                int port_nr,
//                double R,
//                const std::vector<double>& start,
//                const std::vector<double>& stop,
//                int exc_dir,
//                double excite = 0.0,
//                const std::string& PortNamePrefix = "",
//                double delay = 0.0,
//                int priority = 0);
//     virtual void CalcPort(const std::string& sim_path,
//                           const std::vector<double>& freq,
//                           double ref_impedance = -1,
//                           double ref_plane_shift = 0.0,
//                           const std::string& signal_type = "pulse") override;

//     double R;      // Port resistor value
//     int exc_ny;    // Excitation direction index
// };

// //
// // MSLPort: Microstrip transmission line port.
// //
// class MSLPort : public Port {
// public:
//     MSLPort(ContinuousStructure* csx,
//             int port_nr,
//             void* metal_prop, // pointer to metal properties object
//             const std::vector<double>& start,
//             const std::vector<double>& stop,
//             int prop_dir,
//             int exc_dir,
//             double excite = 0.0,
//             double FeedShift = 0.0,
//             double MeasPlaneShift = 0.0,
//             double Feed_R = std::numeric_limits<double>::infinity(),
//             const std::string& PortNamePrefix = "",
//             double delay = 0.0,
//             int priority = 0);
    
//     // Overridden method to read UI data for MSLPort.
//     virtual void ReadUIData(const std::string& sim_path,
//                             const std::vector<double>& freq,
//                             const std::string& signal_type = "pulse") override;

//     // Additional members for MSLPort
//     int exc_ny;
//     int prop_ny;
//     int direction;
//     int upside_down;
//     double feed_shift;
//     double measplane_shift;
//     double measplane_pos;
//     double feed_R;
//     // For voltage/current probe positions
//     std::vector<double> U_delta;
//     std::vector<double> I_delta;
// };

// //
// // WaveguidePort: Base class for waveguide ports
// //
// class WaveguidePort : public Port {
// public:
//     WaveguidePort(ContinuousStructure* csx,
//                   int port_nr,
//                   const std::vector<double>& start,
//                   const std::vector<double>& stop,
//                   int exc_dir,
//                   const std::vector<std::string>& E_WG_func,
//                   const std::vector<std::string>& H_WG_func,
//                   double kc,
//                   double excite = 0.0,
//                   const std::string& PortNamePrefix = "",
//                   double delay = 0.0,
//                   int priority = 0);
    
//     virtual void CalcPort(const std::string& sim_path,
//                           const std::vector<double>& freq,
//                           double ref_impedance = -1,
//                           double ref_plane_shift = 0.0,
//                           const std::string& signal_type = "pulse") override;

//     // Additional members for waveguide ports.
//     int exc_ny;
//     int ny_P;
//     int ny_PP;
//     double kc;     // cutoff parameter
//     std::vector<std::string> E_func;
//     std::vector<std::string> H_func;
//     double ref_index;
//     std::vector<double> ZL;
// };

// //
// // RectWGPort: Rectangular waveguide port using TE modes.
// //
// class RectWGPort : public WaveguidePort {
// public:
//     RectWGPort(ContinuousStructure* csx,
//                int port_nr,
//                const std::vector<double>& start,
//                const std::vector<double>& stop,
//                int exc_dir,
//                double a,
//                double b,
//                const std::string& mode_name,
//                double excite = 0.0,
//                const std::string& PortNamePrefix = "",
//                double delay = 0.0,
//                int priority = 0);

//     // WG_size stores a and b.
//     std::vector<double> WG_size;
//     std::string WG_mode;
//     bool TE;
//     bool TM;
//     double M;
//     double N;
//     double unit;
// };
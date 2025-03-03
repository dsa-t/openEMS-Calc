#pragma once

#include <string>
#include <vector>
#include <array>
#include <complex>
#include <stdexcept>
#include <cmath>
#include <optional>

// Forward declarations (assumed to be defined elsewhere)
class ContinuousStructure;
class CSPropExcitation;
class CSProperties;


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
    Port(ContinuousStructure *csx,
         int port_nr,
         const std::array<double, 3> &start,
         const std::array<double, 3> &stop,
         double excite,
         const std::string &PortNamePrefix = "",
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
        std::optional<std::complex<double>> ref_impedance = std::nullopt,
        double ref_plane_shift = 0.0,
        const std::string& signal_type = "pulse");

    // Member variables
    ContinuousStructure* CSX_ptr;
    int number;
    double excite;
    std::array<double, 3> start;
    std::array<double, 3> stop;
    std::complex<double> Z_ref;
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

    std::vector<std::complex<double>> beta;

    // UI_data objects for U and I.
    std::unique_ptr<UI_data> u_data;
    std::unique_ptr<UI_data> i_data;

protected:
    // A container for properties (excitations, probes, etc.)
    std::vector<CSProperties*> port_props;
};

//
// LumpedPort: A port implementation based on a lumped element
//
class LumpedPort : public Port {
public:
    LumpedPort(ContinuousStructure* csx,
               int port_nr,
               double R,
               const std::array<double, 3>& start,
               const std::array<double, 3>& stop,
               int exc_dir,
               double excite,
               const std::string& PortNamePrefix = "",
               double delay = 0.0,
               int priority = 0);

    virtual void CalcPort(const std::string& sim_path,
        const std::vector<double>& freq,
        std::optional<std::complex<double>> ref_impedance = std::nullopt,
        double ref_plane_shift = 0.0,
        const std::string& signal_type = "pulse") override;

private:
    double R;
    int exc_ny;
    int direction;
};

//
// MSLPort: Microstrip transmission line port.
//
class MSLPort : public Port {
public:
    MSLPort(ContinuousStructure* csx,
            int port_nr,
            CSProperties *metal_prop,
            const std::array<double, 3>& start,
            const std::array<double, 3>& stop,
            int prop_dir,
            int exc_dir,
            double excite = 0.0,
            double FeedShift = 0.0,
            double MeasPlaneShift = std::numeric_limits<double>::signaling_NaN(),
            double Feed_R = std::numeric_limits<double>::infinity(),
            const std::string& PortNamePrefix = "",
            double delay = 0.0,
            int priority = 0);

    MSLPort(ContinuousStructure *cs_ptr, 
        int prio, 
        int port_nr, 
        CSProperties* metal_prop, 
        const std::array<double, 3> &start, 
        const std::array<double, 3> &stop, 
        int dir, const std::array<double, 3> &evec, 
        bool excite_port, 
        double FeedShift = 0.0,
        double MeasPlaneShift = std::numeric_limits<double>::signaling_NaN(),
        double Feed_R = std::numeric_limits<double>::infinity(),
        const std::string& PortNamePrefix = "",
        double delay = 0.0);

    // Overridden method to read UI data for MSLPort.
    virtual void ReadUIData(const std::string& sim_path,
                            const std::vector<double>& freq,
                            const std::string& signal_type = "pulse") override;

    // Additional members for MSLPort
    int exc_ny;
    int prop_ny;
    int direction;
    int upside_down;
    double feed_shift;
    double measplane_shift;
    double measplane_pos;
    double feed_R;
    // For voltage/current probe positions
    std::vector<double> U_delta;
    std::vector<double> I_delta;
};

//
// WaveguidePort: Base class for waveguide ports
//
class WaveguidePort : public Port
{
public:
    WaveguidePort(ContinuousStructure *csx,
                  int port_nr,
                  const std::array<double, 3> &start,
                  const std::array<double, 3> &stop,
                  int exc_dir,
                  const std::array<std::string, 3>& E_WG_func,
                  const std::array<std::string, 3>& H_WG_func,
                  double kc,
                  double excite,
                  const std::string &PortNamePrefix = "",
                  double delay = 0.0,
                  int priority = 0);

    virtual void CalcPort(const std::string& sim_path,
        const std::vector<double>& freq,
        std::optional<std::complex<double>> ref_impedance = std::nullopt,
        double ref_plane_shift = 0.0,
        const std::string& signal_type = "pulse") override;

    // Additional members for waveguide ports.
    int exc_ny;
    int ny_P;
    int ny_PP;
    double kc; // cutoff parameter
    std::array<std::string, 3> E_func;
    std::array<std::string, 3> H_func;
    double ref_index;
    std::vector<std::complex<double>> ZL;
};

class RectWGPort : public WaveguidePort {
    public:
        RectWGPort(ContinuousStructure *csx,
                   int port_nr,
                   const std::array<double, 3> &start,
                   const std::array<double, 3> &stop,
                   int exc_dir,
                   double a,
                   double b,
                   const std::string &mode_name,
                   double excite = 0.0,
                   const std::string &PortNamePrefix = "",
                   double delay = 0.0,
                   int priority = 0);
    };
    

#include "emswrap_ports.h"
#include <emswrap_compat.h>

#include <ContinuousStructure.h>

#include <iostream>
#include <filesystem>
#include <fstream>
#include <complex>
#include <cmath>
#include <vector>
#include <string>

// --- UI_data Implementation ---

// UI_data::UI_data(const std::string &fn,
//                  const std::string &path,
//                  const std::vector<double> &freq,
//                  const std::string &signal_type) : UI_data({fn}, path, freq, signal_type)
// {
// }

// UI_data::UI_data(const std::string &fn,
//                  const std::string &path,
//                  double freq,
//                  const std::string &signal_type) : UI_data({fn}, path, {freq}, signal_type)
// {
// }

UI_data::UI_data(const std::vector<std::string>& fns,
                 const std::string& path,
                 const std::vector<double>& freq,
                 const std::string& signal_type)
    : ui_time(), ui_val(), ui_f_val() {
    // For each filename, read the time-domain data and compute its DFT.
    std::cout << "UI_data constructed with " << fns.size() << " file(s) at path: " << path << std::endl;
    for (const auto& fn : fns) {
        // Construct full file path (assumes '/' as path separator)
        std::string full_path = path + "/" + fn;
        auto abs_path = std::filesystem::absolute(full_path);

        std::cout << "[UI_data] reading" << abs_path.u8string() << std::endl;

        std::ifstream infile(abs_path);
        if (!infile.is_open()) {
            std::cerr << "Error opening file: " << full_path << std::endl;
            
            // Get the error code
            std::ios_base::iostate state = infile.rdstate();

            // Check for specific error bits
            if (state & std::ios_base::eofbit)
            {
                std::cout << "End of file reached." << std::endl;
            }
            if (state & std::ios_base::failbit)
            {
                std::cout << "Non-fatal I/O error occurred." << std::endl;
            }
            if (state & std::ios_base::badbit)
            {
                std::cout << "Fatal I/O error occurred." << std::endl;
            }

            // Print system error message
            std::perror("Error: ");

            throw std::runtime_error("Error opening file");
        }
        std::vector<double> times;
        std::vector<double> values;
        std::string line;
        // Read file line-by-line, skipping comments starting with '%'
        while (std::getline(infile, line)) {
            if (line.empty() || line[0] == '%') continue;
            std::istringstream iss(line);
            double t_val, u_val;
            if (!(iss >> t_val >> u_val)) continue;
            times.push_back(t_val);
            values.push_back(u_val);
        }
        infile.close();
        
        // Store the time and value data for this file.
        ui_time.push_back(times);
        ui_val.push_back(values);
        
        // Compute the DFT for these signals using the helper function.
        std::vector<std::complex<double>> dft = DFT_time2freq(times, values, freq, signal_type);
        ui_f_val.push_back(dft);
    }
}

std::vector<std::complex<double>> UI_data::DFT_time2freq(const std::vector<double>& t,
                                                          const std::vector<double>& u,
                                                          const std::vector<double>& freq,
                                                          const std::string& signal_type) {
    if (t.size() != u.size()) {
        throw std::invalid_argument("t and u must have the same length.");
    }
    if (freq.empty()) {
        throw std::invalid_argument("freq must not be empty.");
    }

    std::vector<std::complex<double>> result(freq.size(), std::complex<double>(0.0, 0.0));
    const double two_pi = 2.0 * M_PI;
    // Compute DFT for each frequency.
    for (size_t k = 0; k < freq.size(); ++k) {
        std::complex<double> sum(0.0, 0.0);
        for (size_t n = 0; n < t.size(); ++n) {
            double phase = -two_pi * freq[k] * t[n];
            sum += u[n] * std::exp(std::complex<double>(0.0, phase));
        }
        result[k] = sum;
    }

    // Apply scaling based on signal type.
    if (signal_type == "pulse") {
        double dt = (t.size() > 1) ? (t[1] - t[0]) : 1.0;
        for (auto& val : result) {
            val *= dt;
        }
    } else if (signal_type == "periodic") {
        for (auto& val : result) {
            val /= static_cast<double>(t.size());
        }
    } else {
        throw std::runtime_error("Unknown signal type: \"" + signal_type + "\"");
    }

    // Return single-sided spectrum.
    for (auto& val : result) {
        val *= 2;
    }
    return result;
}


// --- Port Implementation ---

Port::Port(ContinuousStructure* csx,
           int port_nr,
           const std::vector<double>& start,
           const std::vector<double>& stop,
           double excite,
           const std::string& PortNamePrefix,
           double delay,
           int priority)
    : CSX_ptr(csx),
      number(port_nr),
      excite(excite),
      start(start),
      stop(stop),
      Z_ref(-1),
      delay(delay),
      priority(priority),
      lbl_temp(PortNamePrefix + "port_" + std::to_string(port_nr)),
      u_data(nullptr),
      i_data(nullptr),
      measplane_shift(0.0)
{
    // Filenames (for probes) can be set by derived classes.
    std::cout << "Constructed Port " << number << std::endl;
}

void Port::ReadUIData(const std::string& sim_path,
                        const std::vector<double>& freq,
                        const std::string& signal_type) {
    // Instantiate UI_data for voltage probes.
    if(u_data) { delete u_data; u_data = nullptr; }
    u_data = new UI_data(U_filenames, sim_path, freq, signal_type);
    // Initialize accumulator vectors to zero.
    uf_tot.assign(freq.size(), std::complex<double>(0.0, 0.0));
    ut_tot.assign(freq.size(), 0.0);
    // Sum up frequency-domain and time-domain voltage data from each file.
    for (const auto& data : u_data->ui_f_val) {
        for (size_t i = 0; i < freq.size(); ++i)
            uf_tot[i] += data[i];
    }
    for (const auto& data : u_data->ui_val) {
        for (size_t i = 0; i < freq.size(); ++i)
            ut_tot[i] += data[i];
    }

    // Instantiate UI_data for current probes.
    if(i_data) { delete i_data; i_data = nullptr; }
    i_data = new UI_data(I_filenames, sim_path, freq, signal_type);
    if_tot.assign(freq.size(), std::complex<double>(0.0, 0.0));
    it_tot.assign(freq.size(), 0.0);
    // Sum up frequency-domain and time-domain current data from each file.
    for (const auto& data : i_data->ui_f_val) {
        for (size_t i = 0; i < freq.size(); ++i)
            if_tot[i] += data[i];
    }
    for (const auto& data : i_data->ui_val) {
        for (size_t i = 0; i < freq.size(); ++i)
            it_tot[i] += data[i];
    }
}

void Port::CalcPort(const std::string& sim_path,
                    const std::vector<double>& freq,
                    double ref_impedance,
                    double ref_plane_shift,
                    const std::string& signal_type) {
    // Read UI-data (voltage and current).
    ReadUIData(sim_path, freq, signal_type);

    // Set the port reference impedance.
    if(ref_impedance > 0)
        Z_ref = ref_impedance;
    if(Z_ref < 0)
        throw std::runtime_error("Port Z_ref should not be negative!");

    // Apply reference plane shift if provided.
    if(ref_plane_shift != 0) {
        // In the Python version, the port must have a beta attribute.
        // Here we assume 'beta' has been set by a derived class.
        if(beta.empty()) {
            throw std::runtime_error("Port has no beta attribute!");
        }
        double shift = ref_plane_shift;
        // If a measurement plane shift was previously set, subtract it.
        if(measplane_shift != 0)
            shift -= measplane_shift;
        // Obtain delta unit from the grid (assumed to be provided by CSX_ptr).
        double delta_unit = CSX_ptr->GetGrid()->GetDeltaUnit();
        shift *= delta_unit;
        // Calculate phase shift using the propagation constant (beta).
        double phase = beta[0] * shift;  // using the first frequency point
        std::vector<std::complex<double>> uf_new(freq.size());
        std::vector<std::complex<double>> if_new(freq.size());
        for (size_t i = 0; i < freq.size(); ++i) {
            uf_new[i] = uf_tot[i] * std::cos(-phase)
                        + std::complex<double>(0,1) * if_tot[i] * Z_ref * std::sin(-phase);
            if_new[i] = if_tot[i] * std::cos(-phase)
                        + std::complex<double>(0,1) * uf_tot[i] / Z_ref * std::sin(-phase);
        }
        uf_tot = uf_new;
        if_tot = if_new;
    }

    // Calculate incident voltage and current.
    uf_inc.resize(freq.size());
    if_inc.resize(freq.size());
    for (size_t i = 0; i < freq.size(); ++i) {
        uf_inc[i] = 0.5 * (uf_tot[i] + if_tot[i] * Z_ref);
        if_inc[i] = 0.5 * (if_tot[i] + uf_tot[i] / Z_ref);
    }

    // Calculate reflected voltage and current.
    uf_ref.resize(freq.size());
    if_ref.resize(freq.size());
    for (size_t i = 0; i < freq.size(); ++i) {
        uf_ref[i] = uf_tot[i] - uf_inc[i];
        if_ref[i] = if_inc[i] - if_tot[i];
    }

    // For time-domain signals, compute incident and reflected signals.
    ut_inc.resize(freq.size());
    it_inc.resize(freq.size());
    ut_ref.resize(freq.size());
    it_ref.resize(freq.size());
    for (size_t i = 0; i < freq.size(); ++i) {
        ut_inc[i] = 0.5 * (ut_tot[i] + it_tot[i] * Z_ref);
        it_inc[i] = 0.5 * (it_tot[i] + ut_tot[i] / Z_ref);
        ut_ref[i] = ut_tot[i] - ut_inc[i];
        it_ref[i] = it_tot[i] - it_inc[i];
    }

    // Calculate port power parameters.
    P_inc.resize(freq.size());
    P_ref.resize(freq.size());
    P_acc.resize(freq.size());
    for (size_t i = 0; i < freq.size(); ++i) {
        P_inc[i] = 0.5 * std::real(uf_inc[i] * std::conj(if_inc[i]));
        P_ref[i] = 0.5 * std::real(uf_ref[i] * std::conj(if_ref[i]));
        P_acc[i] = 0.5 * std::real(uf_tot[i] * std::conj(if_tot[i]));
    }
    std::cout << "CalcPort for Port " << number << " using Z_ref: " << Z_ref << std::endl;
}


// // --- LumpedPort Implementation ---

// LumpedPort::LumpedPort(ContinuousStructure* csx,
//                        int port_nr,
//                        double R,
//                        const std::vector<double>& start,
//                        const std::vector<double>& stop,
//                        int exc_dir,
//                        double excite,
//                        const std::string& PortNamePrefix,
//                        double delay,
//                        int priority)
//     : Port(csx, port_nr, start, stop, excite, PortNamePrefix, delay, priority),
//       R(R),
//       exc_ny(exc_dir)
// {
//     std::cout << "Constructed LumpedPort " << number << " with resistor R = " << R << std::endl;
// }

// void LumpedPort::CalcPort(const std::string& sim_path,
//                             const std::vector<double>& freq,
//                             double ref_impedance,
//                             double ref_plane_shift,
//                             const std::string& signal_type) {
//     std::cout << "CalcPort for LumpedPort " << number << std::endl;
//     // Call base class calculation.
//     Port::CalcPort(sim_path, freq, ref_impedance, ref_plane_shift, signal_type);
// }


// // --- MSLPort Implementation ---

// MSLPort::MSLPort(ContinuousStructure* csx,
//                  int port_nr,
//                  void* metal_prop,
//                  const std::vector<double>& start,
//                  const std::vector<double>& stop,
//                  int prop_dir,
//                  int exc_dir,
//                  double excite,
//                  double FeedShift,
//                  double MeasPlaneShift,
//                  double Feed_R,
//                  const std::string& PortNamePrefix,
//                  double delay,
//                  int priority)
//     : Port(csx, port_nr, start, stop, excite, PortNamePrefix, delay, priority),
//       exc_ny(exc_dir),
//       prop_ny(prop_dir),
//       direction((stop[prop_dir] - start[prop_dir]) >= 0 ? 1 : -1),
//       upside_down(0),
//       feed_shift(FeedShift),
//       measplane_shift(MeasPlaneShift),
//       measplane_pos(start[prop_dir] + MeasPlaneShift),
//       feed_R(Feed_R)
// {
//     std::cout << "Constructed MSLPort " << number << std::endl;
// }

// void MSLPort::ReadUIData(const std::string& sim_path,
//                            const std::vector<double>& freq,
//                            const std::string& signal_type) {
//     std::cout << "MSLPort " << number << " reading UI data." << std::endl;
//     // Create UI_data objects. In an actual implementation, U_delta and I_delta would be used.
//     u_data = new UI_data(U_filenames, sim_path, freq, signal_type);
//     // Dummy: assign the second probe voltage if available
//     if (!u_data->ui_f_val.empty() && u_data->ui_f_val.size() > 1) {
//         uf_tot = u_data->ui_f_val[1];
//     } else {
//         uf_tot = std::vector<std::complex<double>>(freq.size(), std::complex<double>(0.0, 0.0));
//     }
//     ut_tot = std::vector<double>(freq.size(), 0.0);

//     i_data = new UI_data(I_filenames, sim_path, freq, signal_type);
//     if (!i_data->ui_f_val.empty() && i_data->ui_f_val.size() > 1) {
//         if_tot = i_data->ui_f_val[1];
//     } else {
//         if_tot = std::vector<std::complex<double>>(freq.size(), std::complex<double>(0.0, 0.0));
//     }
//     it_tot = std::vector<double>(freq.size(), 0.0);
// }


// // --- WaveguidePort Implementation ---

// WaveguidePort::WaveguidePort(ContinuousStructure* csx,
//                              int port_nr,
//                              const std::vector<double>& start,
//                              const std::vector<double>& stop,
//                              int exc_dir,
//                              const std::vector<std::string>& E_WG_func,
//                              const std::vector<std::string>& H_WG_func,
//                              double kc,
//                              double excite,
//                              const std::string& PortNamePrefix,
//                              double delay,
//                              int priority)
//     : Port(csx, port_nr, start, stop, excite, PortNamePrefix, delay, priority),
//       exc_ny(exc_dir),
//       kc(kc),
//       E_func(E_WG_func),
//       H_func(H_WG_func),
//       ref_index(1.0)
// {
//     // Initialize any additional members.
//     ny_P = 0;
//     ny_PP = 0;
//     std::cout << "Constructed WaveguidePort " << number << std::endl;
// }

// void WaveguidePort::CalcPort(const std::string& sim_path,
//                              const std::vector<double>& freq,
//                              double ref_impedance,
//                              double ref_plane_shift,
//                              const std::string& signal_type) {
//     std::cout << "CalcPort for WaveguidePort " << number << std::endl;
//     // Call base class calculation.
//     Port::CalcPort(sim_path, freq, ref_impedance, ref_plane_shift, signal_type);
//     // Dummy assignment of propagation constant and load impedance.
//     beta = std::vector<double>(freq.size(), 1.0);
//     ZL = std::vector<double>(freq.size(), 50.0);
// }


// // --- RectWGPort Implementation ---

// RectWGPort::RectWGPort(ContinuousStructure* csx,
//                        int port_nr,
//                        const std::vector<double>& start,
//                        const std::vector<double>& stop,
//                        int exc_dir,
//                        double a,
//                        double b,
//                        const std::string& mode_name,
//                        double excite,
//                        const std::string& PortNamePrefix,
//                        double delay,
//                        int priority)
//     : WaveguidePort(csx, port_nr, start, stop, exc_dir, std::vector<std::string>{}, std::vector<std::string>{}, 0.0, excite, PortNamePrefix, delay, priority),
//       WG_mode(mode_name),
//       TE(false),
//       TM(false),
//       M(0),
//       N(0),
//       unit(1.0)
// {
//     // Store the rectangular dimensions.
//     WG_size.push_back(a);
//     WG_size.push_back(b);
//     // Set mode flags based on the mode name.
//     if (mode_name.find("TE") != std::string::npos) {
//         TE = true;
//     }
//     else if (mode_name.find("TM") != std::string::npos) {
//         TM = true;
//     }
//     std::cout << "Constructed RectWGPort " << number << " with mode " << WG_mode << std::endl;
// }

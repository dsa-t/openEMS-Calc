
#include "emswrap_ports.h"
#include <emswrap_compat.h>
#include <emswrap_csproperties.h>
#include <emswrap_csrectgrid.h>
#include <emswrap_csx.h>

#include <ContinuousStructure.h>

#include <iostream>
#include <filesystem>
#include <fstream>
#include <complex>
#include <cmath>
#include <vector>
#include <string>
#include <numeric>

// --- UI_data Implementation ---

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

        std::cout << "[UI_data] reading" << abs_path.string() << std::endl;

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

Port::Port(ContinuousStructure *csx,
           int port_nr,
           const std::array<double, 3> &start,
           const std::array<double, 3> &stop,
           double excite,
           const std::string &PortNamePrefix,
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
    u_data = std::make_unique<UI_data>(U_filenames, sim_path, freq, signal_type);
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
    i_data = std::make_unique<UI_data>(I_filenames, sim_path, freq, signal_type);
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
    std::optional<std::complex<double>> ref_impedance,
    double ref_plane_shift,
    const std::string& signal_type) {

    // Read UI-data (voltage and current).
    ReadUIData(sim_path, freq, signal_type);

    // Set the port reference impedance.
    if (ref_impedance)
        Z_ref = *ref_impedance;

    // Apply reference plane shift if provided.
    if (ref_plane_shift != 0) {
        if (beta.empty()) {
            throw std::runtime_error("Port has no beta attribute!");
        }
        double shift = ref_plane_shift;
        if (measplane_shift != 0)
            shift -= measplane_shift;
        double delta_unit = CSX_ptr->GetGrid()->GetDeltaUnit();
        shift *= delta_unit;
        // Calculate phase using the real part of beta[0]
        double phase = std::real(beta[0]) * shift;
        std::vector<std::complex<double>> uf_new(freq.size());
        std::vector<std::complex<double>> if_new(freq.size());
        // Use temporary vectors to update uf_tot and if_tot (do not mix updated values)
        for (size_t i = 0; i < freq.size(); ++i) {
            uf_new[i] = uf_tot[i] * std::cos(-phase)
                + std::complex<double>(0, 1) * if_tot[i] * Z_ref * std::sin(-phase);
            if_new[i] = if_tot[i] * std::cos(-phase)
                + std::complex<double>(0, 1) * uf_tot[i] / Z_ref * std::sin(-phase);
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
    size_t time_size = ut_tot.size();
    ut_inc.resize(time_size);
    it_inc.resize(time_size);
    ut_ref.resize(time_size);
    it_ref.resize(time_size);
    for (size_t i = 0; i < time_size; ++i) {
        auto val0 = 0.5 * (ut_tot[i] + it_tot[i] * Z_ref);
        auto val1 = 0.5 * (it_tot[i] + ut_tot[i] / Z_ref);

        if (val0.imag() != 0.0)
            throw std::runtime_error("Port::CalcPort val0 has imag part");

        if (val1.imag() != 0.0)
            throw std::runtime_error("Port::CalcPort val1 has imag part");

        ut_inc[i] = val0.real();
        it_inc[i] = val1.real();
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


// --- LumpedPort Implementation ---

LumpedPort::LumpedPort(ContinuousStructure* cs_ptr,
                       int port_nr,
                       double R,
                       const std::array<double, 3>& start,
                       const std::array<double, 3>& stop,
                       int exc_dir,
                       double excite,
                       const std::string& PortNamePrefix,
                       double delay,
                       int priority)
    : Port(cs_ptr, port_nr, start, stop, excite, PortNamePrefix, delay, priority),
      R(R),
      exc_ny(CheckNyDir(exc_dir))
{
    _CSX csx(cs_ptr);

    // Determine the excitation direction sign along exc_ny.
    direction = ((stop[exc_ny] - start[exc_ny]) >= 0) ? 1 : -1;
    if (start[exc_ny] == stop[exc_ny]) {
        throw std::runtime_error("LumpedPort: start and stop may not be identical in excitation direction");
    }

    // Add lumped resistor element (or metal if R == 0).
    if (R > 0) {
        auto lumped_R = csx.AddLumpedElement(lbl_temp + "resist", exc_ny, true);
        lumped_R->SetResistance(R);

        port_props.push_back(lumped_R);
        _CSProperties(lumped_R).AddBox(start, stop, priority);
        
    } else if (R == 0) {
        auto lumped_R = csx.AddMetal(lbl_temp + "resist");

        port_props.push_back(lumped_R);
        _CSProperties(lumped_R).AddBox(start, stop, priority);
    }

    // Add excitation if requested.
    if (excite != 0) {
        std::array<double, 3> exc_vec = {0.0, 0.0, 0.0};
        exc_vec[exc_ny] = -1.0 * direction * excite;
        auto exc = csx.AddExcitation(lbl_temp + "excite", 0, exc_vec);
        exc.SetDelay(delay);
        exc.AddBox(start, stop, priority);
        port_props.push_back(exc);
    }

    // Setup voltage probe.
    U_filenames.push_back(lbl_temp + "ut");
    std::array<double, 3> u_start, u_stop;
    for (size_t i = 0; i < 3; ++i) {
        u_start[i] = 0.5 * (start[i] + stop[i]);
        u_stop[i]  = u_start[i];
    }
    u_start[exc_ny] = start[exc_ny];
    u_stop[exc_ny]  = stop[exc_ny];
    auto u_probe = csx.AddProbe(U_filenames.back(), 0); // p_type=0 for voltage probe, weight -1.
    u_probe.SetWeighting(-1);
    u_probe.AddBox(u_start, u_stop);
    port_props.push_back(u_probe);

    // Setup current probe.
    I_filenames.push_back(lbl_temp + "it");
    std::array<double, 3> i_start = start;
    std::array<double, 3> i_stop  = stop;
    i_start[exc_ny] = 0.5 * (start[exc_ny] + stop[exc_ny]);
    i_stop[exc_ny]  = i_start[exc_ny];
    auto i_probe = csx.AddProbe(I_filenames.back(), 1); // p_type=1 for current probe.
    i_probe.SetWeighting(direction);
    i_probe.SetNormalDir(exc_ny);
    i_probe.AddBox(i_start, i_stop);
    port_props.push_back(i_probe);
}

void LumpedPort::CalcPort(const std::string& sim_path,
    const std::vector<double>& freq,
    std::optional<std::complex<double>> ref_impedance,
    double ref_plane_shift,
    const std::string& signal_type)
{
//     if ref_impedance is None:
//         self.Z_ref = self.R
//     if ref_plane_shift is not None:
//         Warning('A lumped port does not support a reference plane shift! Ignoring...')
//     super(LumpedPort, self).CalcPort(sim_path, freq, ref_impedance, ref_plane_shift, signal_type)

    // Set the port reference impedance.
    if (!ref_impedance)
        Z_ref = R;

    Port::CalcPort(sim_path, freq, Z_ref, 0.0, signal_type);
}

// --- MSLPort Implementation ---

// Python-like
MSLPort::MSLPort(ContinuousStructure* cs_ptr,
                 int port_nr,
                 CSProperties* metal_prop,  // expected to be a pointer to a metal property object with an AddBox() method
                 const std::array<double, 3>& start,
                 const std::array<double, 3>& stop,
                 int prop_dir,
                 int exc_dir,
                 double excite,
                 double FeedShift,
                 double MeasPlaneShift,
                 double Feed_R,
                 const std::string& PortNamePrefix,
                 double delay,
                 int priority)
    : Port(cs_ptr, port_nr, start, stop, excite, PortNamePrefix, delay, priority),
      exc_ny(exc_dir),
      prop_ny(prop_dir),
      direction((stop[prop_dir] - start[prop_dir]) >= 0 ? 1 : -1),
      feed_shift(FeedShift),
      measplane_shift((!std::isnan(MeasPlaneShift)) ? MeasPlaneShift : 0.5 * std::abs(stop[prop_dir] - start[prop_dir])),
      feed_R(Feed_R)
{
    _CSX csx(cs_ptr);

    // Determine upside_down based on the excitation direction.
    upside_down = ((stop[exc_ny] - start[exc_ny]) >= 0) ? 1 : -1;

    // Ensure start and stop are not identical in any dimension.
    for (size_t i = 0; i < start.size(); ++i) {
        if (start[i] == stop[i])
            throw std::runtime_error("MSLPort: start and stop may not be identical in all dimensions");
    }
    // Excitation and propagation directions must differ.
    if (exc_ny == prop_ny)
        throw std::runtime_error("MSLPort: Excitation direction must not be equal to propagation direction");

    // Add metal MSL plane.
    std::array<double, 3> MSL_start = start;
    std::array<double, 3> MSL_stop  = stop;
    MSL_stop[exc_ny] = MSL_start[exc_ny];
    _CSProperties(metal_prop).AddBox(MSL_start, MSL_stop, priority);

    // Get grid lines along the propagation direction.
    _CSRectGrid grid = csx->GetGrid();
    std::vector<double> prop_lines = grid.GetLines(prop_ny);
    if (prop_lines.size() <= 5)
        throw std::runtime_error("At least 5 grid lines in propagation direction required!");

    // Determine measurement plane position.
    double measplane_pos = start[prop_ny] + measplane_shift * direction;
    size_t meas_pos_idx = interp1_nearest(measplane_pos, prop_lines);
    meas_pos_idx = std::clamp(meas_pos_idx, 1ULL, prop_lines.size() - 2);

    // Update measplane_shift using the grid.
    measplane_shift = std::abs(start[prop_ny] - prop_lines[meas_pos_idx]);

    // Determine voltage probe positions.
    std::vector<size_t> prope_idx = {meas_pos_idx - 1, meas_pos_idx, meas_pos_idx + 1};
    if (direction < 0)
        std::reverse(prope_idx.begin(), prope_idx.end());
    std::vector<double> u_prope_pos;
    for (auto idx : prope_idx)
        u_prope_pos.push_back(prop_lines[idx]);
    U_delta.clear();
    for (size_t i = 1; i < u_prope_pos.size(); ++i)
        U_delta.push_back(u_prope_pos[i] - u_prope_pos[i - 1]);

    // Create voltage probe boxes.
    const std::vector<std::string> suffix = {"A", "B", "C"};
    // Define the mid-point between start and stop.
    std::array<double, 3> mid_point{};
    for (size_t i = 0; i < start.size(); ++i)
        mid_point[i] = 0.5 * (start[i] + stop[i]);
    for (size_t n = 0; n < u_prope_pos.size(); ++n) {
        std::array<double, 3> u_start = mid_point;
        std::array<double, 3> u_stop  = mid_point;
        u_start[prop_ny] = u_prope_pos[n];
        u_stop[prop_ny]  = u_prope_pos[n];
        u_start[exc_ny] = start[exc_ny];
        u_stop[exc_ny]  = stop[exc_ny];
        std::string u_name = lbl_temp + "ut" + suffix[n];
        U_filenames.push_back(u_name);
        auto u_probe = csx.AddProbe(u_name, 0); // 0 indicates a voltage probe.
        _CSProperties(u_probe).AddBox(u_start, u_stop);
        port_props.push_back(u_probe);
    }

    // Determine current probe positions.
    std::vector<double> i_prope_pos;
    if (u_prope_pos.size() >= 2) {
        for (size_t n = 0; n < 2; ++n) {
            double pos = u_prope_pos[n] + 0.5 * ((n + 1 < u_prope_pos.size() ? u_prope_pos[n + 1] : u_prope_pos[n]) - u_prope_pos[n]);
            i_prope_pos.push_back(pos);
        }
    } else {
        throw std::runtime_error("Not enough voltage probe positions for current probe calculation.");
    }
    I_delta.clear();
    if (i_prope_pos.size() == 2)
        I_delta.push_back(i_prope_pos[1] - i_prope_pos[0]);

    // Create current probe boxes.
    for (size_t n = 0; n < i_prope_pos.size(); ++n) {
        std::array<double, 3> i_start = start;
        std::array<double, 3> i_stop  = stop;
        i_start[prop_ny] = i_prope_pos[n];
        i_stop[prop_ny]  = i_prope_pos[n];
        std::string i_name = lbl_temp + "it" + suffix[n];
        I_filenames.push_back(i_name);
        auto i_probe = csx.AddProbe(i_name, 1); // 1 indicates a current probe.
        i_probe->SetWeighting(direction);
        i_probe->SetNormalDir(prop_ny);
        _CSProperties(i_probe).AddBox(i_start, i_stop);
        port_props.push_back(i_probe);
    }

    // Add excitation if requested.
    if (excite != 0) {
        double target = start[prop_ny] + feed_shift * direction;
        size_t exc_idx = interp1_nearest(target, prop_lines);
        std::array<double, 3> exc_start = start;
        std::array<double, 3> exc_stop  = stop;
        exc_start[prop_ny] = prop_lines[exc_idx];
        exc_stop[prop_ny]  = prop_lines[exc_idx];
        std::array<double, 3> exc_vec{};
        exc_vec[exc_ny] = -1 * upside_down * excite;
        auto exc = csx.AddExcitation(lbl_temp + "excite", 0, exc_vec);
        exc->SetDelay(delay);
        _CSProperties(exc).AddBox(exc_start, exc_stop, priority);
        port_props.push_back(exc);
    }

    // Add feed resistor if applicable.
    if (feed_R >= 0 && !std::isinf(feed_R)) {
        std::array<double, 3> R_start = start;
        std::array<double, 3> R_stop  = stop;
        R_stop[prop_ny] = R_start[prop_ny];
        if (feed_R == 0) {
            _CSProperties(metal_prop).AddBox(R_start, R_stop, priority);
        } else {
            auto lumped_R = csx.AddLumpedElement(lbl_temp + "resist", exc_ny, true);
            lumped_R->SetResistance(feed_R);
            _CSProperties(lumped_R).AddBox(R_start, R_stop, priority);
            port_props.push_back(lumped_R);
        }
    }
}

// MatLab-like
MSLPort::MSLPort(ContinuousStructure* cs_ptr,
                 int priority,
                 int port_nr,
                 CSProperties* metal_prop,
                 const std::array<double, 3>& start,
                 const std::array<double, 3>& stop,
                 int prop_dir,
                 const std::array<double, 3>& evec,
                 bool excite_port,
                 double FeedShift,
                 double MeasPlaneShift,
                 double Feed_R,
                 const std::string& PortNamePrefix,
                 double delay)
    : Port(cs_ptr, port_nr, start, stop, (excite_port ? 1.0 : 0.0), PortNamePrefix, delay, priority),
      prop_ny(prop_dir),
      feed_shift(FeedShift),
      measplane_shift((!std::isnan(MeasPlaneShift)) ? MeasPlaneShift : 0.5 * std::abs(stop[prop_dir] - start[prop_dir])),
      feed_R(Feed_R)
{
    _CSX csx(cs_ptr);

    // Normalize start/stop.
    std::array<double, 3> nstart, nstop;
    for (int i = 0; i < 3; i++) {
        nstart[i] = std::min(start[i], stop[i]);
        nstop[i]  = std::max(start[i], stop[i]);
    }

    // Determine the index for the excitation (height) direction.
    int exc_dir = 0;
    for (int i = 0; i < 3; i++) {
        if (std::abs(evec[i]) > 1e-12) {
            exc_dir = i;
            break;
        }
    }
    // If evec is negative in its nonzero component, reverse the port’s excitation sign.
    if (evec[exc_dir] < 0)
        excite *= -1;

    // Determine the index for the width. (All three indices 0+1+2=3)
    int width_dir = 3 - prop_dir - exc_dir;

    // Determine propagation direction sign.
    int direction = ((stop[prop_dir] - start[prop_dir]) > 0) ? 1 : -1;

    // Create MSL metal plane. Set the metal’s height to the start coordinate.
    std::array<double, 3> MSL_start = start;
    std::array<double, 3> MSL_stop  = stop;
    MSL_stop[exc_dir] = MSL_start[exc_dir];
    _CSProperties(metal_prop).AddBox(MSL_start, MSL_stop, priority);

    // Get grid lines along each direction.
    _CSRectGrid grid = csx->GetGrid();
    std::vector<double> prop_lines  = grid.GetLines(prop_dir);
    std::vector<double> width_lines = grid.GetLines(width_dir);
    std::vector<double> height_lines = grid.GetLines(exc_dir);

    // Determine measurement plane position along propagation (shifted from start).
    double measplane_pos = start[prop_dir] + direction * measplane_shift;
    size_t closest_idx = interp1_nearest(measplane_pos, prop_lines);
    closest_idx = std::clamp(closest_idx, 1ULL, prop_lines.size() - 2);

    // Calculate voltage probe positions (three positions).
    std::vector<double> u_probe_positions = {prop_lines[closest_idx - 1],
                                             prop_lines[closest_idx],
                                             prop_lines[closest_idx + 1]};
    if (direction < 0)
        std::reverse(u_probe_positions.begin(), u_probe_positions.end());

    // Determine the center of the width (using grid lines from normalized start/stop).
    double center_width = 0.5 * (nstart[width_dir] + nstop[width_dir]);
    size_t width_idx = interp1_nearest(center_width, width_lines);
    double mid_width = width_lines[width_idx];

    // Create voltage probe boxes (three boxes: A, B, C).
    const std::array<std::string, 3> suffix = { "A", "B", "C" };
    for (size_t n = 0; n < u_probe_positions.size(); n++) {
        std::array<double, 3> v_start, v_stop;
        // Set midpoint for all coordinates.
        for (int i = 0; i < 3; i++) {
            v_start[i] = 0.5 * (start[i] + stop[i]);
            v_stop[i]  = v_start[i];
        }
        // Set propagation coordinate to voltage probe position.
        v_start[prop_dir] = u_probe_positions[n];
        v_stop[prop_dir]  = u_probe_positions[n];
        // Set width coordinate to grid center.
        v_start[width_dir] = mid_width;
        v_stop[width_dir]  = mid_width;
        // Height probe box runs from start to stop in height.
        v_start[exc_dir] = start[exc_dir];
        v_stop[exc_dir]  = stop[exc_dir];

        std::string u_name = lbl_temp + "ut" + suffix[n];
        U_filenames.push_back(u_name);
        auto u_probe = csx.AddProbe(u_name, 0);
        _CSProperties(u_probe).AddBox(v_start, v_stop, priority);
        port_props.push_back(u_probe);
    }

    // Create current probe positions as the averages of adjacent voltage probe positions.
    std::vector<double> i_probe_positions;
    i_probe_positions.push_back(0.5 * (u_probe_positions[0] + u_probe_positions[1]));
    i_probe_positions.push_back(0.5 * (u_probe_positions[1] + u_probe_positions[2]));

    // Create current probe boxes (two boxes).
    for (size_t n = 0; n < i_probe_positions.size(); n++) {
        std::array<double, 3> i_start_arr, i_stop_arr;
        // Start with midpoint values.
        for (int i = 0; i < 3; i++) {
            i_start_arr[i] = 0.5 * (start[i] + stop[i]);
            i_stop_arr[i]  = i_start_arr[i];
        }
        // Set propagation coordinate.
        i_start_arr[prop_dir] = i_probe_positions[n];
        i_stop_arr[prop_dir]  = i_probe_positions[n];

        // Increase by 0.5 grid cells in the width direction.
        int w_idx0 = interp1_nearest(nstart[width_dir], width_lines);
        int w_idx1 = interp1_nearest(nstop[width_dir], width_lines);
        i_start_arr[width_dir] = width_lines[w_idx0] - (width_lines[w_idx0] - width_lines[w_idx0 - 1]) / 2;
        i_stop_arr[width_dir] = width_lines[w_idx1] + (width_lines[w_idx1 + 1] - width_lines[w_idx1]) / 2;

        // Increase by 1.5 grid cells in the height direction.
        int h_idx = interp1_nearest(start[exc_dir], height_lines);
        i_start_arr[exc_dir] = height_lines[h_idx - 1] - (height_lines[h_idx - 1] - height_lines[h_idx - 2]) / 2;
        i_stop_arr[exc_dir] = height_lines[h_idx + 1] + (height_lines[h_idx + 2] - height_lines[h_idx + 1]) / 2;

        std::string i_name = lbl_temp + "it" + std::string(1, ("AB")[n]);
        I_filenames.push_back(i_name);
        auto i_probe = csx.AddProbe(i_name, 1);
        i_probe->SetWeighting(direction);
        i_probe->SetNormalDir(prop_dir);
        _CSProperties(i_probe).AddBox(i_start_arr, i_stop_arr, priority);
        port_props.push_back(i_probe);
    }

    // Fill I_delta correctly: use the difference between the two current probe positions.
    I_delta.clear();
    if (i_probe_positions.size() == 2)
        I_delta.push_back(i_probe_positions[1] - i_probe_positions[0]);

    // Add excitation if enabled.
    if (excite_port) {
        double exc_target = start[prop_dir] + feed_shift * direction;
        size_t exc_idx = interp1_nearest(exc_target, prop_lines);
        
        std::array<double, 3> ex_start_arr, ex_stop_arr;
        for (int i = 0; i < 3; i++) {
            ex_start_arr[i] = nstart[i];
            ex_stop_arr[i]  = nstop[i];
        }
        ex_start_arr[prop_dir] = prop_lines[exc_idx];
        ex_stop_arr[prop_dir]  = prop_lines[exc_idx];
        std::array<double, 3> exc_vec{};
        // Excitation is applied in the height (e-field) direction.
        exc_vec[exc_dir] = -1 * (((stop[exc_dir] - start[exc_dir]) >= 0) ? 1 : -1) * (excite ? excite : 0);
        auto exc_obj = csx.AddExcitation(lbl_temp + "excite", 0, exc_vec);
        exc_obj->SetDelay(delay);
        _CSProperties(exc_obj).AddBox(ex_start_arr, ex_stop_arr, priority);
        port_props.push_back(exc_obj);
    }

    // Add lumped resistor (or metal if Feed_R==0) at the start of MSL line.
    std::array<double, 3> r_start = start;
    std::array<double, 3> r_stop  = nstop;
    r_stop[prop_dir] = r_start[prop_dir];
    if (feed_R > 0 && !std::isinf(feed_R)) {
        auto lumped_R = csx.AddLumpedElement(lbl_temp + "resist", exc_dir, true);
        lumped_R->SetResistance(feed_R);
        _CSProperties(lumped_R).AddBox(r_start, r_stop, priority);
        port_props.push_back(lumped_R);
    } else if (std::isinf(feed_R)) {
        // Open port: do nothing.
    } else if (feed_R == 0) {
        _CSProperties(metal_prop).AddBox(r_start, r_stop, priority);
    }
}

void MSLPort::ReadUIData(const std::string &sim_path,
                         const std::vector<double> &freq,
                         const std::string &signal_type)
{
    // Create UI_data object for voltage probes.
    u_data = std::make_unique<UI_data>(U_filenames, sim_path, freq, signal_type);
    uf_tot = u_data->ui_f_val[1];
    ut_tot = u_data->ui_val[1];

    // Create UI_data object for current probes.
    i_data = std::make_unique<UI_data>(I_filenames, sim_path, freq, signal_type);

    if_tot.resize(freq.size());
    for (size_t i = 0; i < uf_tot.size(); ++i)
        if_tot[i] = 0.5 * (i_data->ui_f_val[0][i] + i_data->ui_f_val[1][i]);

    it_tot.resize(ut_tot.size());
    for (size_t i = 0; i < ut_tot.size(); ++i)
        it_tot[i] = 0.5 * (i_data->ui_val[0][i] + i_data->ui_val[1][i]);

    // Get grid unit.
    double unit = CSX_ptr->GetGrid()->GetDeltaUnit();

    // Compute dEt = (ui_f_val[2] - ui_f_val[0]) / (sum(|U_delta|)*unit)
    double sumU_delta = std::accumulate(U_delta.begin(), U_delta.end(), 0.0,
        [](double s, double val){ return s + std::abs(val); });
    std::vector<std::complex<double>> dEt(freq.size());
    for (size_t i = 0; i < freq.size(); ++i)
        dEt[i] = (u_data->ui_f_val[2][i] - u_data->ui_f_val[0][i]) / (sumU_delta * unit);

    // Ht is the space-averaged current value.
    auto Ht = if_tot;

    // Compute dHt = (ui_f_val[1] - ui_f_val[0]) / (|I_delta[0]|*unit)
    double absI_delta = std::abs(I_delta[0]);
    std::vector<std::complex<double>> dHt(freq.size());
    for (size_t i = 0; i < freq.size(); ++i)
        dHt[i] = (i_data->ui_f_val[1][i] - i_data->ui_f_val[0][i]) / (absI_delta * unit);

    // Calculate beta = sqrt( - dEt * dHt / (Ht * Et) )
    beta.resize(freq.size());
    for (size_t i = 0; i < freq.size(); ++i) {
        std::complex<double> num = -dEt[i] * dHt[i];
        std::complex<double> den = uf_tot[i] * u_data->ui_f_val[1][i];
        std::complex<double> bval = std::sqrt(num / den);
        if (bval.real() < 0)
            bval = -bval;
        beta[i] = bval;
    }

    // Determine Z_ref = sqrt(Et * dEt / (Ht * dHt)); here we use the first frequency point.
    std::complex<double> ZL = std::sqrt((u_data->ui_f_val[1][0] * dEt[0]) / (Ht[0] * dHt[0]));
    Z_ref = ZL.real();
}

// --- WaveguidePort Implementation ---

WaveguidePort::WaveguidePort(ContinuousStructure* cs_ptr,
                             int port_nr,
                             const std::array<double, 3>& start,
                             const std::array<double, 3>& stop,
                             int exc_dir,
                             const std::array<std::string, 3>& E_WG_func,
                             const std::array<std::string, 3>& H_WG_func,
                             double kc,
                             double excite,
                             const std::string& PortNamePrefix,
                             double delay,
                             int priority)
    : Port(cs_ptr, port_nr, start, stop, excite, PortNamePrefix, delay, priority),
      exc_ny(exc_dir),
      kc(kc),
      E_func(E_WG_func),
      H_func(H_WG_func),
      ref_index(1.0)
{
    _CSX csx(cs_ptr);

    // Determine the other two directions based on the excitation direction.
    ny_P  = (exc_ny + 1) % 3;
    ny_PP = (exc_ny + 2) % 3;
    
    // Calculate measurement plane positions.
    std::array<double, 3> m_start = start;
    std::array<double, 3> m_stop  = stop;
    // Place both voltage and current probes in the plane defined by exc_ny (voltage/current planes share the same box).
    m_start[exc_ny] = m_stop[exc_ny];

    // Add excitation if requested.
    if (excite != 0) {
        std::array<double, 3> e_start = start;
        std::array<double, 3> e_stop  = stop;
        e_stop[exc_ny] = e_start[exc_ny];  // Ensure the excitation is applied in a plane.
        // Initialize an excitation vector (all zeros).
        std::array<double, 3> e_vec = {0.0, 0.0, 0.0};
        // Create the excitation with type 0.
        auto exc = csx.AddExcitation(lbl_temp + "excite", 0, e_vec);
        exc->SetDelay(delay);
        // Set the weighting function for the mode (provided as strings).
        exc.SetWeightFunction(E_func);
        exc.AddBox(e_start, e_stop, priority);
        port_props.push_back(exc);
    }

    // Create the voltage probe using a mode-based probe (p_type 10).
    U_filenames.push_back(lbl_temp + "ut");
    auto u_probe = csx.AddProbe(U_filenames.back(), 10);
    u_probe.SetModeFunction(E_func);
    u_probe.AddBox(m_start, m_stop);
    port_props.push_back(u_probe);

    // Create the current probe using a mode-based probe (p_type 11).
    I_filenames.push_back(lbl_temp + "it");
    auto i_probe = csx.AddProbe(I_filenames.back(), 11);
    i_probe.SetModeFunction(H_func);
    // Set weighting: the sign is determined by the direction of excitation.
    double weight = (stop[exc_ny] - start[exc_ny]) >= 0 ? 1.0 : -1.0;
    i_probe.SetWeighting(weight);
    i_probe.SetNormalDir(exc_ny);
    i_probe.AddBox(m_start, m_stop);
    port_props.push_back(i_probe);
}

void WaveguidePort::CalcPort(const std::string& sim_path,
    const std::vector<double>& freq,
    std::optional<std::complex<double>> ref_impedance,
    double ref_plane_shift,
    const std::string& signal_type)
{
    // Compute k for each frequency: k = 2π*freq*ref_index/C0
    std::vector<double> k_vec(freq.size(), 0.0);
    beta.resize(freq.size());
    ZL.resize(freq.size());
    for (size_t i = 0; i < freq.size(); ++i) {
        k_vec[i] = 2.0 * M_PI * freq[i] * ref_index / C0;
        beta[i] = std::sqrt(k_vec[i] * k_vec[i] - kc * kc);
        ZL[i] = k_vec[i] * Z0 / beta[i]; // Analytic waveguide impedance
    }

    // If no reference impedance was provided (<= 0), default to the analytic value at the first frequency.
    if (!ref_impedance)
    {
        Z_ref = ZL[0];
        std::cout << "WaveguidePort: ref_impedance not set, using " << Z_ref << std::endl;
    }
    else
    {
        Z_ref = *ref_impedance;
    }

    // Call base class CalcPort to perform UI-data reading and further port calculations.
    Port::CalcPort(sim_path, freq, Z_ref, ref_plane_shift, signal_type);
}

// --- RectWGPort Implementation ---


// Helper function to compute the cutoff parameter kc.
static double computeKC(double a, double b, const std::string& mode_name)
{
    // Extract mode numbers M and N from the mode string.
    double M = static_cast<double>(mode_name[2] - '0');
    double N = static_cast<double>(mode_name[3] - '0');
    return std::sqrt(std::pow(M * M_PI / a, 2) + std::pow(N * M_PI / b, 2));
}

// Helper function to create the electric field function strings.
static std::array<std::string, 3> make_EWG(ContinuousStructure* csx,
    const std::array<double, 3>& start,
    int exc_dir,
    double a,
    double b,
    const std::string& mode_name)
{
    double M = static_cast<double>(mode_name[2] - '0');
    double N = static_cast<double>(mode_name[3] - '0');

    // Determine probe coordinate indices.
    int ny_P = (exc_dir + 1) % 3;
    int ny_PP = (exc_dir + 2) % 3;
    char coords[3] = { 'x', 'y', 'z' };

    // Create coordinate name strings.
    std::string name_P = (start[ny_P] != 0.0) ? "(" + std::string(1, coords[ny_P]) + "-" + std::to_string(start[ny_P]) + ")" : std::string(1, coords[ny_P]);
    std::string name_PP = (start[ny_PP] != 0.0) ? "(" + std::string(1, coords[ny_PP]) + "-" + std::to_string(start[ny_PP]) + ")" : std::string(1, coords[ny_PP]);

    std::array<std::string, 3> E_WG_func = { "0", "0", "0" };
    // Build E field functions based on the mode numbers.
    if (N > 0)
    {
        E_WG_func[ny_P] = std::to_string(N / b) + "*cos(" +
            std::to_string(M * M_PI / a) + "*" + name_P + ")*sin(" +
            std::to_string(N * M_PI / b) + "*" + name_PP + ")";
    }
    if (M > 0)
    {
        E_WG_func[ny_PP] = std::to_string(-M / a) + "*sin(" +
            std::to_string(M * M_PI / a) + "*" + name_P + ")*cos(" +
            std::to_string(N * M_PI / b) + "*" + name_PP + ")";
    }
    return E_WG_func;
}

// Helper function to create the magnetic field function strings.
static std::array<std::string, 3> make_HWG(ContinuousStructure* csx,
    const std::array<double, 3>& start,
    int exc_dir,
    double a,
    double b,
    const std::string& mode_name)
{
    double M = static_cast<double>(mode_name[2] - '0');
    double N = static_cast<double>(mode_name[3] - '0');

    int ny_P = (exc_dir + 1) % 3;
    int ny_PP = (exc_dir + 2) % 3;
    char coords[3] = { 'x', 'y', 'z' };

    std::string name_P = (start[ny_P] != 0.0) ? "(" + std::string(1, coords[ny_P]) + "-" + std::to_string(start[ny_P]) + ")" : std::string(1, coords[ny_P]);
    std::string name_PP = (start[ny_PP] != 0.0) ? "(" + std::string(1, coords[ny_PP]) + "-" + std::to_string(start[ny_PP]) + ")" : std::string(1, coords[ny_PP]);

    std::array<std::string, 3> H_WG_func = { "0", "0", "0" };
    if (M > 0)
    {
        H_WG_func[ny_P] = std::to_string(M / a) + "*sin(" +
            std::to_string(M * M_PI / a) + "*" + name_P + ")*cos(" +
            std::to_string(N * M_PI / b) + "*" + name_PP + ")";
    }
    if (N > 0)
    {
        H_WG_func[ny_PP] = std::to_string(N / b) + "*cos(" +
            std::to_string(M * M_PI / a) + "*" + name_P + ")*sin(" +
            std::to_string(N * M_PI / b) + "*" + name_PP + ")";
    }
    return H_WG_func;
}


RectWGPort::RectWGPort(ContinuousStructure *csx,
                       int port_nr,
                       const std::array<double, 3> &start,
                       const std::array<double, 3> &stop,
                       int exc_dir,
                       double a,
                       double b,
                       const std::string &mode_name,
                       double excite,
                       const std::string &PortNamePrefix,
                       double delay,
                       int priority)
    : WaveguidePort(csx,
                    port_nr,
                    start,
                    stop,
                    exc_dir,
                    /* E_WG_func */ make_EWG(csx, start, exc_dir, a, b, mode_name),
                    /* H_WG_func */ make_HWG(csx, start, exc_dir, a, b, mode_name),
                    /* kc */ computeKC(a, b, mode_name),
                    excite,
                    PortNamePrefix,
                    delay,
                    priority)
{
    // Only TE modes are supported.
    if (mode_name.size() != 4)
    {
        throw std::invalid_argument("Invalid mode definition");
    }
    if (!(mode_name[0] == 'T' && mode_name[1] == 'E'))
    {
        throw std::invalid_argument("Currently only TE-modes are supported! Mode provided: " + mode_name);
    }
}
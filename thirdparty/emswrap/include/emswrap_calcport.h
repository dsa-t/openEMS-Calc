#pragma once

#include <emswrap_ports.h>

#include <string>
#include <vector>
#include <complex>
#include <stdexcept>
#include <cmath>

// Forward declarations (assumed to be defined elsewhere)
// Implementation of calcPort based on the MATLAB code

static Port* calcPort(Port* port,
                      const std::string& sim_path,
                      const std::vector<double>& freq,
                      double ref_impedance = -1,
                      double ref_plane_shift = 0.0,
                      const std::string& signal_type = "pulse")
{
    if (!port)
        return port;

    // Call the virtual function; the proper overridden method will be used.
    port->CalcPort(sim_path, freq, ref_impedance, ref_plane_shift, signal_type);

    // Make sure the data arrays have the expected size
    size_t n = freq.size();
    if (port->uf_inc.size() != n || port->if_inc.size() != n ||
        port->uf_ref.size() != n || port->if_ref.size() != n ||
        port->uf_tot.size() != n || port->if_tot.size() != n)
    {
        throw std::runtime_error("Mismatch between frequency vector size and port data size.");
    }

    // Resize power vectors to match frequency vector size.
    port->P_inc.resize(n);
    port->P_ref.resize(n);
    port->P_acc.resize(n);

    // Calculate power parameters:
    //   P_inc = 0.5 * real( uf_inc * conj(if_inc) )
    //   P_ref = 0.5 * real( uf_ref * conj(if_ref) )
    //   P_acc = 0.5 * real( uf_tot * conj(if_tot) )
    for (size_t i = 0; i < n; ++i) {
        port->P_inc[i] = 0.5 * std::real(port->uf_inc[i] * std::conj(port->if_inc[i]));
        port->P_ref[i] = 0.5 * std::real(port->uf_ref[i] * std::conj(port->if_ref[i]));
        port->P_acc[i] = 0.5 * std::real(port->uf_tot[i] * std::conj(port->if_tot[i]));
    }

    return port;
}

// Overload for a vector of Port pointers.
static std::vector<Port*> calcPort(const std::vector<Port*>& ports,
                                   const std::string& sim_path,
                                   const std::vector<double>& freq,
                                   double ref_impedance = -1,
                                   double ref_plane_shift = 0.0,
                                   const std::string& signal_type = "pulse")
{
    std::vector<Port*> updatedPorts;
    updatedPorts.reserve(ports.size());
    for (Port* p : ports) {
        updatedPorts.push_back(calcPort(p, sim_path, freq, ref_impedance, ref_plane_shift, signal_type));
    }
    return updatedPorts;
}
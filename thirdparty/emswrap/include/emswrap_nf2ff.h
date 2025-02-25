#pragma once

#include <vector>
#include <array>
#include <string>
#include <complex>
#include <filesystem>
#include <stdexcept>

// Forward declaration of CSX (assumed to be defined elsewhere)
class ContinuousStructure;
class CSPropDumpBox;

class _nf2ff_results;

class _nf2ff {
public:
    _nf2ff(ContinuousStructure* CSX, const std::string& name, const std::array<double, 3>& start, const std::array<double, 3>& stop,
          const std::vector<bool>& directions = std::vector<bool>(6, true),
          const std::vector<int>& mirror = std::vector<int>(6, 0),
          const std::vector<double>& frequency = {});

    _nf2ff_results CalcNF2FF(const std::string& sim_path,
                            const std::vector<double>& freq,
                            const std::vector<double>& theta,
                            const std::vector<double>& phi,
                            double radius = 1.0,
                            const std::array<double, 3>& center = {0, 0, 0},
                            const std::string& outfile = "",
                            bool read_cached = false,
                            int verbose = 2);

    _nf2ff_results CalcNF2FF(const std::string &sim_path,
                             const std::vector<double> &freq,
                             const std::vector<double> &theta,
                             const std::vector<double> &phi,
                             const std::string &outfile,
                             bool read_cached = false,
                             int verbose = 2);

private:
    ContinuousStructure* csx;
    std::string name;
    std::array<double, 3> start;
    std::array<double, 3> stop;
    std::vector<bool> directions;
    std::vector<int> mirror;
    std::vector<double> freq;
    std::string e_file;
    std::string h_file;
    
    CSPropDumpBox* e_dump;
    CSPropDumpBox* h_dump;

};

class _nf2ff_results {
public:
    explicit _nf2ff_results(const std::string& fn);

    std::string fn;
    std::vector<double> theta;
    std::vector<double> phi;
    double r;
    std::vector<double> freq;
    std::vector<double> Dmax;
    std::vector<double> Prad;
    std::vector<std::vector<std::vector<std::complex<double>>>> E_theta;  // [freq][theta][phi]
    std::vector<std::vector<std::vector<std::complex<double>>>> E_phi;    // [freq][theta][phi]
    std::vector<std::vector<std::vector<double>>> E_norm;                 // [freq][theta][phi]
    std::vector<std::vector<std::vector<std::complex<double>>>> E_cprh;   // [freq][theta][phi]
    std::vector<std::vector<std::vector<std::complex<double>>>> E_cplh;   // [freq][theta][phi]
    std::vector<std::vector<std::vector<double>>> P_rad;                  // [freq][theta][phi]

private:
    void load_hdf5_data();
    // static std::vector<std::complex<double>> read_complex_hdf5(const H5::H5File& file, const std::string& path_real, const std::string& path_imag);
};
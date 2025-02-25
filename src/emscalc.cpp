#include <emswrap.h>

// # Porting examples:
// Python:
/*
    FDTD = openEMS(EndCriteria=1e-5)
    FDTD.SetGaussExcite( f0, fc )
    FDTD.SetBoundaryCond( ['PML_8', 'PML_8', 'PML_8', 'PML_8', 'PML_8', 'PML_8'] )
    ...
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)

    double max_res = C0 / f_stop / unit / 20; # cell size: lambda/20
    double ratio = 1.4;

    mesh.SetLines('x', [-SimBox/2, 0, SimBox/2])
    mesh.SmoothMeshLines('x', max_res, ratio )
    mesh.SetLines('y', mesh.GetLines('x'))
    mesh.SetLines('z', mesh.GetLines('x'))
*/
// C++:
/*
    _openEMS FDTD;
    FDTD.SetEndCriteria(1e-5);
    FDTD.SetGaussExcite(f0, fc);
    FDTD.SetBoundaryCond({ "PML_8","PML_8","PML_8","PML_8","PML_8","PML_8" });
    ...
    _CSX CSX(new ContinuousStructure());
    FDTD.SetCSX(CSX);
    _CSRectGrid mesh(CSX.GetGrid());
    mesh.SetDeltaUnit(unit);

    double max_res = C0 / f_stop / unit / 20; // cell size: lambda/20
    double ratio = 1.4;

    mesh.SetLines('x', { -SimBox/2, 0, SimBox/2 });
    mesh.SmoothMeshLines('x', max_res, ratio );
    mesh.SetLines('y', mesh.GetLines('x'));
    mesh.SetLines('z', mesh.GetLines('x'));
    ...
*/

// Matlab:
/*
    FDTD = InitFDTD();
    FDTD = SetGaussExcite( FDTD, f0, fc );
    BC   = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'};
    FDTD = SetBoundaryCond( FDTD, BC );
    ...
    CSX = InitCSX();
    resolution = c0/(f0+fc)/unit /15; % resolution of lambda/15
    mesh.x = SmoothMeshLines( [-width/2, width/2, -MSL_width/2, MSL_width/2], resolution ); % create smooth lines from fixed lines
    mesh.y = SmoothMeshLines( [linspace(0,MSL_height,5) MSL_height+1 height], resolution );
    mesh.z = SmoothMeshLines( [0 length], resolution );
    CSX = DefineRectGrid( CSX, unit, mesh );
*/
// C++:
/*
    _openEMS FDTD;
    FDTD.SetEndCriteria(1e-5);
    FDTD.SetGaussExcite(f0, fc);
    FDTD.SetBoundaryCond({ "PML_8","PML_8","PML_8","PML_8","PML_8","PML_8" });
    ...
    _CSX CSX(new ContinuousStructure());
    FDTD.SetCSX(CSX);
    _CSRectGrid mesh(CSX.GetGrid());

    double resolution = c0 / (f0 + fc) / unit / 15; // resolution of lambda/15
    mesh.SmoothMeshLines('x', { -width/2, width/2, -MSL_width/2, MSL_width/2 }, resolution ); // create smooth lines from fixed lines
    mesh.SmoothMeshLines('y', { linspace(0,MSL_height,5), MSL_height+1, height }, resolution );
    mesh.SmoothMeshLines('z', { 0, length }, resolution );
    mesh.SetDeltaUnit(unit);
    ...
*/

int main(int argc, char* argv[]) {

    std::string Sim_Path = "RCS_Sphere";

    // Simulation parameters
    double unit = 1e-3; // mm
    double sphere_rad = 200;
    double inc_angle = 0; // deg
    double SimBox = 1200;
    double PW_Box = 750;
    double f_start = 50e6;
    double f_stop  = 1000e6;
    double f0      = 500e6;

    // FDTD setup
    _openEMS FDTD;
    FDTD.SetEndCriteria(1e-5);

    FDTD.DebugCSX();
    FDTD.DebugBox();
    FDTD.DebugMaterial();
    FDTD.DebugOperator();
    FDTD.DebugPEC();

    FDTD.SetGaussExcite(0.5*(f_start+f_stop), 0.5*(f_stop-f_start));
    FDTD.SetBoundaryCond({ "PML_8","PML_8","PML_8","PML_8","PML_8","PML_8" });

    // Geometry & Mesh
    _CSX CSX(new ContinuousStructure());
    FDTD.SetCSX(CSX);
    _CSRectGrid mesh(CSX.GetGrid());
    mesh.SetDeltaUnit(unit);

    std::vector<double> x_lines = {-SimBox/2, 0, SimBox/2};
    mesh.SetLines('x', x_lines);
    double resolution = C0 / f_stop / unit / 20;
    mesh.SmoothMeshLines('x', resolution);
    mesh.SetLines('y', mesh.GetLines('x'));
    mesh.SetLines('z', mesh.GetLines('x'));

    // Create a metal sphere (PEC)
    _CSProperties sphere_metal(CSX.AddMetal("sphere"));
    sphere_metal.AddSphere({0,0,0}, sphere_rad)->SetPriority(10);

    // Plane wave excitation
    double rad = inc_angle * PI/180.0;
    std::array<double,3> k_dir = {cos(rad), sin(rad), 0};
    std::array<double,3> E_dir = {0, 0, 1};

    _CSPropExcitation pw_exc = CSX.AddExcitation("plane_wave", 10, E_dir);
    pw_exc.SetPropagationDir(k_dir);
    pw_exc.SetFrequency(f0);
    std::array<double,3> start_box = {-PW_Box/2, -PW_Box/2, -PW_Box/2};
    std::array<double,3> stop_box  = { PW_Box/2,  PW_Box/2,  PW_Box/2};
    pw_exc.AddBox(start_box, stop_box);


    CSX.Write2XML("RCS_Sphere_debug.xml");

    // Create NF2FF box
    auto nf2ff = FDTD.CreateNF2FFBox();

    // Run the simulation
    FDTD.Run(Sim_Path, true, false, { "numThreads=15" });

    // Postprocessing: get excitation voltage at frequency f0
    UI_data ef({"et"}, Sim_Path, {f0});
    double Pin = 0.5 * (1.0/Z0) * std::pow(std::abs(ef.ui_f_val[0][0]), 2);

    // Calculate NF2FF (for a set of phi angles from -180 to 180 with 2° steps)
    std::vector<double> phi = linspace(-180, 180, static_cast<size_t>((180-(-180))/2 + 1));
    auto nf2ff_res = nf2ff.CalcNF2FF(Sim_Path, {f0}, {90}, phi);
    // Compute RCS for each angle and store in a vector
    std::vector<double> RCS_vals;
    for (size_t i = 0; i < phi.size(); i++) {
        double rcs_val = 4*PI/Pin * nf2ff_res.P_rad[0][0][i];
        RCS_vals.push_back(rcs_val);
    }

    // Write angle and RCS values to a CSV file
    std::ofstream csv_angle_file("RCS_vs_angle.csv");
    if (!csv_angle_file.is_open()) {
        throw std::runtime_error("Unable to open file for writing RCS_vs_angle.csv");
    }
    csv_angle_file << "Angle(deg),RCS\n";
    for (size_t i = 0; i < phi.size(); i++) {
        csv_angle_file << phi[i] << "," << RCS_vals[i] << "\n";
    }
    csv_angle_file.close();
    std::cout << "RCS vs angle data saved to RCS_vs_angle.csv" << std::endl;

    PyPlot("RCS_vs_phi", R"(
data = pd.read_csv('RCS_vs_angle.csv')

phi = np.deg2rad(data['Angle(deg)'])

plt.figure()
ax = plt.subplot(111, projection='polar')
ax.plot(phi, data['RCS'], 'k-', linewidth=2)
ax.grid(True)
ax.set_title('RCS vs Angle (Polar)')
plt.savefig('RCS_vs_phi.png')
    )");


    // // Output RCS at frequency f0
    // std::cout << "RCS at " << f0 << " Hz: " << RCS_val << std::endl;

    // Frequency sweep for RCS calculation
    std::vector<double> freq = linspace(f_start, f_stop, 100);
    UI_data ef_freq({"et"}, Sim_Path, freq);
    std::vector<double> Pin_freq;
    for (size_t i = 0; i < freq.size(); i++) {
        // Since E_dir is unit length, norm(E_dir)^2 = 1.
        double p = 0.5 * (1.0/Z0) * std::pow(std::abs(ef_freq.ui_f_val[0][i]), 2);
        Pin_freq.push_back(p);
    }
    auto nf2ff_res_freq = nf2ff.CalcNF2FF(Sim_Path, freq, {90}, {180+inc_angle}, "back_nf2ff.h5");
    std::vector<double> back_scat;
    for (size_t fn = 0; fn < freq.size(); fn++) {
        double bs = 4*M_PI/Pin_freq[fn] * nf2ff_res_freq.P_rad[fn][0][0];
        back_scat.push_back(bs);
    }

    // Save RCS vs frequency data to CSV file (for external plotting)
    {
        std::ofstream csv_file("RCS_vs_freq.csv");
        if (!csv_file.is_open()) {
            throw std::runtime_error("Unable to open file for writing RCS_vs_freq.csv");
        }
        csv_file << "Frequency(Hz),RCS\n";
        for (size_t i = 0; i < freq.size(); i++) {
            csv_file << freq[i] << "," << back_scat[i] << "\n";
        }
    }
    std::cout << "RCS vs frequency data saved to RCS_vs_freq.csv" << std::endl;

    PyPlot("RCS_vs_freq", R"(
data = pd.read_csv('RCS_vs_freq.csv')

freq = data['Frequency(Hz)']
back_scat = data['RCS']

plt.figure()
plt.plot(freq/1e6, back_scat, linewidth=2)
plt.grid()
plt.xlabel('frequency (MHz)')
plt.ylabel('RCS ($m^2$)')
plt.title('radar cross section')
plt.savefig('RCS_vs_freq.png')
    )");


    PyPlot("normalized_RCS", std::string(R"(
data = pd.read_csv('RCS_vs_freq.csv')

freq = data['Frequency(Hz)']
C0 = )") + std::to_string(299792458.0) + std::string(R"(
sphere_rad = )") + std::to_string(sphere_rad) + std::string(R"(
unit = )") + std::to_string(unit) + std::string(R"(
x = sphere_rad * unit / C0 * freq
norm_factor = np.pi * (sphere_rad * unit)**2

plt.figure()
plt.semilogy(x, data['RCS'] / norm_factor, linewidth=2)
plt.ylim([1e-2, 1e1])
plt.grid()
plt.xlabel('sphere radius / wavelength')
plt.ylabel('RCS / ($\pi a^2$)')
plt.title('normalized radar cross section')
plt.savefig('normalized_RCS.png')
    )"));

    return 0;
}
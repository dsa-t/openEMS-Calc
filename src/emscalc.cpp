#include <RCS_Sphere.h>
#include <MSL_Losses.h>

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
    // Choose one because openEMS doesn't clear m_optionDesc correctly

    //RCS_Sphere(); // Ported from RCS_Sphere.py
    MSL_Losses(); // Ported from MSL_Losses.m
}
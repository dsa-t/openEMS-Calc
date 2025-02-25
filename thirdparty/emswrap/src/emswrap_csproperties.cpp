#include <emswrap_csproperties.h>
#include <emswrap_compat.h>

#include "CSPrimPoint.h"
#include "CSPrimBox.h"
#include "CSPrimMultiBox.h"
#include "CSPrimSphere.h"
#include "CSPrimSphericalShell.h"
#include "CSPrimCylinder.h"
#include "CSPrimCylindricalShell.h"
#include "CSPrimPolygon.h"
#include "CSPrimLinPoly.h"
#include "CSPrimRotPoly.h"
#include "CSPrimPolyhedron.h"
#include "CSPrimPolyhedronReader.h"
#include "CSPrimCurve.h"
#include "CSPrimWire.h"
#include "CSPrimUserDefined.h"

#include <CSProperties.h>
#include <CSPropMaterial.h>
#include <CSPropLumpedElement.h>
#include <CSPropExcitation.h>
#include <CSPropProbeBox.h>
#include <CSPropDumpBox.h>
#include <CSPropConductingSheet.h>

#include <cassert>

// Constructors / Destructor}
_CSProperties::_CSProperties(CSProperties* ptr)
    : _ptr(ptr)
{
    if (!_ptr)
        throw std::invalid_argument("Null CSProperties pointer provided");
}

_CSProperties::~_CSProperties() {
}


// Primitive operations
int _CSProperties::GetQtyPrimitives() const {
    return _ptr->GetQtyPrimitives();
}

CSPrimitives* _CSProperties::GetPrimitive(size_t index) const {
    return _ptr->GetPrimitive(index);
}

std::vector<CSPrimitives*> _CSProperties::GetAllPrimitives() const {
    return _ptr->GetAllPrimitives();
}

// Type and identification
int _CSProperties::GetType() const {
    return _ptr->GetType();
}

std::string _CSProperties::GetTypeString() const {
    return _ptr->GetTypeString();
}

ParameterSet* _CSProperties::GetParameterSet() const {
    return _ptr->GetParameterSet();
}

double _CSProperties::GetMaterial() const {
    return _ptr->GetMaterial();
}

int _CSProperties::GetID() const {
    return _ptr->GetID();
}

// Name and color
void _CSProperties::SetName(const std::string &name) {
    _ptr->SetName(name);
}

std::string _CSProperties::GetName() const {
    return _ptr->GetName();
}

void _CSProperties::SetColor(const std::string &color, int alpha) {
    SetFillColor(color, alpha);
    SetEdgeColor(color, alpha);
}

void _CSProperties::SetFillColor(const std::string &color, int alpha) {
    unsigned int rgb = std::stoul(color.substr(color[0] == '#' ? 1 : 0), nullptr, 16);
    int r = (rgb >> 16) & 0xFF;
    int g = (rgb >> 8) & 0xFF;
    int b = rgb & 0xFF;

    _ptr->SetFillColor(r, g, b, alpha);
}

std::tuple<int,int,int,int> _CSProperties::GetFillColor() const {
    RGBa rgba = _ptr->GetFillColor();
    return std::make_tuple(rgba.R, rgba.G, rgba.B, rgba.a);
}

void _CSProperties::SetEdgeColor(const std::string &color, int alpha) {
    unsigned int rgb = std::stoul(color.substr(color[0] == '#' ? 1 : 0), nullptr, 16);
    int r = (rgb >> 16) & 0xFF;
    int g = (rgb >> 8) & 0xFF;
    int b = rgb & 0xFF;

    _ptr->SetEdgeColor(r, g, b, alpha);
    
}

std::tuple<int,int,int,int> _CSProperties::GetEdgeColor() const {
    RGBa rgba = _ptr->GetEdgeColor();
    return std::make_tuple(rgba.R, rgba.G, rgba.B, rgba.a);
}

// Visibility
bool _CSProperties::GetVisibility() const {
    return _ptr->GetVisibility();
}

void _CSProperties::SetVisibility(bool val) {
    _ptr->SetVisibility(val);
}

// Attributes
bool _CSProperties::ExistAttribute(const std::string &name) const {
    return _ptr->ExistAttribute(name);
}

std::string _CSProperties::GetAttributeValue(const std::string &name) const {
    return _ptr->GetAttributeValue(name);
}

void _CSProperties::SetAttributeValue(const std::string &name, const std::string &val) {
    _ptr->SetAttributeValue(name, val);
}

void _CSProperties::RemoveAttribute(const std::string &name) {
    _ptr->RemoveAttribute(name);
}

std::vector<std::string> _CSProperties::GetAttributeNames() const {
    return _ptr->GetAttributeNames();
}

std::map<std::string, std::string> _CSProperties::GetAttributes() const {
    std::map<std::string, std::string> attr;
    for (const auto& name : GetAttributeNames()) {
        attr[name] = GetAttributeValue(name);
    }
    return attr;
}

// Primitive creation functions
CSPrimitives* _CSProperties::AddPoint(const std::array<double, 3>& coord) {
    auto prim = new CSPrimPoint(GetParameterSet(), _ptr);
    prim->SetCoord(0, coord[0]);
    prim->SetCoord(1, coord[1]);
    prim->SetCoord(2, coord[2]);
    return prim;
}

CSPrimitives* _CSProperties::AddBox(const std::array<double, 3>& start, const std::array<double, 3>& stop) {
    // def SetStart(self, coord):
    //     """ SetStart(coord)

    //     Set the start coordinate for this box primitive.

    //     :param coord: list/array of float -- Set the start point coordinate
    //     """
    //     ptr = <_CSPrimBox*>self.thisptr
    //     assert len(coord)==3, "CSPrimBox:SetStart: length of array needs to be 3"
    //     for n in range(3):
    //         ptr.SetCoord(2*n, coord[n])

    // def SetStop(self, coord):
    //     """ SetStop(coord)

    //     Set the stop coordinate for this box primitive.

    //     :param coord: list/array of float -- Set the stop point coordinate
    //     """
    //     ptr = <_CSPrimBox*>self.thisptr
    //     assert len(coord)==3, "CSPrimBox:SetStop: length of array needs to be 3"
    //     for n in range(3):
    //         ptr.SetCoord(2*n+1, coord[n])

    assert(start.size() == 3);
    assert(stop.size() == 3);

    auto prim = new CSPrimBox(GetParameterSet(), _ptr);
    prim->SetCoord(0, start[0]);
    prim->SetCoord(2, start[1]);
    prim->SetCoord(4, start[2]);

    prim->SetCoord(1, stop[0]);
    prim->SetCoord(3, stop[1]);
    prim->SetCoord(5, stop[2]);

    return prim;
}

CSPrimitives* _CSProperties::AddCylinder(const std::array<double, 3>& start, const std::array<double, 3>& stop, double radius) {
    assert(start.size() == 3);
    assert(stop.size() == 3);

    auto prim = new CSPrimCylinder(GetParameterSet(), _ptr);
    prim->SetCoord(0, start[0]);
    prim->SetCoord(2, start[1]);
    prim->SetCoord(4, start[2]);

    prim->SetCoord(1, stop[0]);
    prim->SetCoord(3, stop[1]);
    prim->SetCoord(5, stop[2]);

    prim->SetRadius(radius);
    return prim;
}

CSPrimitives* _CSProperties::AddCylindricalShell(const std::array<double, 3>& start, const std::array<double, 3>& stop, double radius, double shell_width) {
    assert(start.size() == 3);
    assert(stop.size() == 3);

    auto prim = new CSPrimCylindricalShell(GetParameterSet(), _ptr);
    prim->SetCoord(0, start[0]);
    prim->SetCoord(2, start[1]);
    prim->SetCoord(4, start[2]);

    prim->SetCoord(1, stop[0]);
    prim->SetCoord(3, stop[1]);
    prim->SetCoord(5, stop[2]);

    prim->SetRadius(radius);
    prim->SetShellWidth(shell_width);
    return prim;
}

CSPrimitives* _CSProperties::AddSphere(const std::array<double, 3>& center, double radius) {
    assert(center.size() == 3);

    auto prim = new CSPrimSphere(GetParameterSet(), _ptr);
    prim->SetCenter(center[0], center[1], center[2]);
    prim->SetRadius(radius);
    return prim;
}

CSPrimitives* _CSProperties::AddSphericalShell(const std::array<double, 3>& center, double radius, double shell_width) {
    assert(center.size() == 3);

    auto prim = new CSPrimSphericalShell(GetParameterSet(), _ptr);
    prim->SetCenter(center[0], center[1], center[2]);
    prim->SetRadius(radius);
    prim->SetShellWidth(shell_width);
    return prim;
}

CSPrimitives *_CSProperties::AddPolygon(const std::vector<std::array<double, 2>> &points,
                                        int norm_dir, double elevation)
{
    auto prim = new CSPrimPolygon(GetParameterSet(), _ptr);

    prim->ClearCoords();
    for (const auto &p : points)
    {
        assert(p.size() == 2);
        prim->SetCoord(p[0], p[1]);
    }

    prim->SetNormDir(norm_dir);
    prim->SetElevation(elevation);
    return prim;
}

CSPrimitives *_CSProperties::AddLinPoly(const std::vector<std::array<double, 2>> &points,
                                        int norm_dir, double elevation, double length)
{
    auto prim = new CSPrimLinPoly(GetParameterSet(), _ptr);

    prim->ClearCoords();
    for (const auto &p : points)
    {
        assert(p.size() == 2);
        prim->SetCoord(p[0], p[1]);
    }

    prim->SetNormDir(norm_dir);
    prim->SetElevation(elevation);
    prim->SetLength(length);
    return prim;
}

CSPrimitives *_CSProperties::AddRotPoly(const std::vector<std::array<double, 2>> &points,
                                        int norm_dir, double elevation,
                                        int rot_axis, std::array<double, 2> angle)
{

    // """ SetAngle(a0, a1)

    // Set the start/stop angle (rad) of rotation. Default is (0, 2*pi).

    // :param a0: float -- Start angle (rad) of rotation.
    // :param a1: float -- Stop angle (rad) of rotation.
    // """
    auto prim = new CSPrimRotPoly(GetParameterSet(), _ptr);

    prim->ClearCoords();
    for (const auto &p : points)
    {
        assert(p.size() == 2);
        prim->SetCoord(p[0], p[1]);
    }

    prim->SetNormDir(norm_dir);
    prim->SetElevation(elevation);
    prim->SetRotAxisDir(rot_axis);
    prim->SetAngle(0, angle[0]);
    prim->SetAngle(1, angle[1]);
    return prim;
}

CSPrimitives* _CSProperties::AddCurve(const std::vector<std::array<double, 3>>& points) {
    auto prim = new CSPrimCurve(GetParameterSet(), _ptr);

    prim->ClearPoints();
    for (const auto &p : points)
    {
        assert(p.size() == 3);
        double point[3]{p[0], p[1], p[2]};
        prim->AddPoint(point);
    }

    return prim;
}

CSPrimitives* _CSProperties::AddWire(const std::vector<std::array<double, 3>>& points, double radius) {
    auto prim = new CSPrimWire(GetParameterSet(), _ptr);

    prim->ClearPoints();
    for (const auto &p : points)
    {
        assert(p.size() == 3);
        double point[3]{p[0], p[1], p[2]};
        prim->AddPoint(point);
    }

    prim->SetWireRadius(radius);
    return prim;
}

CSPrimitives* _CSProperties::AddPolyhedron() {
    auto prim = new CSPrimPolyhedron(GetParameterSet(), _ptr);
    return prim;
}

CSPrimitives* _CSProperties::AddPolyhedronReader(const std::string &filename) {
    auto prim = new CSPrimPolyhedronReader(GetParameterSet(), _ptr);
    prim->SetFilename(filename);
    return prim;
}


//////////////////////////////////////////////////////////////////////////
// _CSPropMaterial
//////////////////////////////////////////////////////////////////////////

_CSPropMaterial::_CSPropMaterial(CSPropMaterial* ptr)
    : _CSProperties(ptr)
{
}

void _CSPropMaterial::SetIsotropy(bool isotropic) {
    auto mat = _ptr->ToMaterial();
    if (mat) {
        mat->SetIsotropy(isotropic);
    }
}

bool _CSPropMaterial::GetIsotropy() const {
    auto mat = _ptr->ToMaterial();
    return (mat) ? mat->GetIsotropy() : false;
}

void _CSPropMaterial::SetMaterialProperty(const std::string &prop_name, double value, int ny) {
    ny = CheckNyDir(ny);
    auto mat = _ptr->ToMaterial();
    if (!mat)
        return;
    if (prop_name == "epsilon") {
        mat->SetEpsilon(value, ny);
    } else if (prop_name == "mue") {
        mat->SetMue(value, ny);
    } else if (prop_name == "kappa") {
        mat->SetKappa(value, ny);
    } else if (prop_name == "sigma") {
        mat->SetSigma(value, ny);
    } else if (prop_name == "density") {
        // Density typically does not require a directional index.
        mat->SetDensity(value);
    }
}

double _CSPropMaterial::GetMaterialProperty(const std::string &prop_name, int ny) const {
    ny = CheckNyDir(ny);
    auto mat = _ptr->ToMaterial();
    if (!mat)
        return 0.0;
    if (prop_name == "epsilon") {
        return mat->GetEpsilon(ny);
    } else if (prop_name == "mue") {
        return mat->GetMue(ny);
    } else if (prop_name == "kappa") {
        return mat->GetKappa(ny);
    } else if (prop_name == "sigma") {
        return mat->GetSigma(ny);
    } else if (prop_name == "density") {
        return mat->GetDensity();
    }
    return 0.0;
}

void _CSPropMaterial::SetMaterialWeight(const std::string &prop_name, const std::string &weightFunc, int ny) {
    ny = CheckNyDir(ny);
    auto mat = _ptr->ToMaterial();
    if (!mat)
        return;
    if (prop_name == "epsilon") {
        mat->SetEpsilonWeightFunction(weightFunc, ny);
    } else if (prop_name == "mue") {
        mat->SetMueWeightFunction(weightFunc, ny);
    } else if (prop_name == "kappa") {
        mat->SetKappaWeightFunction(weightFunc, ny);
    } else if (prop_name == "sigma") {
        mat->SetSigmaWeightFunction(weightFunc, ny);
    } else if (prop_name == "density") {
        // Density is usually a scalar without a weight function.
    }
}

std::string _CSPropMaterial::GetMaterialWeight(const std::string &prop_name, int ny) const {
    ny = CheckNyDir(ny);
    auto mat = _ptr->ToMaterial();
    if (!mat)
        return "";
    if (prop_name == "epsilon") {
        return mat->GetEpsilonWeightFunction(ny);
    } else if (prop_name == "mue") {
        return mat->GetMueWeightFunction(ny);
    } else if (prop_name == "kappa") {
        return mat->GetKappaWeightFunction(ny);
    } else if (prop_name == "sigma") {
        return mat->GetSigmaWeightFunction(ny);
    } else if (prop_name == "density") {
        return mat->GetDensityWeightFunction();
    }
    return "";
}


//////////////////////////////////////////////////////////////////////////
// _CSPropExcitation
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// _CSPropExcitation
//////////////////////////////////////////////////////////////////////////

_CSPropExcitation::_CSPropExcitation(CSPropExcitation* ptr)
    : _CSProperties(ptr)
{
}

void _CSPropExcitation::SetExcitType(int val) {
    auto exc = _ptr->ToExcitation();
    if (exc)
        exc->SetExcitType(val);
}

int _CSPropExcitation::GetExcitType() const {
    auto exc = _ptr->ToExcitation();
    return (exc ? exc->GetExcitType() : 0);
}

void _CSPropExcitation::SetEnabled(bool val) {
    auto exc = _ptr->ToExcitation();
    if (exc)
        exc->SetEnabled(val);
}

bool _CSPropExcitation::GetEnabled() const {
    auto exc = _ptr->ToExcitation();
    return (exc ? exc->GetEnabled() : false);
}

void _CSPropExcitation::SetExcitation(const std::array<double, 3>& excitation) {
    auto exc = _ptr->ToExcitation();
    if (exc)
    {
        exc->SetExcitation(excitation[0], 0);
        exc->SetExcitation(excitation[1], 1);
        exc->SetExcitation(excitation[2], 2);
    }
}

std::array<double, 3> _CSPropExcitation::GetExcitation() const {
    std::array<double, 3> vec = {0.0, 0.0, 0.0};
    auto exc = _ptr->ToExcitation();
    if (exc) {
        vec[0] = exc->GetExcitation(0);
        vec[1] = exc->GetExcitation(1);
        vec[2] = exc->GetExcitation(2);
    }
    return vec;
}

void _CSPropExcitation::SetPropagationDir(const std::array<double, 3>& propDir) {
    auto exc = _ptr->ToExcitation();
    if (exc)
    {
        exc->SetPropagationDir(propDir[0], 0);
        exc->SetPropagationDir(propDir[1], 1);
        exc->SetPropagationDir(propDir[2], 2);
    }
}

std::array<double, 3> _CSPropExcitation::GetPropagationDir() const {
    std::array<double, 3> vec = {0.0, 0.0, 0.0};
    auto exc = _ptr->ToExcitation();
    if (exc) {
        vec[0] = exc->GetPropagationDir(0);
        vec[1] = exc->GetPropagationDir(1);
        vec[2] = exc->GetPropagationDir(2);
    }
    return vec;
}

void _CSPropExcitation::SetWeightFunction(const std::array<std::string, 3>& func) {
    auto exc = _ptr->ToExcitation();
    if (exc)
    {
        exc->SetWeightFunction(func[0], 0);
        exc->SetWeightFunction(func[1], 1);
        exc->SetWeightFunction(func[2], 2);
    }
}

std::array<std::string, 3> _CSPropExcitation::GetWeightFunction() const {
    std::array<std::string, 3> funcs = {"", "", ""};
    auto exc = _ptr->ToExcitation();
    if (exc) {
        funcs[0] = exc->GetWeightFunction(0);
        funcs[1] = exc->GetWeightFunction(1);
        funcs[2] = exc->GetWeightFunction(2);
    }
    return funcs;
}

void _CSPropExcitation::SetFrequency(double val) {
    auto exc = _ptr->ToExcitation();
    if (exc)
        exc->SetFrequency(val);
}

double _CSPropExcitation::GetFrequency() const {
    auto exc = _ptr->ToExcitation();
    return (exc ? exc->GetFrequency() : 0.0);
}

void _CSPropExcitation::SetDelay(double val) {
    auto exc = _ptr->ToExcitation();
    if (exc)
        exc->SetDelay(val);
}

double _CSPropExcitation::GetDelay() const {
    auto exc = _ptr->ToExcitation();
    return (exc ? exc->GetDelay() : 0.0);
}

#pragma once

#include <string>
#include <vector>
#include <map>
#include <array>
#include <tuple>
#include <stdexcept>

// Forward declarations for dependencies
class CSProperties;
class CSPrimitives;
class CSPropMaterial;
class CSPropExcitation;
class CSPropProbeBox;
class ParameterSet;

// _CSProperties: wrapper class for a CSProperties pointer.
// It provides all the member functions exposed in the .pyx file.
class _CSProperties {
public:
    // Wrap an existing pointer; no initialization is done.
    _CSProperties(CSProperties* ptr);
    virtual ~_CSProperties();

    // Primitive operations
    int GetQtyPrimitives() const;
    CSPrimitives* GetPrimitive(size_t index) const;
    std::vector<CSPrimitives*> GetAllPrimitives() const;

    // Type and identification
    int GetType() const;
    std::string GetTypeString() const;
    ParameterSet* GetParameterSet() const;
    double GetMaterial() const; // For material properties if applicable
    int GetID() const;

    // Name and color
    void SetName(const std::string &name);
    std::string GetName() const;
    void SetColor(const std::string &color, int alpha = 255);
    void SetFillColor(const std::string &color, int alpha = 255);
    std::tuple<int,int,int,int> GetFillColor() const;
    void SetEdgeColor(const std::string &color, int alpha = 255);
    std::tuple<int,int,int,int> GetEdgeColor() const;

    // Visibility
    bool GetVisibility() const;
    void SetVisibility(bool val);

    // Attributes
    bool ExistAttribute(const std::string &name) const;
    std::string GetAttributeValue(const std::string &name) const;
    void SetAttributeValue(const std::string &name, const std::string &val);
    // Alias for SetAttributeValue
    void AddAttribute(const std::string &name, const std::string &val) { SetAttributeValue(name, val); }
    void RemoveAttribute(const std::string &name);
    std::vector<std::string> GetAttributeNames() const;
    std::map<std::string, std::string> GetAttributes() const;

    // Primitive creation functions
    CSPrimitives *AddPoint(const std::array<double, 3> &coord);
    
    CSPrimitives *AddBox(const std::array<double, 3> &start, const std::array<double, 3> &stop, int priority = 0);
    CSPrimitives *AddBox(int priority, const std::array<double, 3> &start, const std::array<double, 3> &stop);

    CSPrimitives *AddCylinder(const std::array<double, 3> &start, const std::array<double, 3> &stop, double radius);
    CSPrimitives *AddCylindricalShell(const std::array<double, 3> &start, const std::array<double, 3> &stop, double radius, double shell_width);
    CSPrimitives *AddSphere(const std::array<double, 3> &center, double radius);
    CSPrimitives *AddSphericalShell(const std::array<double, 3> &center, double radius, double shell_width);
    CSPrimitives *AddPolygon(const std::vector<std::array<double, 2>> &points, int norm_dir, double elevation);
    CSPrimitives *AddLinPoly(const std::vector<std::array<double, 2>> &points, int norm_dir, double elevation, double length);
    CSPrimitives *AddRotPoly(const std::vector<std::array<double, 2>> &points, int norm_dir, double elevation, int rot_axis, std::array<double, 2> angle);
    CSPrimitives *AddCurve(const std::vector<std::array<double, 3>> &points);
    CSPrimitives *AddWire(const std::vector<std::array<double, 3>> &points, double radius);
    CSPrimitives *AddPolyhedron();
    CSPrimitives *AddPolyhedronReader(const std::string &filename);

    CSProperties *ptr() const { return _ptr; };
    CSProperties *operator->() { return _ptr; };
    operator CSProperties *() const { return _ptr; };

protected:
    CSProperties* _ptr; // underlying C++ CSProperties pointer
};

class _CSPropMaterial : public _CSProperties
{
public:
    _CSPropMaterial(CSPropMaterial *ptr);

    void SetIsotropy(bool isotropic);
    bool GetIsotropy() const;
    void SetMaterialProperty(const std::string &prop_name, double value, int ny = 0);
    double GetMaterialProperty(const std::string &prop_name, int ny = 0) const;
    void SetMaterialWeight(const std::string &prop_name, const std::string &weightFunc, int ny = 0);
    std::string GetMaterialWeight(const std::string &prop_name, int ny = 0) const;

    CSPropMaterial *ptr() const;
    CSPropMaterial *operator->() { return ptr(); };
    operator CSPropMaterial *() const { return ptr(); };
};

class _CSPropExcitation : public _CSProperties {
public:
    _CSPropExcitation(CSPropExcitation* ptr);

    // Set and get the excitation type
    void SetExcitType(int val);
    int GetExcitType() const;

    // Enable or disable the excitation
    void SetEnabled(bool val);
    bool GetEnabled() const;

    // Set and get the excitation vector (3 elements)
    void SetExcitation(const std::array<double, 3>& excitation);
    std::array<double, 3> GetExcitation() const;

    // Set and get the propagation direction vector (3 elements)
    void SetPropagationDir(const std::array<double, 3>& propDir);
    std::array<double, 3> GetPropagationDir() const;

    // Set and get the excitation weighting functions for each axis (3 strings)
    void SetWeightFunction(const std::array<std::string, 3>& func);
    std::array<std::string, 3> GetWeightFunction() const;

    // Set and get the frequency for phase velocity compensation
    void SetFrequency(double val);
    double GetFrequency() const;

    // Set and get the delay for the excitation signal
    void SetDelay(double val);
    double GetDelay() const;

    CSPropExcitation *ptr() const;
    CSPropExcitation *operator->() { return ptr(); };
    operator CSPropExcitation *() const { return ptr(); };
};

class _CSPropProbeBox : public _CSProperties
{
public:
    _CSPropProbeBox(CSPropProbeBox *ptr);

    // Probe type accessors
    void SetProbeType(int val);
    int GetProbeType() const;

    // Weighting factor accessors
    void SetWeighting(double val);
    double GetWeighting() const;

    // Normal direction accessors
    void SetNormalDir(int val);
    int GetNormalDir() const;

    // Frequency domain sample management
    void AddFrequency(double freq);
    void ClearFrequency();
    void SetFrequency(double freq);
    std::vector<double> GetFrequency() const;
    int GetFrequencyCount() const;

    // Mode function setting, stores mode functions as attributes ("ModeFunctionX", etc.)
    void SetModeFunction(const std::array<std::string, 3> &mode_fun);

    CSPropProbeBox *ptr() const;
    CSPropProbeBox *operator->() { return ptr(); };
    operator CSPropProbeBox *() const { return ptr(); };
};

// Wrapper helpers

inline _CSProperties emswrap(CSProperties *ptr)
{
    return _CSProperties(ptr);
}

inline _CSPropMaterial emswrap(CSPropMaterial *ptr)
{
    return _CSPropMaterial(ptr);
}

inline _CSPropExcitation emswrap(CSPropExcitation *ptr)
{
    return _CSPropExcitation(ptr);
}

inline _CSPropProbeBox emswrap(CSPropProbeBox *ptr)
{
    return _CSPropProbeBox(ptr);
}
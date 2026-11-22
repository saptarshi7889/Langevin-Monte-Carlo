#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <langevin/incommensurate.hpp>
#include <langevin/commensurate.hpp>
#include <memory>
#include <mutex>
#include <variant>

namespace {
constexpr const char* CAPSULE_NAME = "langevin_monte_carlo.Simulation";
using Incommensurate = langevin::incommensurate::Simulation;
using Commensurate = langevin::commensurate::Simulation;

struct PythonError {};
struct Decref { void operator()(PyObject* value) const { Py_XDECREF(value); } };
using Object = std::unique_ptr<PyObject, Decref>;

Object owned(PyObject* value) {
    if (!value) throw PythonError{};
    return Object(value);
}

void set_item(PyObject* dictionary, const char* name, PyObject* value) {
    auto item = owned(value);
    if (PyDict_SetItemString(dictionary, name, item.get()) < 0) throw PythonError{};
}

template<class Function>
PyObject* protect(Function function) {
    try { return function(); }
    catch (const PythonError&) { return nullptr; }
    catch (const std::invalid_argument& error) { PyErr_SetString(PyExc_ValueError, error.what()); }
    catch (const std::overflow_error& error) { PyErr_SetString(PyExc_OverflowError, error.what()); }
    catch (const std::bad_alloc&) { return PyErr_NoMemory(); }
    catch (const std::exception& error) { PyErr_SetString(PyExc_RuntimeError, error.what()); }
    return nullptr;
}

// Restore the GIL even if a C++ simulation or file operation throws.
class WithoutGIL {
    PyThreadState* state_ = PyEval_SaveThread();
public:
    ~WithoutGIL() { PyEval_RestoreThread(state_); }
    WithoutGIL() = default;
    WithoutGIL(const WithoutGIL&) = delete;
    WithoutGIL& operator=(const WithoutGIL&) = delete;
};

struct Handle {
    std::variant<Incommensurate, Commensurate> simulation;
    std::mutex mutex;
    template<class Model> explicit Handle(Model model) : simulation(std::move(model)) {}
};

Handle& handle(PyObject* capsule) {
    auto* pointer = static_cast<Handle*>(PyCapsule_GetPointer(capsule, CAPSULE_NAME));
    if (!pointer) throw PythonError{};
    return *pointer;
}

std::unique_lock<std::mutex> lock(Handle& instance) {
    std::unique_lock<std::mutex> guard(instance.mutex, std::try_to_lock);
    if (!guard.owns_lock()) throw std::runtime_error("This simulation is busy in another thread");
    return guard;
}

void destroy(PyObject* capsule) {
    delete static_cast<Handle*>(PyCapsule_GetPointer(capsule, CAPSULE_NAME));
}

PyObject* capsule(std::unique_ptr<Handle> instance) {
    auto result = owned(PyCapsule_New(instance.get(), CAPSULE_NAME, destroy));
    instance.release();
    return result.release();
}

PyObject* item(PyObject* dictionary, const char* name) {
    auto* value = PyDict_GetItemString(dictionary, name);
    if (!value) {
        PyErr_SetString(PyExc_KeyError, name);
        throw PythonError{};
    }
    return value;
}

int integer(PyObject* value) {
    const long result = PyLong_AsLong(value);
    if (PyErr_Occurred()) throw PythonError{};
    if (result < std::numeric_limits<int>::min() || result > std::numeric_limits<int>::max())
        throw std::overflow_error("Integer exceeds the C++ step/size range");
    return static_cast<int>(result);
}

double real(PyObject* value) {
    const double result = PyFloat_AsDouble(value);
    if (PyErr_Occurred()) throw PythonError{};
    return result;
}

langevin::Parameters parameters(PyObject* dictionary) {
    langevin::Parameters p;
    p.L = integer(item(dictionary, "L"));
    p.Beta = integer(item(dictionary, "Beta"));
    p.burn_in = integer(item(dictionary, "burn_in"));
    p.sample_every = integer(item(dictionary, "sample_every"));
    const unsigned long seed = PyLong_AsUnsignedLong(item(dictionary, "seed"));
    if (PyErr_Occurred()) throw PythonError{};
    if (seed > std::numeric_limits<unsigned>::max()) throw std::overflow_error("Seed exceeds the C++ unsigned range");
    p.seed = static_cast<unsigned>(seed);
    p.u = real(item(dictionary, "u"));
    p.K = real(item(dictionary, "K"));
    p.alpha = real(item(dictionary, "alpha"));
    p.s = real(item(dictionary, "s"));
    p.g = real(item(dictionary, "g"));
    p.dt = real(item(dictionary, "dt"));
    return p;
}

PyObject* parameters(const langevin::Parameters& p) {
    auto result = owned(PyDict_New());
    set_item(result.get(), "L", PyLong_FromLong(p.L));
    set_item(result.get(), "Beta", PyLong_FromLong(p.Beta));
    set_item(result.get(), "burn_in", PyLong_FromLong(p.burn_in));
    set_item(result.get(), "sample_every", PyLong_FromLong(p.sample_every));
    set_item(result.get(), "seed", PyLong_FromUnsignedLong(p.seed));
    set_item(result.get(), "u", PyFloat_FromDouble(p.u));
    set_item(result.get(), "K", PyFloat_FromDouble(p.K));
    set_item(result.get(), "alpha", PyFloat_FromDouble(p.alpha));
    set_item(result.get(), "s", PyFloat_FromDouble(p.s));
    set_item(result.get(), "g", PyFloat_FromDouble(p.g));
    set_item(result.get(), "dt", PyFloat_FromDouble(p.dt));
    return result.release();
}

PyObject* buffer(const double* data, Eigen::Index count) {
    return PyBytes_FromStringAndSize(reinterpret_cast<const char*>(data),
                                   static_cast<Py_ssize_t>(count * sizeof(double)));
}

PyObject* create(PyObject*, PyObject* args) {
    int commensurate;
    PyObject* dictionary;
    if (!PyArg_ParseTuple(args, "pO!:create", &commensurate, &PyDict_Type, &dictionary)) return nullptr;
    return protect([&]() {
        const auto p = parameters(dictionary);
        std::unique_ptr<Handle> instance;
        {
            WithoutGIL release;
            if (commensurate) instance = std::make_unique<Handle>(Commensurate(p));
            else instance = std::make_unique<Handle>(Incommensurate(p));
        }
        return capsule(std::move(instance));
    });
}

PyObject* restore(PyObject*, PyObject* args) {
    int commensurate;
    const char* path;
    if (!PyArg_ParseTuple(args, "py:restore", &commensurate, &path)) return nullptr;
    return protect([&]() {
        std::unique_ptr<Handle> instance;
        {
            WithoutGIL release;
            if (commensurate) instance = std::make_unique<Handle>(Commensurate::from_checkpoint(path));
            else instance = std::make_unique<Handle>(Incommensurate::from_checkpoint(path));
        }
        return capsule(std::move(instance));
    });
}

PyObject* advance_to(PyObject*, PyObject* args) {
    PyObject* object;
    PyObject* target;
    if (!PyArg_ParseTuple(args, "OO:advance_to", &object, &target)) return nullptr;
    return protect([&]() -> PyObject* {
        auto& instance = handle(object);
        auto guard = lock(instance);
        const int final_step = integer(target);
        {
            WithoutGIL release;
            std::visit([&](auto& simulation) { simulation.advance_to(final_step); }, instance.simulation);
        }
        Py_RETURN_NONE;
    });
}

PyObject* step(PyObject*, PyObject* object) {
    return protect([&]() -> PyObject* {
        auto& instance = handle(object);
        auto guard = lock(instance);
        std::optional<langevin::Measurement> observation;
        {
            WithoutGIL release;
            observation = std::visit([](auto& simulation) { return simulation.step(); }, instance.simulation);
        }
        if (!observation) Py_RETURN_NONE;
        return Py_BuildValue("(idd)", observation->step, observation->time, observation->order_parameter);
    });
}

PyObject* save(PyObject*, PyObject* args) {
    PyObject* object;
    const char* path;
    if (!PyArg_ParseTuple(args, "Oy:save", &object, &path)) return nullptr;
    return protect([&]() -> PyObject* {
        auto& instance = handle(object);
        auto guard = lock(instance);
        {
            WithoutGIL release;
            std::visit([&](const auto& simulation) { simulation.save_checkpoint(path); }, instance.simulation);
        }
        Py_RETURN_NONE;
    });
}

PyObject* snapshot(PyObject*, PyObject* object) {
    return protect([&]() {
        auto& instance = handle(object);
        auto guard = lock(instance);
        return std::visit([](const auto& simulation) {
            const auto& state = simulation.state();
            auto result = owned(PyDict_New());
            const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> field = state.phi;
            set_item(result.get(), "parameters", parameters(state.p));
            set_item(result.get(), "commensurate", PyBool_FromLong(state.commensurate));
            set_item(result.get(), "step", PyLong_FromLong(state.step));
            set_item(result.get(), "samples", PyLong_FromLong(state.samples));
            set_item(result.get(), "order_sum", PyFloat_FromDouble(state.order_sum));
            set_item(result.get(), "phi", buffer(field.data(), field.size()));
            set_item(result.get(), "rough_x", buffer(state.rough_x.data(), state.rough_x.size()));
            set_item(result.get(), "rough_tau", buffer(state.rough_tau.data(), state.rough_tau.size()));
            return result.release();
        }, instance.simulation);
    });
}

PyObject* coefficients(PyObject*, PyObject* object) {
    return protect([&]() {
        auto& instance = handle(object);
        auto guard = lock(instance);
        return std::visit([](const auto& simulation) {
            const auto& values = simulation.coefficients();
            return buffer(values.data(), values.size());
        }, instance.simulation);
    });
}

PyObject* progress(PyObject*, PyObject* object) {
    return protect([&]() {
        auto& instance = handle(object);
        auto guard = lock(instance);
        return std::visit([](const auto& simulation) {
            const auto& state = simulation.state();
            return Py_BuildValue("(iid)", state.step, state.samples, state.step * state.p.dt);
        }, instance.simulation);
    });
}

PyMethodDef methods[] = {
    {"create", create, METH_VARARGS, "Create a model from its parameter dictionary."},
    {"restore", restore, METH_VARARGS, "Restore a model checkpoint."},
    {"advance_to", advance_to, METH_VARARGS, "Advance to an absolute step."},
    {"step", step, METH_O, "Advance one step, optionally returning a measurement."},
    {"save", save, METH_VARARGS, "Write a complete checkpoint."},
    {"snapshot", snapshot, METH_O, "Copy the field and accumulated observables."},
    {"coefficients", coefficients, METH_O, "Copy the kernel coefficients."},
    {"progress", progress, METH_O, "Return step, sample count and time."},
    {nullptr, nullptr, 0, nullptr}
};

PyModuleDef module = {
    PyModuleDef_HEAD_INIT, "_core", "Native bindings to the Langevin C++ libraries.",
    -1, methods, nullptr, nullptr, nullptr, nullptr
};
} // namespace

PyMODINIT_FUNC PyInit__core() { return PyModule_Create(&module); }

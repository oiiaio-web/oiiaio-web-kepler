#ifdef _MSC_VER
#pragma execution_character_set("utf-8")
#endif

#include <Python.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ========== 原有功能 ==========

// 计算两个物体之间的引力
static PyObject* compute_gravity(PyObject* self, PyObject* args) {
    double m1, m2, distance;
    if (!PyArg_ParseTuple(args, "ddd", &m1, &m2, &distance)) return NULL;

    if (distance <= 0) {
        PyErr_SetString(PyExc_ValueError, "distance must be positive");
        return NULL;
    }

    double G = 6.67430e-11;
    double force = G * m1 * m2 / (distance * distance);
    
    return Py_BuildValue("d", force);
}

// 计算轨道的“速度”
static PyObject* orbital_velocity(PyObject* self, PyObject* args) {
    double mass, radius;
    if (!PyArg_ParseTuple(args, "dd", &mass, &radius)) return NULL;

    if (radius <= 0) {
        PyErr_SetString(PyExc_ValueError, "radius must be positive");
        return NULL;
    }

    double G = 6.67430e-11;
    double v = sqrt(G * mass / radius);
    
    return Py_BuildValue("d", v);
}

// 霍曼转移轨道计算
static PyObject* hohmann_transfer(PyObject* self, PyObject* args) {
    double r1, r2, mu;
    if (!PyArg_ParseTuple(args, "ddd", &r1, &r2, &mu)) return NULL;

    if (r1 <= 0 || r2 <= 0 || mu <= 0) {
        PyErr_SetString(PyExc_ValueError, "r1, r2, mu must be positive");
        return NULL;
    }

    double a_transfer = (r1 + r2) / 2.0;
    double v1 = sqrt(mu * (2.0 / r1 - 1.0 / a_transfer));
    double v2 = sqrt(mu * (2.0 / r2 - 1.0 / a_transfer));
    double v_circular1 = sqrt(mu / r1);
    double v_circular2 = sqrt(mu / r2);
    double delta_v1 = v1 - v_circular1;
    double delta_v2 = v_circular2 - v2;
    double total_delta_v = fabs(delta_v1) + fabs(delta_v2);
    double transfer_time = M_PI * sqrt(a_transfer * a_transfer * a_transfer / mu);
    
    return Py_BuildValue("(dddd)", delta_v1, delta_v2, total_delta_v, transfer_time);
}

// ========== 1.4.0 新增功能 ==========

// 1. 开普勒第三定律验证（带误差显示）
static PyObject* kepler_third_law(PyObject* self, PyObject* args) {
    double a1, T1, a2, T2;
    if (!PyArg_ParseTuple(args, "dddd", &a1, &T1, &a2, &T2)) return NULL;

    if (a1 <= 0 || a2 <= 0 || T1 <= 0 || T2 <= 0) {
        PyErr_SetString(PyExc_ValueError, "all parameters must be positive");
        return NULL;
    }

    double left = (T1 * T1) / (T2 * T2);
    double right = (a1 * a1 * a1) / (a2 * a2 * a2);
    double diff = fabs(left - right);

    if (diff < 1e-6) {
        char result[256];
        snprintf(result, sizeof(result),
            "✅ 开普勒第三定律验证通过！误差: %.2e | ✅ Kepler's third law verified! Error: %.2e",
            diff, diff);
        return Py_BuildValue("s", result);
    } else {
        char result[256];
        snprintf(result, sizeof(result),
            "❌ 开普勒第三定律验证失败！误差: %.2e 过大 | ❌ Kepler's third law verification failed! Error: %.2e too large",
            diff, diff);
        return Py_BuildValue("s", result);
    }
}

// 2. 轨道离心率
static PyObject* orbital_eccentricity(PyObject* self, PyObject* args) {
    double v, r, mu;
    if (!PyArg_ParseTuple(args, "ddd", &v, &r, &mu)) return NULL;

    if (r <= 0 || mu <= 0) {
        PyErr_SetString(PyExc_ValueError, "r and mu must be positive");
        return NULL;
    }

    double epsilon = (v * v) / 2.0 - mu / r;
    double h = r * v;
    double inner = 1.0 + (2.0 * epsilon * h * h) / (mu * mu);

    if (inner < 0) {
        PyErr_SetString(PyExc_ValueError, "Invalid orbit: negative under sqrt");
        return NULL;
    }

    return Py_BuildValue("d", sqrt(inner));
}

// 3. 最近/最远点距离
static PyObject* apsis_distances(PyObject* self, PyObject* args) {
    double a, e;
    if (!PyArg_ParseTuple(args, "dd", &a, &e)) return NULL;

    if (a <= 0) {
        PyErr_SetString(PyExc_ValueError, "a must be positive");
        return NULL;
    }

    if (e >= 1.0) {
        PyErr_SetString(PyExc_ValueError, "e >= 1: not a closed orbit");
        return NULL;
    }

    double r_p = a * (1.0 - e);
    double r_a = a * (1.0 + e);

    return Py_BuildValue("(dd)", r_p, r_a);
}

// 4. 拉格朗日点（L1~L5）
static PyObject* lagrange_point(PyObject* self, PyObject* args) {
    double m1, m2, r;
    int point;
    if (!PyArg_ParseTuple(args, "dddi", &m1, &m2, &r, &point)) return NULL;

    if (m1 <= 0 || m2 <= 0 || r <= 0) {
        PyErr_SetString(PyExc_ValueError, "m1, m2, r must be positive");
        return NULL;
    }

    if (point < 1 || point > 5) {
        PyErr_SetString(PyExc_ValueError, "point must be 1, 2, 3, 4, or 5");
        return NULL;
    }

    double mu = m2 / (m1 + m2);

    // L4 和 L5：等边三角形顶点
    if (point == 4 || point == 5) {
        double x = r * (0.5 - mu);
        double y = r * sqrt(3.0) / 2.0;
        if (point == 5) y = -y;
        return Py_BuildValue("(dd)", x, y);
    }

    // L1、L2、L3：牛顿迭代
    double gamma = pow(mu / 3.0, 1.0 / 3.0);
    double x = 0.0;

    if (point == 1) {
        x = r * (1.0 - gamma);
    } else if (point == 2) {
        x = r * (1.0 + gamma);
    } else {
        x = -r * (1.0 + 5.0 * mu / 12.0);
    }

    for (int i = 0; i < 50; i++) {
        double f, df;
        if (point == 1) {
            double r1 = x;
            double r2 = r - x;
            f = m1 / (r1 * r1) - m2 / (r2 * r2) - x * (m1 + m2) / (r * r * r);
            df = -2.0 * m1 / (r1 * r1 * r1) - 2.0 * m2 / (r2 * r2 * r2) - (m1 + m2) / (r * r * r);
        } else if (point == 2) {
            double r1 = x;
            double r2 = x - r;
            f = m1 / (r1 * r1) + m2 / (r2 * r2) - x * (m1 + m2) / (r * r * r);
            df = -2.0 * m1 / (r1 * r1 * r1) - 2.0 * m2 / (r2 * r2 * r2) - (m1 + m2) / (r * r * r);
        } else {
            double r1 = -x;
            double r2 = r - x;
            f = -m1 / (r1 * r1) - m2 / (r2 * r2) - x * (m1 + m2) / (r * r * r);
            df = -2.0 * m1 / (r1 * r1 * r1) - 2.0 * m2 / (r2 * r2 * r2) - (m1 + m2) / (r * r * r);
        }

        if (fabs(df) < 1e-12) break;
        double dx = f / df;
        x -= dx;
        if (fabs(dx) < 1e-10) break;
    }

    return Py_BuildValue("(dd)", x, 0.0);
}

// 5. 多体引力合力
static PyObject* multi_body_gravity(PyObject* self, PyObject* args) {
    PyObject* bodies_list;
    if (!PyArg_ParseTuple(args, "O", &bodies_list)) return NULL;

    if (!PyList_Check(bodies_list)) {
        PyErr_SetString(PyExc_TypeError, "Expected a list of bodies");
        return NULL;
    }

    Py_ssize_t n = PyList_Size(bodies_list);
    if (n < 2) {
        return Py_BuildValue("(dd)", 0.0, 0.0);
    }

    double G = 6.67430e-11;
    double fx = 0.0, fy = 0.0;

    for (Py_ssize_t i = 0; i < n; i++) {
        for (Py_ssize_t j = i + 1; j < n; j++) {
            PyObject* body_i = PyList_GetItem(bodies_list, i);
            PyObject* body_j = PyList_GetItem(bodies_list, j);

            double mi = PyFloat_AsDouble(PyTuple_GetItem(body_i, 0));
            double xi = PyFloat_AsDouble(PyTuple_GetItem(body_i, 1));
            double yi = PyFloat_AsDouble(PyTuple_GetItem(body_i, 2));

            double mj = PyFloat_AsDouble(PyTuple_GetItem(body_j, 0));
            double xj = PyFloat_AsDouble(PyTuple_GetItem(body_j, 1));
            double yj = PyFloat_AsDouble(PyTuple_GetItem(body_j, 2));

            double dx = xj - xi;
            double dy = yj - yi;
            double dist2 = dx * dx + dy * dy;
            double dist = sqrt(dist2);

            if (dist < 1e-6) continue;

            double F = G * mi * mj / dist2;
            fx += F * dx / dist;
            fy += F * dy / dist;
        }
    }

    return Py_BuildValue("(dd)", fx, fy);
}

// ========== 版本更新信息 ==========

const char* get_update_info(const char* version) {
    if (strcmp(version, "1.4.1") == 0) {
        return "【中文】1.4.1 版本：新增 macOS Intel、macOS Apple Silicon、Linux aarch64 三个平台的预编译包，全平台制霸。【English】Version 1.4.1: Added prebuilt wheels for macOS Intel, macOS Apple Silicon, and Linux aarch64. Now available on all platforms.";
    }
    if (strcmp(version, "1.4.0") == 0) {
        return "🚀 史前大更新 | 1.4.0 版本：一夜之间新增开普勒第三定律、轨道离心率、近远点距离、拉格朗日点、多体引力合力五大物理模块！从此升级为航天模拟器引擎！ | 🚀 Prehistoric Mega Update | Version 1.4.0: Added 5 physics modules overnight — Kepler's third law, orbital eccentricity, apsis distances, Lagrange points, and multi-body gravity! Now a spaceflight simulator engine!";
    }
    if (strcmp(version, "1.3.2") == 0) {
        return "【中文】1.3.2 版本：新增 macOS x86_64、macOS arm64、Linux aarch64 平台的预编译包，修正 Homepage 链接。【English】Version 1.3.2: Added prebuilt wheels for macOS x86_64, macOS arm64, and Linux aarch64, fixed Homepage link.";
    }
    if (strcmp(version, "1.3.1") == 0) {
        return "【中文】1.3.1 版本：添加 python_requires、classifiers、project_urls，优化 setup.py 支持 Linux 和 Windows 双平台编译。【English】Version 1.3.1: Added python_requires, classifiers, project_urls, and improved setup.py to support both Linux and Windows compilation.";
    }
    if (strcmp(version, "1.3.0") == 0) {
        return "【中文】1.3.0 版本：新增霍曼转移轨道计算函数 hohmann_transfer()，并提供别名 hohmann()。【English】Version 1.3.0: Added hohmann_transfer() function for Hohmann transfer orbit calculation, with alias hohmann().";
    }
    if (strcmp(version, "1.2.2") == 0) {
        return "【中文】1.2.2 版本：移除了 play_desktop() 功能，不再依赖 pygame-ce，纯轨道计算器。【English】Version 1.2.2: Removed play_desktop() function, no longer depends on pygame-ce, pure orbital mechanics calculator.";
    }
    if (strcmp(version, "1.2.1") == 0) {
        return "【中文】1.2.1 版本：更新了 abi3 稳定接口，明确支持 Python 3.8~3.14 版本。【English】Version 1.2.1: Updated abi3 stable interface, clearly supports Python 3.8~3.14 versions.";
    }
    if (strcmp(version, "1.2.0.post1") == 0) {
        return "【中文】1.2.0.post1 版本：修复了 PyPI 页面显示更新日志的问题。【English】Version 1.2.0.post1: Fixed the issue of displaying changelog on the PyPI page.";
    }
    if (strcmp(version, "1.2.0") == 0) {
        return "【中文】1.2.0 版本：正式发布，修复了部分用户无法安装的 bug。【English】Version 1.2.0: Officially released, fixed installation bugs for some users.";
    }
    if (strcmp(version, "1.2.0rc1") == 0) {
        return "【中文】1.2.0rc1 版本：新增游戏模块，可打开自制桌面系统。【English】Version 1.2.0rc1: Added game module, can open custom desktop system.";
    }
    if (strcmp(version, "1.2.0b1") == 0) {
        return "【中文】1.2.0b1 版本：Beta 测试版，使用 abi3 稳定接口，支持所有 Python 3.8+ 版本。【English】Version 1.2.0b1: Beta version, uses abi3 stable interface, supports all Python 3.8+ versions.";
    }
    if (strcmp(version, "1.2.0a2.post3") == 0) {
        return "【中文】1.2.0a2.post3 版本：修复必要 bug。【English】Version 1.2.0a2.post3: Fixed necessary bugs.";
    }
    if (strcmp(version, "1.2.0a2.post2") == 0) {
        return "【中文】1.2.0a2.post2 版本：使用 MSVC 编译器重新编译，修复 DLL 依赖问题。【English】Version 1.2.0a2.post2: Recompiled with MSVC compiler, fixed DLL dependency issues.";
    }
    if (strcmp(version, "1.2.0a2.post1") == 0) {
        return "【中文】1.2.0a2.post1 版本：使用 abi3 稳定接口，支持所有 Python 3.8+ 版本。【English】Version 1.2.0a2.post1: Uses abi3 stable interface, supports all Python 3.8+ versions.";
    }
    if (strcmp(version, "1.2.0a2") == 0) {
        return "【中文】1.2.0a2 版本：在 1.2.0a1 的基础上，增加了 abi3 稳定接口，支持所有 Python 3.8+ 版本。【English】Version 1.2.0a2: Added abi3 stable interface on top of 1.2.0a1, supports all Python 3.8+ versions.";
    }
    if (strcmp(version, "1.2.0a1") == 0) {
        return "【中文】1.2.0a1 版本：首次使用 abi3 稳定接口，支持所有 Python 3.8+ 版本。【English】Version 1.2.0a1: First use of abi3 stable interface, supports all Python 3.8+ versions.";
    }
    if (strcmp(version, "1.1.5.post2") == 0) {
        return "【中文】1.1.5.post2 版本：更新了 Apache License 2.0 许可证中的版权信息。【English】Version 1.1.5.post2: Updated copyright information in the Apache License 2.0.";
    }
    if (strcmp(version, "1.1.5.post1") == 0) {
        return "【中文】1.1.5.post1 版本：修正了 README 中 1.0.0 和 1.0.1 的发布日期。【English】Version 1.1.5.post1: Fixed the release dates of versions 1.0.0 and 1.0.1 in the README.";
    }
    if (strcmp(version, "1.1.5") == 0) {
        return "【中文】1.1.5 版本：修复了网页显示的 bug。【English】Version 1.1.5: Fixed the webpage display bug.";
    }
    if (strcmp(version, "1.1.4") == 0) {
        return "【中文】1.1.4 版本：升级到新版本，新增功能。【English】Version 1.1.4: Upgraded to new version, added new features.";
    }
    if (strcmp(version, "1.1.3") == 0) {
        return "【中文】1.1.3 版本：修复了版本更新信息查询功能的编码错误。【English】Version 1.1.3: Fixed encoding error in version update info query function.";
    }
    if (strcmp(version, "1.1.2") == 0) {
        return "【中文】1.1.2 版本：完善了项目描述，添加了安装和使用说明。【English】Version 1.1.2: Improved project description, added installation and usage instructions.";
    }
    if (strcmp(version, "1.1.1") == 0) {
        return "【中文】1.1.1 版本：新增了 update_information 功能，更新了许可证为 Apache 2.0。【English】Version 1.1.1: Added update_information function, updated license to Apache 2.0.";
    }
    if (strcmp(version, "1.1.0") == 0) {
        return "【中文】1.1.0 版本：新增了 update_information 功能。【English】Version 1.1.0: Added update_information function.";
    }
    if (strcmp(version, "1.0.1") == 0) {
        return "【中文】1.0.1 版本：更新了 README 和许可证。【English】Version 1.0.1: Updated README and license.";
    }
    if (strcmp(version, "1.0.0") == 0) {
        return "【中文】1.0.0 版本：首次发布，支持计算引力和轨道速度。【English】Version 1.0.0: Initial release with gravity and orbital velocity calculations.";
    }
    return NULL;
}

static PyObject* update_information(PyObject* self, PyObject* args) {
    const char* version;
    if (!PyArg_ParseTuple(args, "s", &version)) return NULL;

    const char* info = get_update_info(version);
    if (info == NULL) {
        PyErr_SetString(PyExc_ValueError, "版本号不存在！");
        return NULL;
    }

    return Py_BuildValue("s", info);
}

// ========== 方法表 ==========

static PyMethodDef KeplerMethods[] = {
    {"update_information", update_information, METH_VARARGS, "查询版本更新信息"},
    {"compute_gravity", compute_gravity, METH_VARARGS, "计算引力"},
    {"orbital_velocity", orbital_velocity, METH_VARARGS, "计算轨道速度"},
    {"hohmann_transfer", hohmann_transfer, METH_VARARGS, "霍曼转移轨道计算"},
    {"hohmann", hohmann_transfer, METH_VARARGS, "霍曼转移轨道计算（别名）"},
    {"kepler_third_law", kepler_third_law, METH_VARARGS, "开普勒第三定律验证"},
    {"orbital_eccentricity", orbital_eccentricity, METH_VARARGS, "轨道离心率"},
    {"apsis_distances", apsis_distances, METH_VARARGS, "最近/最远点距离"},
    {"lagrange_point", lagrange_point, METH_VARARGS, "拉格朗日点计算（point=1~5）"},
    {"multi_body_gravity", multi_body_gravity, METH_VARARGS, "多体引力合力计算"},
    {NULL, NULL, 0, NULL}
};

// ========== 模块定义 ==========

static struct PyModuleDef kepler_module = {
    PyModuleDef_HEAD_INIT,
    "oiiaio_web_kepler",
    NULL,
    -1,
    KeplerMethods
};

// ========== 初始化 ==========

PyMODINIT_FUNC PyInit_oiiaio_web_kepler(void) {
    return PyModule_Create(&kepler_module);
}
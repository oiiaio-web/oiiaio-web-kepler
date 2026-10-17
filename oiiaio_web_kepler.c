#ifdef _MSC_VER
#pragma execution_character_set("utf-8")
#endif

#include <Python.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165271201909145648566923460348610454326648213393607260249141273724587006606315588174881520920962829254091715364367892590360011330530548820466521384146951941511609433057270365759591953092186117381932611793105118548074462379962749567351885752724891227938183011949128831426076900422421902267105562632111110937054421750694165896040807198403850962455444362981230987879927244284909188845801561660979191338754992005240636899125607176060588611646710940507754100225698315520005593572972571636269561882670428252483600823257530420752963450
#endif

// ========== 星球数据库 / Planet database ==========

typedef struct {
    const char* name_cn;
    const char* name_en;
    double radius;            // 米 / meters
    double mass;              // 千克 / kilograms
    double gravity;           // m/s^2
    double escape_velocity;   // 米/秒 / m/s
    double orbital_period;    // 秒（绕太阳公转）/ seconds (around Sun)
    double rotation_period;   // 秒（自转）/ seconds (self rotation)
    double atmosphere_height; // 米，0 表示无大气 / meters, 0 = no atmosphere
    double ocean_depth;       // 米，0 表示无海洋 / meters, 0 = no ocean
} PlanetInfo;

static const PlanetInfo PLANETS[] = {
    {"水星", "Mercury", 2439700, 3.301e23, 3.7, 4250, 7600521, 5067000, 0, 0},
    {"金星", "Venus", 6051800, 4.867e24, 8.87, 10360, 19414146, 20997000, 250000, 0},
    {"地球", "Earth", 6371000, 5.972e24, 9.8, 11186, 31558149, 86164, 100000, 5000},
    {"火星", "Mars", 3389500, 6.417e23, 3.71, 5030, 59354294, 88642, 100000, 0},
    {"木星", "Jupiter", 69911000, 1.898e27, 24.79, 59500, 374335776, 35730, 500000, 0},
    {"土星", "Saturn", 58232000, 5.683e26, 10.44, 35500, 929596608, 38360, 600000, 0},
    {"天王星", "Uranus", 25362000, 8.681e25, 8.87, 21300, 2651486400, 62064, 300000, 0},
    {"海王星", "Neptune", 24622000, 1.024e26, 11.15, 23500, 5197910400, 57996, 300000, 0},
    {"月球", "Moon", 1737400, 7.342e22, 1.62, 2380, 2360591, 2360591, 0, 0},
    {"太阳", "Sun", 696340000, 1.989e30, 274, 617700, 0, 2192832, 0, 0},
};

#define PLANET_COUNT (sizeof(PLANETS) / sizeof(PLANETS[0]))

// 根据名字找星球质量与半径。成功返回 1，失败返回 0
static int get_planet_data(const char* planet_name, double* out_mass, double* out_radius) {
    for (size_t i = 0; i < PLANET_COUNT; i++) {
        if (strcmp(planet_name, PLANETS[i].name_cn) == 0 ||
            strcmp(planet_name, PLANETS[i].name_en) == 0) {
            *out_mass = PLANETS[i].mass;
            *out_radius = PLANETS[i].radius;
            return 1;
        }
    }
    return 0;
}

// ========== 基础物理 / Basic physics ==========

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

// ========== 霍曼转移 / Hohmann transfer ==========

static PyObject* hohmann_transfer(PyObject* self, PyObject* args) {
    double r1 = 0, r2 = 0, mu = 0;

    const char* planet1;
    const char* planet2;

    // 模式一：两个星球名 / Mode 1: two planet names
    if (PyArg_ParseTuple(args, "ss", &planet1, &planet2)) {
        double mass1, radius1, mass2, radius2;
        if (!get_planet_data(planet1, &mass1, &radius1)) {
            PyErr_SetString(PyExc_ValueError, "未知的第一颗星球！ / Unknown first planet!");
            return NULL;
        }
        if (!get_planet_data(planet2, &mass2, &radius2)) {
            PyErr_SetString(PyExc_ValueError, "未知的第二颗星球！ / Unknown second planet!");
            return NULL;
        }
        // 以第一颗星球为中心天体 / Use first planet as central body
        r1 = radius1;
        r2 = radius2;
        mu = 6.67430e-11 * mass1;
    }
    // 模式二：三个数字 / Mode 2: three numbers
    else if (PyArg_ParseTuple(args, "ddd", &r1, &r2, &mu)) {
        // 数字模式，直接使用 / numeric mode
    }
    else {
        PyErr_SetString(PyExc_TypeError,
            "hohmann_transfer 支持两种调用方式 / supports two forms:\n"
            "  1) hohmann_transfer(r1_m, r2_m, mu_m3_s2)\n"
            "  2) hohmann_transfer('Earth', 'Mars')");
        return NULL;
    }

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

// ========== 轨道分析 / Orbital analysis ==========

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
    char result[256];
    if (diff < 1e-6) {
        snprintf(result, sizeof(result),
            "✅ 开普勒第三定律验证通过！误差: %.2e | ✅ Kepler's third law verified! Error: %.2e",
            diff, diff);
    } else {
        snprintf(result, sizeof(result),
            "❌ 开普勒第三定律验证失败！误差: %.2e 过大 | ❌ Kepler's third law verification failed! Error: %.2e too large",
            diff, diff);
    }
    return Py_BuildValue("s", result);
}

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
    if (point == 4 || point == 5) {
        double x = r * (0.5 - mu);
        double y = r * sqrt(3.0) / 2.0;
        if (point == 5) y = -y;
        return Py_BuildValue("(dd)", x, y);
    }
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

// ========== 星球参数查询 / Planet info ==========

static PyObject* planet_information(PyObject* self, PyObject* args) {
    const char* planet_name;
    if (!PyArg_ParseTuple(args, "s", &planet_name)) return NULL;

    for (size_t i = 0; i < PLANET_COUNT; i++) {
        if (strcmp(planet_name, PLANETS[i].name_cn) == 0 ||
            strcmp(planet_name, PLANETS[i].name_en) == 0) {
            char result[2048];
            snprintf(result, sizeof(result),
                "【中文】星球：%s\n"
                "  半径：%.0f 米\n"
                "  质量：%.3e 千克\n"
                "  表面重力：%.2f m/s²\n"
                "  逃逸速度：%.0f m/s\n"
                "  公转周期：%.0f 秒\n"
                "  自转周期：%.0f 秒\n"
                "  大气层高度：%.0f 米\n"
                "  海洋深度：%.0f 米\n"
                "【English】Planet: %s\n"
                "  Radius: %.0f m\n"
                "  Mass: %.3e kg\n"
                "  Surface Gravity: %.2f m/s²\n"
                "  Escape Velocity: %.0f m/s\n"
                "  Orbital Period: %.0f s\n"
                "  Rotation Period: %.0f s\n"
                "  Atmosphere Height: %.0f m\n"
                "  Ocean Depth: %.0f m",
                PLANETS[i].name_cn,
                PLANETS[i].radius,
                PLANETS[i].mass,
                PLANETS[i].gravity,
                PLANETS[i].escape_velocity,
                PLANETS[i].orbital_period,
                PLANETS[i].rotation_period,
                PLANETS[i].atmosphere_height,
                PLANETS[i].ocean_depth,
                PLANETS[i].name_en,
                PLANETS[i].radius,
                PLANETS[i].mass,
                PLANETS[i].gravity,
                PLANETS[i].escape_velocity,
                PLANETS[i].orbital_period,
                PLANETS[i].rotation_period,
                PLANETS[i].atmosphere_height,
                PLANETS[i].ocean_depth);
            return Py_BuildValue("s", result);
        }
    }

    PyErr_SetString(PyExc_ValueError,
        "未知星球！支持：水星、金星、地球、火星、木星、土星、天王星、海王星、月球、太阳\n"
        "Unknown planet! Supported: Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Moon, Sun");
    return NULL;
}

// ========== 版本更新信息 / Version info ==========

const char* get_update_info(const char* version) {
    if (strcmp(version, "1.6.0") == 0) {
        return "【中文】1.6.0 版本：hohmann_transfer 支持直接输入星球名称（如 'Earth', 'Mars'），自动计算星际转移轨道。【English】Version 1.6.0: hohmann_transfer now accepts planet names (e.g., 'Earth', 'Mars') for automatic interplanetary transfer calculations.";
    }
    if (strcmp(version, "1.5.0") == 0) {
        return "【中文】1.5.0 版本：新增 planet_information() / ptmation() 函数，内置 10 个星球参数。【English】Version 1.5.0: Added planet_information() / ptmation() functions with built-in data for 10 planets.";
    }
    if (strcmp(version, "1.4.2") == 0) {
        return "【中文】1.4.2 版本：新增 Python 3.15 兼容。【English】Version 1.4.2: Added Python 3.15 compatibility.";
    }
    if (strcmp(version, "1.4.1") == 0) {
        return "【中文】1.4.1 版本：新增 macOS、Linux aarch64、musllinux 预编译包，全平台制霸。【English】Version 1.4.1: Added macOS, Linux aarch64, and musllinux wheels. Now available on all platforms.";
    }
    if (strcmp(version, "1.4.0") == 0) {
        return "【中文】1.4.0 版本：新增开普勒第三定律、轨道离心率、近远点距离、拉格朗日点、多体引力合力五大物理模块。【English】Version 1.4.0: Added 5 physics modules — Kepler's third law, orbital eccentricity, apsis distances, Lagrange points, and multi-body gravity.";
    }
    if (strcmp(version, "1.3.2") == 0) {
        return "【中文】1.3.2 版本：新增 macOS x86_64、macOS arm64、Linux aarch64 平台的预编译包。【English】Version 1.3.2: Added prebuilt wheels for macOS x86_64, macOS arm64, and Linux aarch64.";
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
    if (strcmp(version, "1.2.0") == 0) {
        return "【中文】1.2.0 版本：正式发布，修复了部分用户无法安装的 bug。【English】Version 1.2.0: Officially released, fixed installation bugs for some users.";
    }
    if (strcmp(version, "1.1.5") == 0) {
        return "【中文】1.1.5 版本：修复了网页显示的 bug。【English】Version 1.1.5: Fixed the webpage display bug.";
    }
    if (strcmp(version, "1.1.1") == 0) {
        return "【中文】1.1.1 版本：新增了 update_information 功能，更新了许可证为 Apache 2.0。【English】Version 1.1.1: Added update_information function, updated license to Apache 2.0.";
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
        PyErr_SetString(PyExc_ValueError, "版本号不存在！ / Version not found!");
        return NULL;
    }
    return Py_BuildValue("s", info);
}

// ========== 方法表 / Method table (中英双语 docstring) ==========

static PyMethodDef KeplerMethods[] = {
    {
        "update_information",
        update_information,
        METH_VARARGS,
        "update_information(version) -> str\n\n"
        "查询指定版本的更新说明（中英文双语）。\n"
        "Query the release notes for a specific version (bilingual CN/EN).\n\n"
        "参数 / Parameters:\n"
        "    version : str，版本号字符串 / version string, e.g. '1.6.0', '1.5.0', '1.0.0'\n\n"
        "返回 / Returns:\n"
        "    str，该版本的更新说明 / release notes；版本不存在时抛出 ValueError。"
    },
    {
        "compute_gravity",
        compute_gravity,
        METH_VARARGS,
        "compute_gravity(m1, m2, distance) -> float\n\n"
        "计算两个物体之间的万有引力（牛顿）。\n"
        "Calculate the gravitational force between two bodies (Newtons).\n\n"
        "参数 / Parameters:\n"
        "    m1       : float，物体1质量，单位：千克(kg) / mass of body 1, kg\n"
        "    m2       : float，物体2质量，单位：千克(kg) / mass of body 2, kg\n"
        "    distance : float，两物体中心距离，单位：米(m)，必须 > 0\n"
        "               distance between centers, meters, must be > 0\n\n"
        "返回 / Returns:\n"
        "    float，引力大小，单位：牛顿(N) / gravitational force, Newtons\n\n"
        "示例 / Example:\n"
        "    compute_gravity(5.972e24, 1000, 6371000)  # 地球对1000kg物体的引力"
    },
    {
        "orbital_velocity",
        orbital_velocity,
        METH_VARARGS,
        "orbital_velocity(mass, radius) -> float\n\n"
        "计算圆轨道上的环绕速度。\n"
        "Calculate circular orbital velocity.\n\n"
        "参数 / Parameters:\n"
        "    mass   : float，中心天体质量，单位：千克(kg) / central body mass, kg\n"
        "    radius : float，轨道半径（从中心天体中心算起），单位：米(m)，必须 > 0\n"
        "             orbital radius from center, meters, must be > 0\n\n"
        "返回 / Returns:\n"
        "    float，环绕速度，单位：米/秒(m/s) / orbital velocity, m/s\n\n"
        "示例 / Example:\n"
        "    orbital_velocity(5.972e24, 6800000)  # 地球低轨道环绕速度"
    },
    {
        "hohmann_transfer",
        hohmann_transfer,
        METH_VARARGS,
        "hohmann_transfer(...) -> tuple\n\n"
        "霍曼转移轨道计算。支持两种调用方式。\n"
        "Hohmann transfer orbit calculation. Supports two calling forms.\n\n"
        "方式1（数字）/ Form 1 (numeric): hohmann_transfer(r1, r2, mu)\n"
        "    r1 : float，起始圆轨道半径，单位：米(m)，必须 > 0\n"
        "         initial circular orbit radius, meters, must be > 0\n"
        "    r2 : float，目标圆轨道半径，单位：米(m)，必须 > 0\n"
        "         target circular orbit radius, meters, must be > 0\n"
        "    mu : float，中心天体引力常数 G*M，单位：m³/s²，必须 > 0\n"
        "         gravitational parameter G*M, m^3/s^2, must be > 0\n\n"
        "方式2（星球名）/ Form 2 (planet names): hohmann_transfer('Earth', 'Mars')\n"
        "    以第一颗星球为中心天体，使用其质量计算 mu，半径取两颗星球各自的半径。\n"
        "    支持中英文名。\n"
        "    Uses the first planet as central body; radius = each planet's radius.\n"
        "    Supports both Chinese and English names.\n\n"
        "返回 / Returns:\n"
        "    tuple(delta_v1, delta_v2, total_delta_v, transfer_time)\n"
        "    delta_v1      : float，第一次点火速度增量，单位：m/s / first burn delta-v, m/s\n"
        "    delta_v2      : float，第二次点火速度增量，单位：m/s / second burn delta-v, m/s\n"
        "    total_delta_v : float，总速度增量，单位：m/s / total delta-v, m/s\n"
        "    transfer_time : float，转移飞行时间，单位：秒(s) / transfer time, seconds\n\n"
        "示例 / Example:\n"
        "    hohmann_transfer(6371000, 42164000, 3.986e14)\n"
        "    hohmann_transfer('Earth', 'Mars')"
    },
    {
        "hohmann",
        hohmann_transfer,
        METH_VARARGS,
        "hohmann(...) -> tuple\n\n"
        "hohmann_transfer 的别名，参数和返回值完全相同。\n"
        "Alias of hohmann_transfer; same parameters and return value."
    },
    {
        "kepler_third_law",
        kepler_third_law,
        METH_VARARGS,
        "kepler_third_law(a1, T1, a2, T2) -> str\n\n"
        "验证开普勒第三定律：T² ∝ a³。\n"
        "Verify Kepler's third law: T^2 proportional to a^3.\n\n"
        "参数 / Parameters:\n"
        "    a1 : float，天体1半长轴，单位：米(m)，必须 > 0 / semi-major axis 1, meters\n"
        "    T1 : float，天体1周期，单位：秒(s)，必须 > 0 / period 1, seconds\n"
        "    a2 : float，天体2半长轴，单位：米(m)，必须 > 0 / semi-major axis 2, meters\n"
        "    T2 : float，天体2周期，单位：秒(s)，必须 > 0 / period 2, seconds\n\n"
        "返回 / Returns:\n"
        "    str，验证结果及误差（中英文双语）/ verification result and error (bilingual)"
    },
    {
        "orbital_eccentricity",
        orbital_eccentricity,
        METH_VARARGS,
        "orbital_eccentricity(v, r, mu) -> float\n\n"
        "由速度、半径、引力常数计算轨道离心率。\n"
        "Calculate orbital eccentricity from velocity, radius and mu.\n\n"
        "参数 / Parameters:\n"
        "    v  : float，当前速度，单位：米/秒(m/s) / current speed, m/s\n"
        "    r  : float，当前到中心天体距离，单位：米(m)，必须 > 0 / distance, meters\n"
        "    mu : float，中心天体 G*M，单位：m³/s²，必须 > 0 / gravitational parameter\n\n"
        "返回 / Returns:\n"
        "    float，离心率 e。e<1 椭圆，e=1 抛物线，e>1 双曲线。\n"
        "    eccentricity e: e<1 ellipse, e=1 parabola, e>1 hyperbola."
    },
    {
        "apsis_distances",
        apsis_distances,
        METH_VARARGS,
        "apsis_distances(a, e) -> tuple\n\n"
        "计算椭圆轨道的近点距离和远点距离。\n"
        "Calculate periapsis and apoapsis distances of an elliptic orbit.\n\n"
        "参数 / Parameters:\n"
        "    a : float，半长轴，单位：米(m)，必须 > 0 / semi-major axis, meters\n"
        "    e : float，离心率，必须满足 0 <= e < 1 / eccentricity, 0 <= e < 1\n\n"
        "返回 / Returns:\n"
        "    tuple(r_periapsis, r_apoapsis)，单位：米(m) / meters"
    },
    {
        "lagrange_point",
        lagrange_point,
        METH_VARARGS,
        "lagrange_point(m1, m2, r, point) -> tuple\n\n"
        "计算两个天体系统的拉格朗日点位置（以两天体连线中点为原点）。\n"
        "Calculate Lagrange point position (origin at midpoint of the two bodies).\n\n"
        "参数 / Parameters:\n"
        "    m1    : float，主天体质量，单位：千克(kg)，必须 > 0 / primary mass, kg\n"
        "    m2    : float，次天体质量，单位：千克(kg)，必须 > 0 / secondary mass, kg\n"
        "    r     : float，两天体间距，单位：米(m)，必须 > 0 / separation, meters\n"
        "    point : int，拉格朗日点编号，1~5 / Lagrange point index, 1~5\n\n"
        "返回 / Returns:\n"
        "    tuple(x, y)，单位：米(m) / meters"
    },
    {
        "multi_body_gravity",
        multi_body_gravity,
        METH_VARARGS,
        "multi_body_gravity(bodies) -> tuple\n\n"
        "计算多体系统中某点受到的总引力合力。\n"
        "Calculate the total gravitational force vector in a multi-body system.\n\n"
        "参数 / Parameters:\n"
        "    bodies : list，元素为 (mass, x, y) 三元组\n"
        "             list of (mass, x, y) tuples\n"
        "             mass 单位：千克(kg) / mass in kg\n"
        "             x, y 单位：米(m) / positions in meters\n\n"
        "返回 / Returns:\n"
        "    tuple(fx, fy)，合力分量，单位：牛顿(N) / force components in Newtons"
    },
    {
        "planet_information",
        planet_information,
        METH_VARARGS,
        "planet_information(name) -> str\n\n"
        "查询内置星球参数（中英文双语）。\n"
        "Query built-in planet data (bilingual CN/EN).\n\n"
        "参数 / Parameters:\n"
        "    name : str，星球名称，支持中文或英文 / planet name, Chinese or English\n"
        "           支持 / supported: 水星/金星/地球/火星/木星/土星/天王星/海王星/月球/太阳\n"
        "           Mercury/Venus/Earth/Mars/Jupiter/Saturn/Uranus/Neptune/Moon/Sun\n\n"
        "返回 / Returns:\n"
        "    str，包含半径(m)、质量(kg)、表面重力(m/s²)、逃逸速度(m/s)、\n"
        "         公转周期(s)、自转周期(s)、大气层高度(m)、海洋深度(m)\n"
        "    radius(m), mass(kg), surface gravity(m/s^2), escape velocity(m/s),\n"
        "    orbital period(s), rotation period(s), atmosphere height(m), ocean depth(m)"
    },
    {
        "ptmation",
        planet_information,
        METH_VARARGS,
        "ptmation(name) -> str\n\n"
        "planet_information 的别名，参数和返回值完全相同。\n"
        "Alias of planet_information; same parameters and return value."
    },
    {NULL, NULL, 0, NULL}
};

// ========== 模块定义 / Module definition ==========

static struct PyModuleDef kepler_module = {
    PyModuleDef_HEAD_INIT,
    "oiiaio_web_kepler",
    "轨道力学计算器 / Orbital Mechanics Calculator\n\n"
    "提供引力、轨道速度、霍曼转移、开普勒定律、拉格朗日点、多体引力、"
    "星球参数查询等功能。\n"
    "Provides gravity, orbital velocity, Hohmann transfer, Kepler's laws, "
    "Lagrange points, multi-body gravity, and planet data queries.\n\n"
    "所有物理量使用国际单位制（SI）：\n"
    "All quantities use SI units:\n"
    "    长度 length : 米 m\n"
    "    质量 mass   : 千克 kg\n"
    "    时间 time   : 秒 s\n"
    "    速度 speed  : 米/秒 m/s\n"
    "    力 force    : 牛顿 N",
    -1,
    KeplerMethods
};

// ========== 初始化 / Module init ==========

PyMODINIT_FUNC PyInit_oiiaio_web_kepler(void) {
    return PyModule_Create(&kepler_module);
}
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>
#include <cstdio>
#include <thread>
using namespace std;

namespace AlgorithmParams {
    int D = 1;
    double G = 1;
    double EPS = 0.0000000000001;
    double T = 0.1;
    double theta = 1;
};

struct Point {
    double M{};
    double V_X{}, V_Y{}, V_Z{};
    double F_X{}, F_Y{}, F_Z{};
    double X{}, Y{}, Z{};
    Point(double x, double y, double z, double m, double v_x, double v_y, double v_z) {
        X = x;
        Y = y;
        Z = z;
        M = m;
        V_X = v_x;
        V_Y = v_y;
        V_Z = v_z;
    }

    double dist(Point p) const {
        return sqrt((X - p.X) * (X - p.X) + (Y - p.Y) * (Y - p.Y) + (Z - p.Z) * (Z - p.Z));
    }
    double dist(double px, double py, double pz) const {
        return sqrt((X - px) * (X - px) + (Y - py) * (Y - py) + (Z - pz) * (Z - pz));
    }
    void AddPowerFromPoint(Point *p) const {
        double d = dist(*p);
        if (abs(d) < AlgorithmParams::EPS) {
            return;
        }
        p->F_X -= AlgorithmParams::G * M * p->M * (p->X - X) / pow(d,3);
        p->F_Y -= AlgorithmParams::G * M * p->M * (p->Y - Y) / pow(d,3);
        p->F_Z -= AlgorithmParams::G * M * p->M * (p->Z - Z) / pow(d,3);
    }
    void UpdatedPointState() {
        double A_X, A_Y, A_Z;
        A_X = F_X / M;
        A_Y = F_Y / M;
        A_Z = F_Z / M;
        X = X + V_X * AlgorithmParams::T + A_X * AlgorithmParams::T * AlgorithmParams::T / 2;
        Y = Y + V_Y * AlgorithmParams::T + A_Y * AlgorithmParams::T * AlgorithmParams::T / 2;
        Z = Z + V_Z * AlgorithmParams::T + A_Z * AlgorithmParams::T * AlgorithmParams::T / 2;
        V_X = V_X + AlgorithmParams::T * A_X;
        V_Y = V_Y + AlgorithmParams::T * A_Y;
        V_Z = V_Z + AlgorithmParams::T * A_Z;
    }

};

struct Field {
    Field* child_fields[8] = {};
    vector<Point> points;

    // field params
    double X_AVG = 0, Y_AVG = 0, Z_AVG = 0;
    double M_SUM = 0, M_X_SUM = 0, M_Y_SUM = 0, M_Z_SUM = 0;
    double size = -1;

    void AddPowerFromField(Point* p) {
        double P_X = AlgorithmParams::G * p->M * (Multipole_X_1(p));
        double P_Y = AlgorithmParams::G * p->M * (Multipole_Y_1(p));
        double P_Z = AlgorithmParams::G * p->M * (Multipole_Z_1(p));
        p->F_X -= P_X;// + Multipole_X_2(p) + Multipole_X_3(p) + Multipole_X_4(p));
        p->F_Y -= P_Y;// + Multipole_Y_2(p) + Multipole_Y_3(p) + Multipole_Y_4(p));
        p->F_Z -= P_Z;// + Multipole_Z_2(p) + Multipole_Z_3(p) + Multipole_Z_4(p));

        Point pp = *p;
        pp.F_X = 0; pp.F_Y = 0; pp.F_Z = 0;
        for (auto ppp: points) {
            ppp.AddPowerFromPoint(&pp);
        }
        if (abs(pp.F_X + P_X) > 5) {
            cerr << pp.F_X << " " << P_X << "???\n";
        }
    }

    // X
    double Multipole_X_1(Point* p) {
        return M_SUM * (p->X - X_AVG) / pow(p->dist(X_AVG, Y_AVG, Z_AVG), 3);
    }
    double Multipole_X_2(Point* p) {
        return (M_X_SUM - M_SUM * X_AVG) * -(-2 * p->X * p->X + 4 * p->X * X_AVG + p->Y * p->Y - 2 * p->Y * Y_AVG + p->Z * p->Z - 2 * p->Z * Z_AVG - 2 * X_AVG * X_AVG + Y_AVG * Y_AVG + Z_AVG * Z_AVG) / pow(p->dist(X_AVG, Y_AVG, Z_AVG), 5);
    }
    double Multipole_X_3(Point* p) {
        return (M_Y_SUM - M_SUM * Y_AVG) * (3 * (p->X - X_AVG) * (p->Y - Y_AVG) / pow(p->dist(X_AVG, Y_AVG, Z_AVG), 5));
    }
    double Multipole_X_4(Point* p) {
        return (M_Z_SUM - M_SUM * Z_AVG) * (3 * (p->X - X_AVG) * (p->Z - Z_AVG) / pow(p->dist(X_AVG, Y_AVG, Z_AVG), 5));
    }

    // Y
    double Multipole_Y_1(Point* p) {
        return M_SUM * (p->Y - Y_AVG) / pow(p->dist(X_AVG, Y_AVG, Z_AVG), 3);
    }
    double Multipole_Y_2(Point* p) {
        return (M_X_SUM - M_SUM * X_AVG) * (3 * (p->X - X_AVG) * (p->Y - Y_AVG)) / pow(p->dist(X_AVG, Y_AVG, Z_AVG), 5);
    }
    double Multipole_Y_3(Point* p) {
        return (M_Y_SUM - M_SUM * Y_AVG) * -((p->X * p->X - 2 * p->X * X_AVG - 2 * p->Y * p->Y + 4 * p->Y * Y_AVG + p->Z * p->Z - 2 * p->Z * Z_AVG + X_AVG * X_AVG - 2 * Y_AVG * Y_AVG + Z_AVG * Z_AVG))/ pow(p->dist(X_AVG, Y_AVG, Z_AVG), 5);
    }
    double Multipole_Y_4(Point* p) {
        return (M_Z_SUM - M_SUM * Z_AVG) * (3 * (p->Y - Y_AVG) * (p->Z - Z_AVG) / pow(p->dist(X_AVG, Y_AVG, Z_AVG), 5));
    }

    // Z
    double Multipole_Z_1(Point* p) {
        return M_SUM * (p->Z - Z_AVG) / pow(p->dist(X_AVG, Y_AVG, Z_AVG), 3);
    }
    double Multipole_Z_2(Point* p) {
        return (M_X_SUM - M_SUM * X_AVG) * (3 * (p->X - X_AVG) * (p->Z - Z_AVG)) / pow(p->dist(X_AVG, Y_AVG, Z_AVG), 5);
    }
    double Multipole_Z_3(Point* p) {
        return (M_Y_SUM - M_SUM * Y_AVG) * (3 * (p->Y - Y_AVG) * (p->Z - Z_AVG) / pow(p->dist(X_AVG, Y_AVG, Z_AVG), 5));
    }
    double Multipole_Z_4(Point* p) {
        return (M_Z_SUM - M_SUM * Z_AVG) * -((p->X * p->X - 2 * p->X * X_AVG + p->Y * p->Y - 2 * p->Y * Y_AVG - 2 * p->Z * p->Z + 4 * p->Z * Z_AVG + X_AVG * X_AVG + Y_AVG * Y_AVG - 2 * Z_AVG * Z_AVG)) / pow(p->dist(X_AVG, Y_AVG, Z_AVG), 5);
    }

    double get_size() {
        if (size != -1) {
            return size;
        }
        double min_x = 1e18, max_x = -1e18, min_y = 1e18, max_y = -1e18;
        for (auto p: points) {
            min_x = min(min_x, p.X);
            max_x = max(max_x, p.X);
            min_y = min(min_y, p.Y);
            max_y = max(max_y, p.Y);
        }
        size = min((max_x - min_x), (max_y - min_y));
        return size;
    }
};

vector<Point> GetOctet(int n, double x_0, double y_0, double z_0, const vector<Point> &points) {
    vector<Point> octet;
    for (const auto &p : points) {
        bool ox{}, oy{}, oz{};
        ox = p.X > x_0;
        oy = p.Y > y_0;
        oz = p.Z > z_0;
        if (((ox * 4) + (oy * 2) + oz) == n) {
            octet.push_back(p);
        }
    }
    return octet;
}

bool CheckInOctet(int n, double x_0, double y_0, double z_0, Point* p) {
    bool ox{}, oy{}, oz{};
    ox = p->X > x_0;
    oy = p->Y > y_0;
    oz = p->Z > z_0;
    if (((ox * 4) + (oy * 2) + oz) == n) {
        return true;
    }
    return false;
}

void field_initializer(Field* f) {
    // partition field
    for (auto &p : f->points) {
        f->M_SUM += p.M;
        f->M_X_SUM += p.M * p.X;
        f->M_Y_SUM += p.M * p.Y;
        f->M_Z_SUM += p.M * p.Z;
    }
    f->X_AVG = f->M_X_SUM / f->M_SUM; f->Y_AVG = f->M_Y_SUM / f->M_SUM; f->Z_AVG = f->M_Z_SUM / f->M_SUM;
    if (f->points.size() > AlgorithmParams::D) {
        for (int octet = 0; octet < 8; octet++) {
            vector<Point> points = GetOctet(octet, f->X_AVG, f->Y_AVG, f->Z_AVG, f->points);
            if (!points.empty()) {
                f->child_fields[octet] = new Field;
                f->child_fields[octet]->points = points;
            }
        }

        for (auto & child_field : f->child_fields) {
            if (child_field != nullptr) {
                field_initializer(child_field);
            }
        }
    }
}

void calc_power_from_field_to_point(Field* f, Point *p) {
    bool last_field = true;
    for (int octet = 0; octet < 8; octet++) {
        double DD = 0, r = 1e18;
        if (f->child_fields[octet] != nullptr) {
            DD = f->child_fields[octet]->get_size();
            r = p->dist(f->child_fields[octet]->X_AVG, f->child_fields[octet]->Y_AVG, f->child_fields[octet]->Z_AVG);
        }
        // cerr << DD / r << "\n";
        if (!CheckInOctet(octet, f->X_AVG, f->Y_AVG, f->Z_AVG, p) && f->child_fields[octet] != nullptr && DD / r < 0.05) {

            f->child_fields[octet]->AddPowerFromField(p);
        } else if (f->child_fields[octet] != nullptr) {
            last_field = false;
            calc_power_from_field_to_point(f->child_fields[octet], p);
        }
    }
    if (last_field) {
        for (auto& point: f->points) {
            point.AddPowerFromPoint(p);
        }
    }
}

vector<Point> direct_calc_power(vector<Point> points) {
    vector<Point> updated_points;
    for (auto cur_p: points) {
        for (auto p: points) {
            p.AddPowerFromPoint(&cur_p);
        }
        cur_p.UpdatedPointState();
        updated_points.emplace_back(cur_p);
    }
    return updated_points;
}

void cerr_errors(vector<Point> correct, vector<Point> my_answer) {
    double F_X_ERR{}, F_Y_ERR{}, F_Z_ERR{};
    double F_ERR;
    for (int i = 0; i < correct.size(); i++) {
        double myF = sqrt(pow(my_answer[i].F_X, 2) + pow(my_answer[i].F_Y, 2) + pow(my_answer[i].F_Z, 2));
        double correctF = sqrt(pow(correct[i].F_X, 2) + pow(correct[i].F_Y, 2) + pow(correct[i].F_Z, 2));
        if (i < 100) {
            cerr << myF << " " << correctF << "\n";
        }
        F_ERR += abs(myF - correctF) / abs(correctF);
        F_X_ERR += abs(correct[i].F_X - my_answer[i].F_X) / abs(correct[i].F_X);
        F_Y_ERR += abs(correct[i].F_Y - my_answer[i].F_Y) / abs(correct[i].F_Y);
        F_Z_ERR += abs(correct[i].F_Z - my_answer[i].F_Z) / abs(correct[i].F_Z);
    }
    F_ERR /= correct.size();
    F_X_ERR /= correct.size();
    F_Y_ERR /= correct.size();
    F_Z_ERR /= correct.size();

    cerr << "F error = " << F_ERR * 100 << "%\n";
    cerr << "F on OX error = " << F_X_ERR * 100 << "%\n";
    cerr << "F on OY error = " << F_Y_ERR * 100 << "%\n";
    cerr << "F on OZ error = " << F_Z_ERR * 100 << "%\n";
}

int main() {
    freopen("C:\\Diplom\\input.txt", "r", stdin);

    int N;
    cin >> N;

    Field *f = new Field;
    vector<Point> points;
    for (int i = 0; i < N; ++i) {
        double m, v_x, v_y, v_z, x, y, z;
        cin >> x >> y >> z >> m >> v_x >> v_y >> v_z;
        f->points.emplace_back(x, y, z, m, v_x, v_y, v_z);
    }
    points = f->points;
    
    
    freopen("C:\\Diplom\\output.txt", "w", stdout);
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    
    field_initializer(f);
    
    for (int i = 0; i < N; i++) {
        calc_power_from_field_to_point(f, &points[i]);
    }
    
    cout << N << '\n';
    for (int i = 0; i < N; i++) {
        points[i].UpdatedPointState();
        cout << 
        points[i].X << " " << 
        points[i].Y << " " << 
        points[i].Z << " " << 
        points[i].M << " " << 
        points[i].V_X << " " << 
        points[i].V_Y << " " << 
        points[i].V_Z << "\n";
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cerr << "Execution Time = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1e6 << "[s]" << std::endl;

    vector<Point> correct_points = direct_calc_power(f->points);

    cerr_errors(correct_points, points);
    return 0;
}

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>
#include <cstdio>
#include <thread>
#include <complex>
#include <map>
#include <thread>
#include <mutex>

std::mutex g_lock;
using namespace std;

// ALGORITHM PARAMS:
const long double EPS = 1e-15;
const int P = 10;
const int D = 1 << 15;
const int MIN_Y = 0, MAX_Y = (1 << 29);
const int MIN_X = 0, MAX_X = (1 << 29);
const int thread_count = 16;
//


// FUNCTIONS:
bool USED_I[400];
complex<long double> RES_I[400];
const complex<long double> I(0, 1);
complex<long double> pow_i(int n) {
    if (USED_I[n + 3 * P]) {
        return RES_I[n + 3 * P];
    }
    USED_I[n + 3 * P] = 1;
    RES_I[n + 3 * P] = pow(I, n);
    return RES_I[n + 3 * P];
}

long double fact[200];
long double factorial(int n) {
    if (n == 0) {
        return 1;
    }
    if (fact[n] != 0) {
        return fact[n];
    }
    fact[n] = n * factorial(n - 1);
    return fact[n];
}

bool USED_A_m_n[200][200];
long double RES_A_m_n[200][200];
long double A_m_n(int m, int n) {
    if (n - m < 0 || n + m < 0) {
        return 0;
    }
    if (USED_A_m_n[m + 2 * P][n + 2 * P]) {
        return RES_A_m_n[m + 2 * P][n + 2 * P];
    }
    USED_A_m_n[m + 2 * P][n + 2 * P] = 1;
    RES_A_m_n[m + 2 * P][n + 2 * P] = pow(-1, n & 1) / sqrt(factorial(n - m) * factorial(n + m));
    return RES_A_m_n[m + 2 * P][n + 2 * P];
}


bool USED_P_m_n[200][200];
long double RES_P_m_n[200][200];
std::mutex p_lock;
long double P_m_n(int m, int n) {
    long double res = 0;
    if (m < 0) {
        exit(-1);
    }
    if (((n - m) & 1) != 0) {
        return 0;
    }

    if (USED_P_m_n[m + 2 * P][n + 2 * P]) {
        return RES_P_m_n[m + 2 * P][n + 2 * P];
    }
    USED_P_m_n[m + 2 * P][n + 2 * P] = 1;
    int k = (n - m) / 2;
    res = pow(-1, k & 1) * factorial(2 * n - 2 * k) / factorial(k) / factorial(n - k) / factorial(n - m - 2 * k);
    res = pow(-1, m & 1) * res / (1ll << n);
    RES_P_m_n[m + 2 * P][n + 2 * P] = res;
    return res;
}

struct Point {
    long double x{}, y{};
    long double q{};
    long double r{}, phi{};
    long double pot{}, F_x{}, F_y{};
    Point(long double xx, long double yy, long double q_p) {
        q = q_p;
        x = xx;
        y = yy;
        r = sqrt(x * x + y * y); // don't use 0, 0, 0
        if (r == 0) {
            phi = 0;
            return;
        }
        if (x == 0 && y == 0) {
            phi = 0;
            return;
        }
        if (y < 0) {
            phi = -1 * acos(x / sqrt(x * x + y * y));
        } else {
            phi = acos(x / sqrt(x * x + y * y));
        }
    }
    long double dist(Point *p) const {
        return sqrt((x - p->x) * (x - p->x) + (y - p->y) * (y - p->y));
    }
    long double dist(const Point &p) const {
        return sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y));
    }
    int get_octet(const Point &p) {
        int b1 = p.x < x;
        int b2 = p.y < y;
        return 2 * b1 + b2;
    }

    Point() {}
};

bool USED_Y_m_n[200][200];
long double RES_Y_m_n[200][200];
const complex<long double> Y_m_n(int m, int n, const Point &p) {
    if (abs(m) > n) {
        return 0;
    }
    if (USED_Y_m_n[m + 2 * P][n + 2 * P]) {
        return RES_Y_m_n[m + 2 * P][n + 2 * P] * complex<long double>(cos(m * p.phi), sin(m * p.phi));
    }
    USED_Y_m_n[m + 2 * P][n + 2 * P] = 1;
    RES_Y_m_n[m + 2 * P][n + 2 * P] = sqrt(factorial(n - abs(m)) / factorial(n + abs(m))) * P_m_n(abs(m), n);
    return RES_Y_m_n[m + 2 * P][n + 2 * P] * complex<long double>(cos(m * p.phi), sin(m * p.phi));
}

complex<long double> Y_m_n_phi(int m, int n, Point p) {
    complex<long double> z(m * -sin(m * p.phi), m * cos(m * p.phi));
    return RES_Y_m_n[m + 2 * P][n + 2 * P] * z;
}

//

vector<Point> input() {
    freopen("/Users/kuzalexey/fmm2D_SH/input.txt", "r", stdin);
    int N;
    cin >> N;
    vector<Point> points;
    for (int i = 0; i < N; ++i) {
        long double x, y, m;
        cin >> x >> y >>  m;
        points.push_back(Point(x, y, m));
    }
    return points;
}

struct Field {
    vector<Point> points;
    vector<vector<complex<long double> > > M, L;
    vector<Field*> neighbours, interaction_list;
    Point x_c;
    Field* child_fields[4] = {nullptr, nullptr, nullptr, nullptr};
    Field* parent;
    int level{};
    int l_x{}, r_x{}, l_y{}, r_y{};

    bool is_last() {
        return points.empty() || r_x - l_x == D;
    }

    complex<long double> get_M_m_n(int m, int n) {
        return M[m + P][n];
    }

    void init_M() {
        M.resize(2 * P);
        for (int i = 0; i < 2 * P; i++) {
            M[i].resize(P);
        }
        for (int k = 0; k < points.size(); k++) {
            Point shift_p(points[k].x - x_c.x, points[k].y - x_c.y, points[k].q);
            long double r_n = 1;
            for (int n = 0; n < P; n++) {
                for (int m = -n; m <= n; m++) {
                    M[m + P][n] += shift_p.q * r_n * Y_m_n(-m, n, shift_p);
                }
                r_n *= shift_p.r;
            }
        }
    }

    complex<long double> get_Q(Point *x) {
        complex<long double> res(0, 0);
        Point shift_p(x->x - x_c.x, x->y - x_c.y, 1);
        long double r_n = 1;
        for (int n = 0; n < P; n++) {
            r_n *= shift_p.r;
            for (int m = -n; m <= n; m++) {
                res += get_M_m_n(m, n) / r_n * Y_m_n(m, n, shift_p);
            }
        }
        return res;
    }

    complex<long double> get_F_x(Point x) {
        complex<long double> res_r(0, 0), res_phi(0, 0);
        Point shift_p(x.x - x_c.x, x.y - x_c.y, 1);
        x = shift_p;
        for (int n = 0; n < P; n++) {
            for (int m = -n; m <= n; m++) {
                res_r += get_M_m_n(m, n) * (long double)(-n - 1) / pow(x.r, n + 2) * Y_m_n(m, n, x);
                res_phi += get_M_m_n(m, n) / pow(x.r, n + 1) * Y_m_n_phi(m, n, x);
            }
        }
        return -(res_r * x.x / x.r + res_phi * (-x.y/(x.x * x.x + x.y * x.y)));
    }

    void init_M2L(Field *f) {
        if (f->L.size() == 0) {
            f->L.resize(2 * P);
            for (int i = 0; i < 2 * P; i++) {
                f->L[i].resize(P);
            }
        }
        const Point x_dif = Point(x_c.x - f->x_c.x, x_c.y - f->x_c.y, 0);
        long double r_j_1 = 1, r_n;
        for (int j = 0; j < P; j++) {
            r_j_1 *= x_dif.r;
            for (int k = -j; k <= j; k++) {
                r_n = 1;
                for (int n = 0; n < P; n++) {
                    for (int m = -n; m <= n; m++) {
                        f->L[k + P][j] += (M[m + P][n] * pow_i(abs(k - m) - abs(k) - abs(m)) * A_m_n(m, n) * A_m_n(k, j) *
                                Y_m_n(m - k, j + n, x_dif)) /
                                        (pow(-1, n & 1) * A_m_n(m - k, j + n) * r_j_1 * r_n);
                    }
                    r_n *= x_dif.r;
                }
            }
        }
    }

    void init_L2L(Field *f) {
        if (f->L.size() == 0) {
            f->L.resize(2 * P);
            for (int i = 0; i < 2 * P; i++) {
                f->L[i].resize(P);
            }
        }
        if (L.size() == 0) {
            return;
        }
        const Point x_dif = Point(x_c.x - f->x_c.x, x_c.y - f->x_c.y, 0);
        long double rr = 1;
        for (int k = -P; k <= P; k++) {
            for (int j = 0; j < P; j++) {
                if (abs(k) > j) {
                    continue;
                }
                rr = 1;
                for (int n = j; n < P; n++) {
                    for (int m = -n; m <= n; m++) {
                        f->L[k + P][j] += L[m + P][n] * pow_i(abs(m) - abs(m - k) - abs(k)) * A_m_n(m - k, n - j) * A_m_n(k, j) * Y_m_n(m - k, n - j, x_dif) * rr / (pow(-1, (n + j) & 1) * A_m_n(m, n));
                    }
                    rr *= x_dif.r;
                }
            }
        }
    }

    complex<long double> get_q_local(Point* x) {
        complex<long double> res(0, 0);
        Point shift_p(x->x - x_c.x, x->y - x_c.y, 1);
        long double r_n = 1;
        for (int n = 0; n < P; n++) {
            for (int m = -n; m <= n; m++) {
                res += L[m + P][n] * r_n * Y_m_n(m, n, shift_p);
            }
            r_n *= shift_p.r;
        }
        return res;
    }

    pair<long double, long double> get_f_local(Point* x) {
        long double f_x, f_y; f_x = 0; f_y = 0;
        long double f_r, f_phi; f_r = 0; f_phi = 0;
        Point shift_p(x->x - x_c.x, x->y - x_c.y, 1);
        long double r_n = 1;
        long double r_n_1 = 1;
        for (int n = 0; n < P; n++) {
            for (int m = -n; m <= n; m++) {
                if (n >= 1) {
                    f_r += n * (L[m + P][n] * r_n_1 * Y_m_n(m, n, shift_p)).real();
                }
                f_phi += (L[m + P][n] * r_n * Y_m_n_phi(m, n, shift_p)).real();
            }
            if (n != 0) {
                r_n_1 *= shift_p.r;
            }
            r_n *= shift_p.r;
        }
        f_x = -(f_r * shift_p.x / shift_p.r + f_phi * (-shift_p.y/(shift_p.x * shift_p.x + shift_p.y * shift_p.y)));
        f_y = -(f_r * shift_p.y / shift_p.r + f_phi * (shift_p.x/(shift_p.x * shift_p.x + shift_p.y * shift_p.y)));
        return {f_x, f_y};
    }

    // DEBUG

    complex<long double> get_direct_Q(Point* x) {
        complex<long double> res(0, 0);
        for (const auto &p: points) {
            res += p.q / p.dist(x);
        }
        return res;
    }
    pair<long double, long double> get_direct_F(Point *x) {
        long double res_x = 0, res_y = 0, dist = 0;
        for (int i = 0; i < points.size(); i++) {
            dist = points[i].dist(x);
            if (dist < EPS) {
                continue;
            }
            res_x += points[i].q * (x->x - points[i].x) / (dist * dist * dist);
            res_y += points[i].q * (x->y - points[i].y) / (dist * dist * dist);
        }
        return {res_x, res_y};
    }
};

vector<map<pair<int, int>, Field* > > MAP(40);
vector<vector<Field*>>fields(40);

void init_KDTREE(Field* f) {
    // init f
    if (f->parent == nullptr) {
        f->level = 0;
    } else {
        f->level = f->parent->level + 1;
    }
    f->x_c = Point((f->r_x + f->l_x) / 2, (f->r_y + f->l_y) / 2, 0);
    MAP[f->level][{f->l_x, f->l_y}] = f;
    fields[f->level].push_back(f);

    if (f->is_last()) {
        return;
    }

    // init childs
    for (auto point: f->points) {
        int octet = point.get_octet(f->x_c);
        if (f->child_fields[octet] == nullptr) {
            f->child_fields[octet] = new Field;
            f->child_fields[octet]->parent = f;
            f->child_fields[octet]->l_x = f->l_x;
            f->child_fields[octet]->r_x = f->r_x;
            f->child_fields[octet]->l_y = f->l_y;
            f->child_fields[octet]->r_y = f->r_y;
            if (octet & 1) {
                f->child_fields[octet]->l_y = f->x_c.y;
            } else {
                f->child_fields[octet]->r_y = f->x_c.y;
            }
            if (octet & 2) {
                f->child_fields[octet]->l_x = f->x_c.x;
            } else {
                f->child_fields[octet]->r_x = f->x_c.x;
            }
        }
        f->child_fields[octet]->points.push_back(point);
    }
    for (int i = 0; i < 4; i++) {
        if (f->child_fields[i] != nullptr) {
            init_KDTREE(f->child_fields[i]);
        }
    }
}

void init_near() {
    for (int level = 0; level < fields.size(); level++) {
        for (auto f: fields[level]) {
            Field* parent = f->parent;
            if (parent == nullptr) {
                continue;
            }
            int D_parent = parent->r_x - parent->l_x;
            int D_f = f->r_x - f->l_x;
            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    int nl_x = parent->l_x + i * D_parent;
                    int nl_y = parent->l_y + j * D_parent;
                    Field* neigh_parent = MAP[parent->level][{nl_x, nl_y}];
                    if (neigh_parent == nullptr) {
                        continue;
                    }
                    for (auto neigh: neigh_parent->child_fields) {
                        if (neigh == nullptr || neigh == f) {
                            continue;
                        }
                        int dx = abs(neigh->l_x - f->l_x);
                        int dy = abs(neigh->l_y - f->l_y);
                        if (dx >= 2 * D_f || dy >= 2 * D_f) {
                            f->interaction_list.push_back(neigh);
                        } else {
                            f->neighbours.push_back(neigh);
                        }
                    }
                }
            }
        }
    }
}

void init_M() {
    for (int level = 0; level < fields.size(); level++) {
        thread threads[thread_count];
        int D = (fields[level].size() + thread_count - 1) / thread_count;
        for (int cnt = 0; cnt < thread_count; cnt++) {
            threads[cnt] = thread([D,cnt,level](){
                for (int i = cnt * D; i < min(int(fields[level].size()), cnt * D + D);i++) {
                    fields[level][i]->init_M();
                }
            });
        }
        for (int cnt = 0; cnt < thread_count; cnt++) {
            threads[cnt].join();
        }
    }
}

void init_L() {
    for (int level = 0; level < fields.size(); level++) {
        thread threads[thread_count];
        int D = (fields[level].size() + thread_count - 1) / thread_count;
        for (int cnt = 0; cnt < thread_count; cnt++) {
            threads[cnt] = thread([D,cnt,level](){
                for (int i = cnt * D; i < min(int(fields[level].size()), cnt * D + D);i++) {
                    // interaction_list A -> B
                    for (auto near: fields[level][i]->interaction_list) {
                        near->init_M2L(fields[level][i]);
                    }

                    // B -> B to childs
                    for (auto child: fields[level][i]->child_fields) {
                        if (child == nullptr) {
                            continue;
                        }
                        fields[level][i]->init_L2L(child);
                    }
                }
            });
        }
        for (int cnt = 0; cnt < thread_count; cnt++) {
            threads[cnt].join();
        }
    }
}

void cacl_pf(Field *f, Point *p) {
    if (f->is_last()) {
        auto f_local = f->get_f_local(p);
        p->F_x += get<0>(f_local);
        p->F_y += get<1>(f_local);
        p->pot += f->get_q_local(p).real();
        pair<long double, long double> F;
        for (const auto& neigh: f->neighbours) {
            F = neigh->get_direct_F(p);
            p->F_x += F.first;
            p->F_y += F.second;
            p->pot += neigh->get_direct_Q(p).real();
        }
        long double r_p = 0;
        for (auto point: f->points) {
            if (point.dist(p) < EPS) {
                continue;
            }
            r_p = point.dist(p);
            p->pot += point.q / r_p;
            p->F_x += point.q * (p->x - point.x) / (r_p * r_p * r_p);
            p->F_y += point.q * (p->y - point.y) / (r_p * r_p * r_p);
        }
    } else {
        cacl_pf(f->child_fields[p->get_octet(f->x_c)], p);
    }
}

vector<Point> direct_calc(vector<Point> points) {
    vector<Point> updated_points;
    thread threads[thread_count];
    int D = (points.size() + thread_count - 1) / thread_count;
    D = 1e5;
    for (int cnt = 0; cnt < thread_count; cnt++) {
        threads[cnt] = thread([cnt,&points,D,&updated_points](){
            for (int i = cnt * D; i < min(int(points.size()), cnt * D + D);i++) {
                Point cur_p = points[i];
                for (auto p: points) {
                    if (cur_p.dist(&p) < EPS) {
                        continue;
                    }
                    cur_p.pot += p.q / cur_p.dist(&p);
                    cur_p.F_x += p.q * (cur_p.x - p.x) / pow(cur_p.dist(&p), 3);
                    cur_p.F_y += p.q * (cur_p.y - p.y) / pow(cur_p.dist(&p), 3);
                }
                g_lock.lock();
                updated_points.emplace_back(cur_p);
                g_lock.unlock();
            }
        });
    }
    for (int cnt = 0; cnt < thread_count; cnt++) {
        threads[cnt].join();
    }
    return updated_points;
}

void cerr_errors(vector<Point> correct, vector<Point> my_answer) {
    long double POT_ERR = 0, FX_ERR = 0, FY_ERR = 0;
    for (int i = 0; i < correct.size(); i++) {
        if (i < 50) {
            //            cerr << correct[i].F_x << " " << my_answer[i].F_x << "\n";
            //            cerr << correct[i].F_y << " " << my_answer[i].F_y << "\n";
        }
        POT_ERR += abs(correct[i].pot - my_answer[i].pot) / abs(correct[i].pot);
        FX_ERR += (abs(correct[i].F_x - my_answer[i].F_x)) / (abs(correct[i].F_x));
        FY_ERR += (abs(correct[i].F_y - my_answer[i].F_y)) / (abs(correct[i].F_y));
    }
    POT_ERR /= correct.size();
    FX_ERR /= correct.size();
    FY_ERR /= correct.size();
    cerr << "POT relative error = " << POT_ERR * 100 << "%\n";
    cerr << "F_x relative error = " << FX_ERR * 100 << "%\n";
    cerr << "F_y relative error = " << FY_ERR * 100 << "%\n";
}

int main() {
    // init functions
    for (int i = 0; i < 200; i++) {
        factorial(i);
    }
    for (int i = -3 * P; i < 100; i++) {
        pow_i(i);
    }
    for (int i = -2 * P; i < 100; i++) {
        for (int j = -2 * P; j < 100; j++) {
            A_m_n(i, j);
        }
    }
    for (int i = 0; i < 4 * P; i++) {
        for (int j = 0; j < 4 * P; j++) {
            P_m_n(i, j);
        }
    }
    for (int i = -2 * P; i < 4 * P; i++) {
        for (int j = -2 * P; j < 4 * P; j++) {
            Y_m_n(i, j, Point(10, 10, 10));
        }
    }




    cout.precision(20);

    Field *f = new Field;
    vector<Point> points = input();
    int N = points.size();
    f->points = points;
    f->l_y = MIN_Y; f->r_y = MAX_Y; f->l_x = MIN_X; f->r_x = MAX_X;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now(), beginALL =  std::chrono::steady_clock::now();
    init_KDTREE(f);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cerr << "Init KD Tree Time = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1e6 << "[s]" << std::endl;

    begin = std::chrono::steady_clock::now();
    init_near();
    end = std::chrono::steady_clock::now();
    std::cerr << "Init Near Time = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1e6 << "[s]" << std::endl;

    begin = std::chrono::steady_clock::now();
    init_M();
    end = std::chrono::steady_clock::now();
    std::cerr << "Init M-Multipoles Time = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1e6 << "[s]" << std::endl;

    begin = std::chrono::steady_clock::now();
    init_L();
    end = std::chrono::steady_clock::now();
    std::cerr << "Init L-Multipoles Time = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1e6 << "[s]" << std::endl;

    begin = std::chrono::steady_clock::now();
    thread threads[thread_count];
    int D = (N + thread_count - 1) / thread_count;
    for (int cnt = 0; cnt < thread_count; cnt++) {
        threads[cnt] = thread([cnt, &D, &N, &points, &f](){
            for (int i = cnt * D; i < min(N, cnt * D + D); i++) {
                cacl_pf(f, &points[i]);
            }
        });
    }
    for (int i = 0; i < thread_count; i++){
        threads[i].join();
    }
    end = std::chrono::steady_clock::now();
    std::cerr << "Calc potentials Time = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1e6 << "[s]" << std::endl;
    std::cerr << "ALL FMM Time = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - beginALL).count()) / 1e6 << "[s]" << std::endl;

    begin = std::chrono::steady_clock::now();
    vector<Point> correct_points = direct_calc(f->points);
    end = std::chrono::steady_clock::now();
    std::cerr << "Direct calc potentials Time = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1e6 << "[s]" << std::endl;
    cerr_errors(correct_points, points);
    return 0;
}

#include "gnss.h"

namespace gnssins {

// gnss.cpp
// ---------------------------------------------------------------------------
// 本文件实现 GNSS 和 RTK 解算。
//
// 处理链条：
// 1. readRinexObs：读取流动站/基准站 RINEX 观测文件；
// 2. readRinexNav：读取广播星历；
// 3. satelliteState：由广播星历计算卫星位置和钟差；
// 4. sppSolve：用伪距单点定位得到初始坐标；
// 5. runRtkPass / solveDdEpoch：构造站间单差、站星双差观测方程，
//    通过 Kalman 滤波估计基线和模糊度浮点解；
// 6. lambdaIntegerSearch / 固定质量检查：尝试模糊度固定；
// 7. runRtkFloat：组合前后向 RTK 结果，输出 GnssSolution 序列。
//
// 输出的 GnssSolution 是 GNSS/INS 松组合的量测输入，主要包含：
// GPST、ECEF/LLH 坐标、NED 速度、质量标志 q、标准差和文字模式。
// ---------------------------------------------------------------------------
constexpr int kRtkPosStates = 3;
constexpr int kRtkVelStates = 3;
constexpr int kRtkClockStates = 2;
constexpr int kRtkAmbStart = kRtkPosStates + kRtkVelStates + kRtkClockStates;
static std::map<int, int> gGloFreqChannel;

// 将 RINEX 卫星系统字符转换为内部枚举。
static System systemFromChar(char c) {
    if (c == 'G') return System::GPS;
    if (c == 'C') return System::BDS;
    if (c == 'R') return System::GLO;
    return System::OTHER;
}

// 将内部卫星号转换为 RINEX 风格字符串，例如 G08、C19。
std::string satToString(const SatId &sat) {
    char c = sat.sys == System::GPS ? 'G' : (sat.sys == System::BDS ? 'C' : (sat.sys == System::GLO ? 'R' : '?'));
    std::ostringstream oss;
    oss << c << std::setw(2) << std::setfill('0') << sat.prn;
    return oss.str();
}

// 安全截取固定宽度字段，避免 RINEX 行长度不足时越界。
static std::string safeSubstr(const std::string &s, size_t pos, size_t len) {
    if (pos >= s.size()) return "";
    return s.substr(pos, std::min(len, s.size() - pos));
}

// 根据系统和观测类型判断频点下标；本工程主要使用双频 0/1。
static int freqIndex(System sys, const std::string &type) {
    if (type.size() < 2) return -1;
    const char kind = type[0];
    const char band = type[1];
    if (kind != 'C' && kind != 'L' && kind != 'D' && kind != 'S') return -1;
    if (sys == System::GPS) {
        if (band == '1') return 0;
        if (band == '2') return 1;
    } else if (sys == System::GLO) {
        if (band == '1') return 0;
        if (band == '2') return 1;
    } else if (sys == System::BDS) {
        if (band == '2' || band == '1') return 0;
        if (band == '6') return 1;
    }
    return -1;
}

// 同一频点可能有多个码型，使用优先级选择更常用、更稳定的观测类型。
static int obsPriority(System sys, const std::string &type) {
    if (type.size() < 3) return 0;
    if (sys == System::GPS) {
        if (type == "C1C" || type == "L1C" || type == "D1C" || type == "S1C") return 40;
        if (type == "C2W" || type == "L2W" || type == "D2W" || type == "S2W") return 40;
        if (type[1] == '1' || type[1] == '2') return 10;
    } else if (sys == System::GLO) {
        if (type == "C1C" || type == "L1C" || type == "D1C" || type == "S1C") return 40;
        if (type == "C2P" || type == "L2P" || type == "D2P" || type == "S2P") return 40;
        if (type == "C1P" || type == "L1P" || type == "D1P" || type == "S1P") return 30;
        if (type == "C2C" || type == "L2C" || type == "D2C" || type == "S2C") return 20;
        if (type[1] == '1' || type[1] == '2') return 10;
    } else if (sys == System::BDS) {
        if (type == "C2I" || type == "L2I" || type == "D2I" || type == "S2I") return 40;
        if (type == "C1I" || type == "L1I" || type == "D1I" || type == "S1I") return 30;
        if (type == "C6I" || type == "L6I" || type == "D6I" || type == "S6I") return 40;
        if (type == "C7I" || type == "L7I" || type == "D7I" || type == "S7I") return 20;
    }
    return 1;
}

// 解析 RINEX 头文件中的 GLONASS 频点号映射。
static void parseGloSlotFreqLine(const std::string &line) {
    // GLONASS 频点号在 RINEX 里需要单独读，这里只做基础映射缓存。
    std::istringstream iss(safeSubstr(line, 0, 60));
    std::string tok;
    if (!(iss >> tok)) return;
    while (iss >> tok) {
        if (tok.size() < 2 || tok[0] != 'R') continue;
        const int prn = parseInt(tok.substr(1), -1);
        int channel = 0;
        if (prn > 0 && (iss >> channel)) gGloFreqChannel[prn] = channel;
    }
}

// 读取 RINEX 观测文件，输出近似坐标、观测类型和逐历元观测值。
RinexObsFile readRinexObs(const std::string &path, int maxEpochs) {
    // 读取 RINEX 观测文件，尽量兼容普通 RINEX 和部分厂商格式差异。
    // 输出的 ObsEpoch 序列是后续 SPP/RTK 的原始观测来源。
    // maxEpochs 只限制流动站读取数量，便于现场把模式改成小数据量快速演示。
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open observation file: " + path);

    RinexObsFile file;
    bool isSeptentrio = false;  // detected from header comments
    std::string line;
    while (std::getline(in, line)) {
        const std::string label = safeSubstr(line, 60, 20);
        // Detect Septentrio RINEX format (uses space-separated "value LLI" tokens)
        // Check any header line for Septentrio markers
        if (!isSeptentrio &&
            (line.find("SEPT") != std::string::npos ||
             line.find("Septentrio") != std::string::npos ||
             line.find("ssrcrin") != std::string::npos)) {
            isSeptentrio = true;
        }
        if (label.find("APPROX POSITION XYZ") != std::string::npos) {
            file.approxPos.x = parseDouble(safeSubstr(line, 0, 14));
            file.approxPos.y = parseDouble(safeSubstr(line, 14, 14));
            file.approxPos.z = parseDouble(safeSubstr(line, 28, 14));
        } else if (label.find("SYS / # / OBS TYPES") != std::string::npos && !line.empty()) {
            const System sys = systemFromChar(line[0]);
            if (sys == System::GPS || sys == System::BDS || sys == System::GLO) {
                int n = parseInt(safeSubstr(line, 3, 3));
                auto &types = file.obsTypes[sys];
                int read = 0;
                while (read < n) {
                    for (int i = 0; i < 13 && read < n; ++i) {
                        std::string t = trim(safeSubstr(line, 7 + i * 4, 3));
                        if (!t.empty()) {
                            types.push_back(t);
                            ++read;
                        }
                    }
                    if (read < n) std::getline(in, line);
                }
            }
        } else if (label.find("GLONASS SLOT / FRQ #") != std::string::npos) {
            parseGloSlotFreqLine(line);
        } else if (line.find("END OF HEADER") != std::string::npos) {
            break;
        }
    }

    while (std::getline(in, line)) {
        if (line.empty() || line[0] != '>') continue;
        std::istringstream iss(line.substr(1));
        int y, mo, d, h, mi, flag, nsat;
        double sec;
        if (!(iss >> y >> mo >> d >> h >> mi >> sec >> flag >> nsat)) continue;
        if (flag != 0 && flag != 1) {
            for (int i = 0; i < nsat; ++i) std::getline(in, line);
            continue;
        }
        ObsEpoch epoch;
        epoch.time = epochToGps(y, mo, d, h, mi, sec);
        for (int i = 0; i < nsat; ++i) {
            if (!std::getline(in, line) || line.size() < 3) continue;
            SatId sat{systemFromChar(line[0]), parseInt(safeSubstr(line, 1, 2))};
            if (sat.sys != System::GPS && sat.sys != System::BDS && sat.sys != System::GLO) continue;
            auto it = file.obsTypes.find(sat.sys);
            if (it == file.obsTypes.end()) continue;
            ObsValue ov;
            int codePri[3] = {};
            int phasePri[3] = {};
            int dopPri[3] = {};
            int snrPri[3] = {};
            const auto &types = it->second;

            if (false && isSeptentrio) {
                // Disabled: Septentrio RINEX is still fixed 16-char OBS fields.
                // Token parsing loses blank fields and shifts C6I/L6I into C7I/L7I.
                std::vector<std::string> tokens;
                {
                    std::istringstream tstream(line.substr(3));
                    std::string t;
                    while (tstream >> t) tokens.push_back(t);
                }
                size_t tokIdx = 0;
                for (size_t k = 0; k < types.size(); ++k) {
                    while (tokIdx < tokens.size() && tokens[tokIdx].size() == 1 &&
                           std::isdigit(static_cast<unsigned char>(tokens[tokIdx][0]))) {
                        ++tokIdx;
                    }
                    if (tokIdx >= tokens.size()) break;
                    const double val = parseDouble(tokens[tokIdx]);
                    ++tokIdx;
                    if (val == 0.0) continue;
                    const int fi = freqIndex(sat.sys, types[k]);
                    if (fi < 0 || fi >= 3) continue;
                    const int pri = obsPriority(sat.sys, types[k]);
                    if (types[k][0] == 'C' && pri > codePri[fi]) {
                        ov.code[fi] = val; ov.hasCode[fi] = true; codePri[fi] = pri;
                    } else if (types[k][0] == 'L' && pri > phasePri[fi]) {
                        ov.phase[fi] = val; ov.hasPhase[fi] = true; phasePri[fi] = pri;
                        if (tokIdx < tokens.size() && tokens[tokIdx].size() == 1 &&
                            std::isdigit(static_cast<unsigned char>(tokens[tokIdx][0]))) {
                            ov.lli[fi] = tokens[tokIdx][0] - '0'; ++tokIdx;
                        }
                    } else if (types[k][0] == 'D' && pri > dopPri[fi]) {
                        ov.doppler[fi] = val; ov.hasDoppler[fi] = true; dopPri[fi] = pri;
                    } else if (types[k][0] == 'S' && pri > snrPri[fi]) {
                        ov.snr[fi] = val; snrPri[fi] = pri;
                    }
                }
            } else {
                // Standard RINEX 16-character fixed-field format
                // 逐卫星逐观测类型提取码伪距、载波相位、多普勒和 SNR。
                for (size_t k = 0; k < types.size(); ++k) {
                    std::string field = safeSubstr(line, 3 + k * 16, 14);
                    bool hasNumeric = false;
                    for (char c : field) {
                        if (c == '.' || c == 'D' || c == 'd' || c == 'E' || c == 'e') {
                            hasNumeric = true;
                            break;
                        }
                    }
                    if (!hasNumeric) {
                        continue;
                    }
                    const double val = parseDouble(field);
                    if (val == 0.0) continue;
                    const int fi = freqIndex(sat.sys, types[k]);
                    if (fi < 0 || fi >= 3) continue;
                    const int pri = obsPriority(sat.sys, types[k]);
                    if (types[k][0] == 'C' && pri > codePri[fi]) {
                        ov.code[fi] = val; ov.hasCode[fi] = true; codePri[fi] = pri;
                    } else if (types[k][0] == 'L' && pri > phasePri[fi]) {
                        ov.phase[fi] = val; ov.hasPhase[fi] = true; phasePri[fi] = pri;
                        // LLI is at position 14 within the 16-char field (the char AFTER our 14-char field)
                        if (line.size() > 3 + k * 16 + 14) {
                            const char lli = line[3 + k * 16 + 14];
                            if (lli >= '0' && lli <= '9') ov.lli[fi] = lli - '0';
                        }
                    } else if (types[k][0] == 'D' && pri > dopPri[fi]) {
                        ov.doppler[fi] = val; ov.hasDoppler[fi] = true; dopPri[fi] = pri;
                    } else if (types[k][0] == 'S' && pri > snrPri[fi]) {
                        ov.snr[fi] = val; snrPri[fi] = pri;
                    }
                }
            }
            epoch.obs[sat] = ov;
        }
        file.epochs.push_back(epoch);
        if (maxEpochs > 0 && static_cast<int>(file.epochs.size()) >= maxEpochs) break;
    }
    return file;
}

// 读取一个 GPS/BDS 广播星历块的后 7 行数值字段。
static bool readNavBlock(std::ifstream &in, std::vector<double> &v) {
    v.clear();
    std::string line;
    for (int i = 0; i < 7; ++i) {
        if (!std::getline(in, line)) return false;
        for (int j = 0; j < 4; ++j) v.push_back(parseDouble(safeSubstr(line, 4 + j * 19, 19)));
    }
    return true;
}

// 读取 RINEX 广播星历文件，输出 GPS/BDS/GLO 星历集合。
NavData readRinexNav(const std::string &path) {
    // 读取广播星历，并统一组织成 Eph 列表。
    // RTK 构造观测方程前，会通过 satelliteState 查询卫星 ECEF 坐标和钟差。
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open navigation file: " + path);
    std::string line;
    while (std::getline(in, line)) {
        if (line.find("END OF HEADER") != std::string::npos) break;
    }
    NavData nav;
    while (std::getline(in, line)) {
        if (line.size() < 23) continue;
        SatId sat{systemFromChar(line[0]), parseInt(safeSubstr(line, 1, 2))};
        int y = parseInt(safeSubstr(line, 4, 4));
        int mo = parseInt(safeSubstr(line, 9, 2));
        int d = parseInt(safeSubstr(line, 12, 2));
        int h = parseInt(safeSubstr(line, 15, 2));
        int mi = parseInt(safeSubstr(line, 18, 2));
        double sec = parseDouble(safeSubstr(line, 21, 2));
        if (sat.sys == System::GLO) {
            Eph eph;
            eph.sat = sat;
            eph.toc = epochToGps(y, mo, d, h, mi, sec);
            eph.af0 = -parseDouble(safeSubstr(line, 23, 19));
            eph.af1 = parseDouble(safeSubstr(line, 42, 19));
            eph.toe = parseDouble(safeSubstr(line, 61, 19));
            std::vector<double> v;
            v.clear();
            for (int i = 0; i < 3; ++i) {
                if (!std::getline(in, line)) break;
                for (int j = 0; j < 4; ++j) v.push_back(parseDouble(safeSubstr(line, 4 + j * 19, 19)));
            }
            if (v.size() < 12) break;
            eph.gloPosKm = {v[0], v[4], v[8]};
            eph.gloVelKmS = {v[1], v[5], v[9]};
            eph.gloAccKmS2 = {v[2], v[6], v[10]};
            eph.week = eph.toc.week;
            eph.toeTime = eph.toc;
            nav.ephs.push_back(eph);
            continue;
        }
        if (sat.sys != System::GPS && sat.sys != System::BDS) {
            const int skipLines = (line[0] == 'R' || line[0] == 'S') ? 3 : 7;
            for (int i = 0; i < skipLines; ++i) {
                if (!std::getline(in, line)) break;
            }
            continue;
        }
        Eph eph;
        eph.sat = sat;
        eph.toc = epochToGps(y, mo, d, h, mi, sec);
        eph.af0 = parseDouble(safeSubstr(line, 23, 19));
        eph.af1 = parseDouble(safeSubstr(line, 42, 19));
        eph.af2 = parseDouble(safeSubstr(line, 61, 19));
        std::vector<double> v;
        if (!readNavBlock(in, v) || v.size() < 28) break;
        eph.iode = v[0]; eph.crs = v[1]; eph.dn = v[2]; eph.m0 = v[3];
        eph.cuc = v[4]; eph.e = v[5]; eph.cus = v[6]; eph.sqrtA = v[7];
        eph.toe = v[8]; eph.cic = v[9]; eph.omg0 = v[10]; eph.cis = v[11];
        eph.i0 = v[12]; eph.crc = v[13]; eph.omg = v[14]; eph.omgDot = v[15];
        eph.idot = v[16]; eph.week = static_cast<int>(v[18]);
        if (eph.week <= 0) eph.week = eph.toc.week;
        eph.toes = eph.toe;
        if (sat.sys == System::BDS) {
            // RTKLIB converts BDT to GPST by adding 14 s; BDT week starts at 2006-01-01.
            const GTime bdt0 = epochToGps(2006, 1, 1, 0, 0, 0.0);
            double full = bdt0.week * 604800.0 + bdt0.tow + eph.week * 604800.0 + eph.toe + 14.0;
            eph.toeTime.week = static_cast<int>(std::floor(full / 604800.0));
            eph.toeTime.tow = full - eph.toeTime.week * 604800.0;
            eph.toeTime = normalizeTime(eph.toeTime);
            eph.toc.tow += 14.0;
            eph.toc = normalizeTime(eph.toc);
        } else {
            eph.toeTime = normalizeTime({eph.week, eph.toe});
        }
        nav.ephs.push_back(eph);
    }
    return nav;
}

// 从广播星历列表中选择距离当前时刻最近的一组星历。
static const Eph *selectEph(const NavData &nav, const SatId &sat, const GTime &t) {
    const Eph *best = nullptr;
    double bestDt = 1e99;
    for (const auto &e : nav.ephs) {
        if (e.sat.sys != sat.sys || e.sat.prn != sat.prn) continue;
        const double dt = std::abs(timediff(t, e.toeTime));
        if (dt < bestDt) {
            bestDt = dt;
            best = &e;
        }
    }
    return best;
}

static Vec3 propagateGloPosition(const Eph &e, double dt);

// 由广播星历计算卫星 ECEF 坐标和卫星钟差。
bool satelliteState(const NavData &nav, const SatId &sat, const GTime &t, SatState &state) {
    const Eph *e = selectEph(nav, sat, t);
    if (!e) return false;
    if (sat.sys == System::GLO) {
        const double dt = timediff(t, e->toeTime);
        state.pos = propagateGloPosition(*e, dt);
        state.clk = e->af0 + e->af1 * dt;
        state.valid = true;
        return true;
    }
    if (e->sqrtA <= 0.0) return false;
    const double mu = sat.sys == System::BDS ? kMuBds : kMuGps;
    const double a = e->sqrtA * e->sqrtA;
    const double n0 = std::sqrt(mu / (a * a * a));
    const double tk = timediff(t, e->toeTime);
    const double n = n0 + e->dn;
    double m = e->m0 + n * tk;
    double ek = m;
    for (int i = 0; i < 12; ++i) ek = m + e->e * std::sin(ek);
    const double sinv = std::sqrt(1.0 - e->e * e->e) * std::sin(ek) / (1.0 - e->e * std::cos(ek));
    const double cosv = (std::cos(ek) - e->e) / (1.0 - e->e * std::cos(ek));
    const double v = std::atan2(sinv, cosv);
    const double phi = v + e->omg;
    const double du = e->cus * std::sin(2.0 * phi) + e->cuc * std::cos(2.0 * phi);
    const double dr = e->crs * std::sin(2.0 * phi) + e->crc * std::cos(2.0 * phi);
    const double di = e->cis * std::sin(2.0 * phi) + e->cic * std::cos(2.0 * phi);
    const double u = phi + du;
    const double r = a * (1.0 - e->e * std::cos(ek)) + dr;
    const double inc = e->i0 + di + e->idot * tk;
    const double xop = r * std::cos(u);
    const double yop = r * std::sin(u);
    int prn = sat.prn;
    if (sat.sys == System::BDS && prn <= 5) {
        const double omg = e->omg0 + e->omgDot * tk - kOmegaEarth * e->toes;
        const double xg = xop * std::cos(omg) - yop * std::cos(inc) * std::sin(omg);
        const double yg = xop * std::sin(omg) + yop * std::cos(inc) * std::cos(omg);
        const double zg = yop * std::sin(inc);
        const double sino = std::sin(kOmegaEarth * tk);
        const double coso = std::cos(kOmegaEarth * tk);
        const double sin5 = -0.0871557427476582;
        const double cos5 = 0.9961946980917456;
        state.pos.x = xg * coso + yg * sino * cos5 + zg * sino * sin5;
        state.pos.y = -xg * sino + yg * coso * cos5 + zg * coso * sin5;
        state.pos.z = -yg * sin5 + zg * cos5;
    } else {
        const double omg = e->omg0 + (e->omgDot - kOmegaEarth) * tk - kOmegaEarth * e->toes;
        state.pos.x = xop * std::cos(omg) - yop * std::cos(inc) * std::sin(omg);
        state.pos.y = xop * std::sin(omg) + yop * std::cos(inc) * std::cos(omg);
        state.pos.z = yop * std::sin(inc);
    }
    const double tc = timediff(t, e->toc);
    state.clk = e->af0 + e->af1 * tc + e->af2 * tc * tc + kFRel * e->e * e->sqrtA * std::sin(ek);
    state.valid = true;
    return true;
}

// 计算接收机到卫星的几何距离和视线方向单位向量。
static bool geometric(const Vec3 &rx, const Vec3 &sat, double &rho, Vec3 &los) {
    Vec3 dr = sat - rx;
    rho = norm(dr);
    if (rho <= 0.0) return false;
    los = dr / rho;
    return true;
}

// 地球自转改正：信号传播时间内地球旋转会改变卫星 ECEF 表达。
static Vec3 earthRotate(const Vec3 &sat, double tau) {
    const double a = kOmegaEarth * tau;
    const double ca = std::cos(a);
    const double sa = std::sin(a);
    return {ca * sat.x + sa * sat.y, -sa * sat.x + ca * sat.y, sat.z};
}

// GLONASS 数值积分所需的加速度模型，包含 J2 摄动和外加加速度。
static void gloAccel(const Vec3 &r, const Vec3 &extAcc, Vec3 &a) {
    constexpr double mu = 3.9860044e14;
    constexpr double ae = 6378136.0;
    constexpr double j2 = 1.0826257e-3;
    const double r2 = dot(r, r);
    const double rn = std::sqrt(r2);
    const double z2 = r.z * r.z;
    const double c = -mu / (rn * rn * rn);
    const double f = 1.5 * j2 * mu * ae * ae / (rn * rn * rn * rn * rn);
    const double zf = 5.0 * z2 / r2;
    a.x = (c - f * (1.0 - zf)) * r.x + extAcc.x + kOmegaEarth * kOmegaEarth * r.x;
    a.y = (c - f * (1.0 - zf)) * r.y + extAcc.y + kOmegaEarth * kOmegaEarth * r.y;
    a.z = (c - f * (3.0 - zf)) * r.z + extAcc.z;
}

// GLONASS 轨道积分微分方程，返回位置和速度导数。
static void gloDerivative(const Vec3 &r, const Vec3 &v, const Vec3 &extAcc, Vec3 &dr, Vec3 &dv) {
    dr = v;
    gloAccel(r, extAcc, dv);
    dv.x += 2.0 * kOmegaEarth * v.y;
    dv.y -= 2.0 * kOmegaEarth * v.x;
}

// GLONASS 轨道四阶 Runge-Kutta 积分一步。
static void rk4GloStep(Vec3 &r, Vec3 &v, const Vec3 &extAcc, double h) {
    Vec3 k1r, k1v, k2r, k2v, k3r, k3v, k4r, k4v;
    gloDerivative(r, v, extAcc, k1r, k1v);
    gloDerivative(r + k1r * (0.5 * h), v + k1v * (0.5 * h), extAcc, k2r, k2v);
    gloDerivative(r + k2r * (0.5 * h), v + k2v * (0.5 * h), extAcc, k3r, k3v);
    gloDerivative(r + k3r * h, v + k3v * h, extAcc, k4r, k4v);
    r = r + (k1r + k2r * 2.0 + k3r * 2.0 + k4r) * (h / 6.0);
    v = v + (k1v + k2v * 2.0 + k3v * 2.0 + k4v) * (h / 6.0);
}

// 从 GLONASS 星历历元向目标时刻积分得到卫星位置。
static Vec3 propagateGloPosition(const Eph &e, double dt) {
    Vec3 r = e.gloPosKm * 1000.0;
    Vec3 v = e.gloVelKmS * 1000.0;
    Vec3 acc = e.gloAccKmS2 * 1000.0;
    const double dir = dt >= 0.0 ? 1.0 : -1.0;
    double rem = std::abs(dt);
    while (rem > 1e-9) {
        const double h = dir * std::min(60.0, rem);
        rk4GloStep(r, v, acc, h);
        rem -= std::abs(h);
    }
    return r;
}

// 用观测伪距估计发射时刻，并计算带卫星钟差改正的几何量。
static bool observedGeometry(const NavData &nav, const SatId &sat, const GTime &recvTime, double code,
                             const Vec3 &rx, double &rho, Vec3 &los) {
    if (code <= 0.0) return false;
    GTime tx = recvTime;
    tx.tow -= code / kClight;
    tx = normalizeTime(tx);
    SatState clkState;
    if (satelliteState(nav, sat, tx, clkState)) {
        tx.tow -= clkState.clk;
        tx = normalizeTime(tx);
    }
    SatState ss;
    if (!satelliteState(nav, sat, tx, ss)) return false;
    ss.pos = earthRotate(ss.pos, code / kClight);
    if (!geometric(rx, ss.pos, rho, los)) return false;
    rho -= kClight * ss.clk;
    return true;
}

// 单点定位：用伪距最小二乘估计接收机坐标和系统钟差。
bool sppSolve(const ObsEpoch &epoch, const NavData &nav, Vec3 initial, GnssSolution &sol) {
    Vec3 x = norm(initial) > 1.0 ? initial : Vec3{-2267639.0, 5009266.0, 3221240.0};
    double cbGps = 0.0, cbBds = 0.0;
    int used = 0;
    Matrix q;
    for (int iter = 0; iter < 8; ++iter) {
        std::vector<std::array<double, 5>> rows;
        std::vector<double> v;
        bool hasBds = false;
        for (const auto &kv : epoch.obs) {
            const SatId sat = kv.first;
            if (sat.sys != System::GPS && sat.sys != System::BDS) continue;
            const ObsValue &o = kv.second;
            if (!o.hasCode[0]) continue;
            GTime tx = epoch.time;
            tx.tow -= o.code[0] / kClight;
            tx = normalizeTime(tx);
            SatState clkState;
            if (satelliteState(nav, sat, tx, clkState)) {
                tx.tow -= clkState.clk;
                tx = normalizeTime(tx);
            }
            SatState ss;
            if (!satelliteState(nav, sat, tx, ss)) continue;
            ss.pos = earthRotate(ss.pos, o.code[0] / kClight);
            double rho;
            Vec3 e;
            if (!geometric(x, ss.pos, rho, e)) continue;
            const double cb = sat.sys == System::GPS ? cbGps : cbBds;
            const double pred = rho + cb - kClight * ss.clk;
            const double res = o.code[0] - pred;
            rows.push_back({-e.x, -e.y, -e.z, sat.sys == System::GPS ? 1.0 : 0.0,
                            sat.sys == System::BDS ? 1.0 : 0.0});
            v.push_back(res);
            if (sat.sys == System::BDS) hasBds = true;
        }
        used = static_cast<int>(v.size());
        const int nx = hasBds ? 5 : 4;
        if (used < nx) return false;
        Matrix h(used, nx);
        for (int i = 0; i < used; ++i) for (int j = 0; j < nx; ++j) h(i, j) = rows[i][j];
        std::vector<double> dx;
        if (!leastSquares(h, v, dx, &q)) return false;
        x.x += dx[0]; x.y += dx[1]; x.z += dx[2]; cbGps += dx[3]; cbBds += hasBds ? dx[4] : 0.0;
        if (std::sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]) < 1e-4) break;
    }
    sol.time = epoch.time;
    sol.ecef = x;
    sol.llh = ecefToLlh(x);
    sol.clock = cbGps / kClight;
    sol.ns = used;
    sol.q = 5;
    sol.mode = "SPP";
    sol.stdN = sol.stdE = 3.0;
    sol.stdD = 8.0;
    return true;
}

static double carrierWavelength(const SatId &sat, int fi);

// 使用多普勒观测估计 ECEF 速度，为 RTK 和松组合提供速度量测。
static bool dopplerVelocitySolve(const ObsEpoch &epoch, const NavData &nav, const Vec3 &rx,
                                 Vec3 &velEcef, int &used, double &rms) {
    used = 0;
    rms = 0.0;
    std::vector<std::array<double, 5>> rows;
    std::vector<double> v;
    bool hasBds = false;
    for (const auto &kv : epoch.obs) {
        const SatId sat = kv.first;
        if (sat.sys != System::GPS && sat.sys != System::BDS) continue;
        const ObsValue &o = kv.second;
        if (!o.hasCode[0] || !o.hasDoppler[0]) continue;
        GTime tx = epoch.time;
        tx.tow -= o.code[0] / kClight;
        tx = normalizeTime(tx);
        SatState clkState;
        if (satelliteState(nav, sat, tx, clkState)) {
            tx.tow -= clkState.clk;
            tx = normalizeTime(tx);
        }
        SatState ss1, ss0;
        if (!satelliteState(nav, sat, tx, ss1)) continue;
        GTime tx0 = tx;
        tx0.tow -= 1.0;
        tx0 = normalizeTime(tx0);
        if (!satelliteState(nav, sat, tx0, ss0)) continue;
        const double tau = o.code[0] / kClight;
        ss1.pos = earthRotate(ss1.pos, tau);
        ss0.pos = earthRotate(ss0.pos, tau);
        double rho = 0.0;
        Vec3 los;
        if (!geometric(rx, ss1.pos, rho, los)) continue;
        const Vec3 satVel = ss1.pos - ss0.pos;
        const double lam = carrierWavelength(sat, 0);
        if (lam <= 0.0) continue;
        const double rangeRate = -o.doppler[0] * lam;
        rows.push_back({los.x, los.y, los.z, sat.sys == System::GPS ? 1.0 : 0.0,
                        sat.sys == System::BDS ? 1.0 : 0.0});
        v.push_back(dot(los, satVel) - rangeRate);
        if (sat.sys == System::BDS) hasBds = true;
    }
    used = static_cast<int>(v.size());
    const int nx = hasBds ? 5 : 4;
    if (used < nx) return false;
    Matrix h(used, nx);
    for (int i = 0; i < used; ++i) {
        for (int j = 0; j < nx; ++j) h(i, j) = rows[i][j];
    }
    std::vector<double> dx;
    if (!leastSquares(h, v, dx, nullptr)) return false;
    velEcef = {dx[0], dx[1], dx[2]};
    double sum2 = 0.0;
    for (int i = 0; i < used; ++i) {
        double fit = 0.0;
        for (int j = 0; j < nx; ++j) fit += h(i, j) * dx[j];
        sum2 += sqr(v[i] - fit);
    }
    rms = std::sqrt(sum2 / std::max(1, used - nx));
    return std::isfinite(velEcef.x) && std::isfinite(velEcef.y) && std::isfinite(velEcef.z) && rms < 3.0;
}

// 在另一站观测序列中寻找与当前历元时间一致的观测，用于基准站/流动站配对。
static const ObsEpoch *findEpochNear(const std::vector<ObsEpoch> &epochs, const GTime &t, size_t &idx) {
    while (idx + 1 < epochs.size() && timediff(epochs[idx].time, t) < -1e-4) ++idx;
    if (idx < epochs.size() && std::abs(timediff(epochs[idx].time, t)) < 1e-4) return &epochs[idx];
    if (idx > 0 && std::abs(timediff(epochs[idx - 1].time, t)) < 1e-4) return &epochs[idx - 1];
    return nullptr;
}

// 根据卫星系统和频点返回载波波长，载波相位米化和模糊度固定都依赖它。
static double carrierWavelength(const SatId &sat, int fi) {
    if (sat.sys == System::GPS) {
        if (fi == 0) return kClight / 1575.42e6;
        if (fi == 1) return kClight / 1227.60e6;
    } else if (sat.sys == System::GLO) {
        const int k = gGloFreqChannel.count(sat.prn) ? gGloFreqChannel[sat.prn] : 0;
        if (fi == 0) return kClight / ((1602.0 + 0.5625 * k) * 1e6);
        if (fi == 1) return kClight / ((1246.0 + 0.4375 * k) * 1e6);
    } else if (sat.sys == System::BDS) {
        if (fi == 0) return kClight / 1561.098e6;
        if (fi == 1) return kClight / 1268.52e6;
    }
    return kClight / 1575.42e6;
}

struct DdSat {
    SatId sat;
    Vec3 los;
    double elev = 0.0;
    double roverCode[2] = {};
    double baseCode[2] = {};
    double roverPhaseM[2] = {};
    double basePhaseM[2] = {};
    double sdCode[2] = {};
    double sdPhaseM[2] = {};
    int roverLli[2] = {};
    int baseLli[2] = {};
    bool ok[2] = {};
};

struct RtkFilterState {
    bool initialized = false;
    GTime time;
    std::vector<double> x;        // [pos(3), vel(3), clkGps(1), clkBds(1), amb_sd...]
    Matrix p;
    std::map<std::string, int> ambIndex;  // key->state index
    std::map<int, int> ambLock;           // state index->lock count
    std::map<int, int> ambSlipFlag;       // state index->slip detected in this epoch
    std::map<std::string, int> lastRoverLli;
    std::map<std::string, int> lastBaseLli;
    std::map<std::string, double> lastRoverGf;
    std::map<std::string, double> lastBaseGf;
    std::map<std::string, double> lastRoverMw;
    std::map<std::string, double> lastBaseMw;
    std::map<std::pair<System, int>, int> refPrn;
    Vec3 basePosition;                    // base station ECEF (fixed)
    bool acceptedSolution = false;
    bool hasLastAccepted = false;
    GTime lastAcceptedTime;
    Vec3 lastAcceptedPos;
    int fixedHoldAge = 1000000;
};

static constexpr int kClkGpsIdx = kRtkPosStates + kRtkVelStates;   // = 6
static constexpr int kClkBdsIdx = kClkGpsIdx + 1;                   // = 7
static constexpr int kRtkMinState = kClkBdsIdx + 1;                 // = 8

// 单差模糊度状态键：卫星号 + 频点。
static std::string ambStateKey(const SatId &sat, int f) {
    return satToString(sat) + ":" + std::to_string(f);
}

// 双差模糊度诊断键：目标星、参考星和频点共同确定。
static std::string ddAmbStateKey(const SatId &sat, const SatId &ref, int f) {
    return satToString(sat) + "-" + satToString(ref) + ":" + std::to_string(f);
}

// 初始化 RTK Kalman 滤波器，状态包含位置、速度、系统钟差和后续动态添加的模糊度。
static void initRtkFilter(RtkFilterState &flt, const GTime &time, const Vec3 &pos, const Vec3 &basePos) {
    flt.initialized = true;
    flt.time = time;
    flt.basePosition = basePos;
    // State: pos(3), vel(3), clkGps(1), clkBds(1)
    flt.x.resize(kRtkMinState);
    flt.x[0] = pos.x; flt.x[1] = pos.y; flt.x[2] = pos.z;
    for (int i = 3; i < kRtkMinState; ++i) flt.x[i] = 0.0;
    flt.p = Matrix(kRtkMinState, kRtkMinState);
    flt.p(0, 0) = flt.p(1, 1) = flt.p(2, 2) = 10000.0;
    flt.p(3, 3) = flt.p(4, 4) = flt.p(5, 5) = 100.0;
    flt.p(kClkGpsIdx, kClkGpsIdx) = 1e6;   // GPS clock: large initial uncertainty
    flt.p(kClkBdsIdx, kClkBdsIdx) = 1e6;   // BDS clock: large initial uncertainty
    flt.ambIndex.clear();
    flt.ambLock.clear();
    flt.ambSlipFlag.clear();
    flt.lastRoverLli.clear();
    flt.lastBaseLli.clear();
    flt.lastRoverGf.clear();
    flt.lastBaseGf.clear();
    flt.lastRoverMw.clear();
    flt.lastBaseMw.clear();
    flt.refPrn.clear();
    flt.acceptedSolution = false;
    flt.hasLastAccepted = false;
    flt.lastAcceptedTime = time;
    flt.lastAcceptedPos = pos;
    flt.fixedHoldAge = 1000000;
}

// 确保某颗卫星某频点的单差模糊度状态存在；不存在时动态扩展状态向量和协方差。
static int ensureAmbState(RtkFilterState &flt, const std::string &key, double initialCycles) {
    auto it = flt.ambIndex.find(key);
    if (it != flt.ambIndex.end()) return it->second;

    const int oldN = static_cast<int>(flt.x.size());
    Matrix np(oldN + 1, oldN + 1);
    for (int r = 0; r < oldN; ++r) {
        for (int c = 0; c < oldN; ++c) np(r, c) = flt.p(r, c);
    }
    np(oldN, oldN) = 100.0;
    flt.p = np;
    flt.x.push_back(initialCycles);
    flt.ambIndex[key] = oldN;
    flt.ambLock[oldN] = 0;
    return oldN;
}

// 周跳或质量异常时重置单个模糊度状态，清空相关协方差并重新给初值。
static void resetAmbState(RtkFilterState &flt, int idx, double valueCycles, double variance) {
    if (idx < 0 || idx >= static_cast<int>(flt.x.size())) return;
    flt.x[idx] = valueCycles;
    for (int i = 0; i < flt.p.rows(); ++i) {
        flt.p(idx, i) = 0.0;
        flt.p(i, idx) = 0.0;
    }
    flt.p(idx, idx) = variance;
    flt.ambLock[idx] = 0;
}

// 清空所有模糊度状态，仅保留位置、速度和系统钟差，用于大失锁后的重新收敛。
static void resetRtkAmbiguities(RtkFilterState &flt) {
    if (static_cast<int>(flt.x.size()) <= kRtkMinState) return;
    flt.x.resize(kRtkMinState);
    Matrix np(kRtkMinState, kRtkMinState);
    for (int r = 0; r < kRtkMinState; ++r) {
        for (int c = 0; c < kRtkMinState; ++c) np(r, c) = flt.p(r, c);
    }
    flt.p = np;
    flt.ambIndex.clear();
    flt.ambLock.clear();
    flt.ambSlipFlag.clear();
    flt.refPrn.clear();
    flt.fixedHoldAge = 1000000;
}

// 对已有单差模糊度做系统内共同偏移校正，使相位-伪距初值更一致。
static void alignSdAmbiguityOffsets(RtkFilterState &flt, const std::vector<DdSat> &sats) {
    for (int f = 0; f < 2; ++f) {
        for (System sys : {System::GPS, System::BDS, System::GLO}) {
            double offset = 0.0;
            int nExisting = 0;
            std::vector<std::pair<int, double>> available;
            for (const auto &s : sats) {
                if (!s.ok[f] || s.sat.sys != sys) continue;
                const double lam = carrierWavelength(s.sat, f);
                if (lam <= 0.0 || s.sdCode[f] == 0.0 || s.sdPhaseM[f] == 0.0) continue;
                const double phaseCodeBias = (s.sdPhaseM[f] - s.sdCode[f]) / lam;
                const std::string key = ambStateKey(s.sat, f);
                auto it = flt.ambIndex.find(key);
                if (it == flt.ambIndex.end() || it->second < 0 || it->second >= static_cast<int>(flt.x.size())) {
                    continue;
                }
                available.push_back({it->second, phaseCodeBias});
                offset += phaseCodeBias - flt.x[it->second];
                ++nExisting;
            }
            if (nExisting <= 0) continue;
            offset /= static_cast<double>(nExisting);
            for (const auto &p : available) flt.x[p.first] += offset;
        }
    }
}

// RTK Kalman 预测步：位置按速度外推，钟差和模糊度采用常值模型。
static void predictRtkFilter(RtkFilterState &flt, const GTime &time, const RtkOptions &opt) {
    if (!flt.initialized) return;
    double dt = timediff(time, flt.time);
    double qdt = std::abs(dt);
    if (qdt <= 0.0 || qdt > 30.0) {
        dt = (dt < 0.0) ? -1.0 : 1.0;
        qdt = 1.0;
    }
    const int n = static_cast<int>(flt.x.size());
    Matrix f = Matrix::identity(n);
    if (n >= kRtkAmbStart) {
        for (int i = 0; i < 3; ++i) {
            flt.x[i] += flt.x[i + 3] * dt;
            f(i, i + 3) = dt;
        }
    }
    // Clocks and ambiguities: constant model.
    flt.p = f * flt.p * f.transpose();

    // Process noise
    const double qpos = sqr(opt.processPos) * qdt;
    const double qvel = sqr(opt.processAccel) * qdt;
    const double qclk = sqr(10.0) * qdt;  // 10 m/s clock drift noise
    for (int i = 0; i < 3 && i < n; ++i) flt.p(i, i) += qpos;
    for (int i = 3; i < kClkGpsIdx && i < n; ++i) flt.p(i, i) += qvel;
    // Clock process noise
    if (kClkGpsIdx < n) flt.p(kClkGpsIdx, kClkGpsIdx) += qclk;
    if (kClkBdsIdx < n) flt.p(kClkBdsIdx, kClkBdsIdx) += qclk;
    // Ambiguity process noise (very small)
    for (int i = kRtkMinState; i < n; ++i) flt.p(i, i) += 1e-8 * qdt;

    flt.time = time;
}

// Scalar Kalman update: single observation at a time
// h: row of H matrix (length n)
// v: innovation scalar = obs - predicted
// r: measurement variance
// 单观测标量 Kalman 更新，用于伪量测、保持约束和部分 RTK 观测。
static bool scalarKalmanUpdate(RtkFilterState &flt, const std::vector<double> &h, double v, double r) {
    const int n = static_cast<int>(flt.x.size());
    if (static_cast<int>(h.size()) != n || r <= 0.0) return false;

    // Compute s = h*P*h' + r
    double s = r;
    std::vector<double> Ph(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) sum += flt.p(i, j) * h[j];
        Ph[i] = sum;
        s += sum * h[i];
    }
    if (s < 1e-12 || !std::isfinite(s)) return false;

    // Kalman gain k = P*h' / s
    std::vector<double> k(n);
    for (int i = 0; i < n; ++i) {
        k[i] = Ph[i] / s;
        if (!std::isfinite(k[i])) return false;
    }

    // State update: x += k * v
    for (int i = 0; i < n; ++i) flt.x[i] += k[i] * v;

    // Covariance update: P = P - K * H * P. Here Ph = P*H^T, and P is symmetric.
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double pij = flt.p(i, j) - k[i] * Ph[j];
            flt.p(i, j) = pij;
            if (i != j) flt.p(j, i) = pij;
        }
    }
    for (int i = 0; i < n; ++i) if (flt.p(i, i) < 1e-12) flt.p(i, i) = 1e-12;
    return true;
}

// 位置伪量测更新，用 SPP/先验结果约束 RTK 滤波器不发生大漂移。
static bool applyPositionPseudoMeasurement(RtkFilterState &flt, const Vec3 &pos, double sigma) {
    if (sigma <= 0.0 || !std::isfinite(sigma) ||
        !std::isfinite(pos.x) || !std::isfinite(pos.y) || !std::isfinite(pos.z)) {
        return false;
    }
    const int n = static_cast<int>(flt.x.size());
    if (n < 3) return false;
    const double r = sigma * sigma;
    for (int i = 0; i < 3; ++i) {
        std::vector<double> h(n, 0.0);
        h[i] = 1.0;
        const double target = i == 0 ? pos.x : (i == 1 ? pos.y : pos.z);
        if (!scalarKalmanUpdate(flt, h, target - flt.x[i], r)) return false;
    }
    return true;
}

// 速度伪量测更新，用多普勒速度辅助 RTK 状态中的速度分量。
static bool applyVelocityPseudoMeasurement(RtkFilterState &flt, const Vec3 &vel, double sigma) {
    if (sigma <= 0.0 || !std::isfinite(sigma) ||
        !std::isfinite(vel.x) || !std::isfinite(vel.y) || !std::isfinite(vel.z)) {
        return false;
    }
    const int n = static_cast<int>(flt.x.size());
    if (n < kClkGpsIdx) return false;
    const double r = sigma * sigma;
    for (int i = 0; i < 3; ++i) {
        std::vector<double> h(n, 0.0);
        h[i + 3] = 1.0;
        const double target = i == 0 ? vel.x : (i == 1 ? vel.y : vel.z);
        if (!scalarKalmanUpdate(flt, h, target - flt.x[i + 3], r)) return false;
    }
    return true;
}

// 双差模糊度保持：固定成功后用整数值反向约束浮点滤波器。
static bool applyDdAmbiguityHold(RtkFilterState &flt, int satIdx, int refIdx, double fixedCycles, double sigmaCycles) {
    if (satIdx < 0 || refIdx < 0 || satIdx >= static_cast<int>(flt.x.size()) ||
        refIdx >= static_cast<int>(flt.x.size()) || sigmaCycles <= 0.0) {
        return false;
    }
    std::vector<double> h(flt.x.size(), 0.0);
    h[satIdx] = 1.0;
    h[refIdx] = -1.0;
    const double innovation = fixedCycles - (flt.x[satIdx] - flt.x[refIdx]);
    return scalarKalmanUpdate(flt, h, innovation, sigmaCycles * sigmaCycles);
}

// 多观测批量 Kalman 更新，用于同一历元内成组双差观测。
static bool batchKalmanUpdate(RtkFilterState &flt, const Matrix &h, const std::vector<double> &v,
                              const Matrix &r) {
    const int n = static_cast<int>(flt.x.size());
    const int m = static_cast<int>(v.size());
    if (m <= 0 || h.rows() != m || h.cols() != n || r.rows() != m || r.cols() != m || flt.p.rows() != n) {
        return false;
    }
    Matrix hp = h * flt.p;
    Matrix s = hp * h.transpose() + r;
    Matrix sinv;
    try {
        sinv = inverse(s);
    } catch (...) {
        return false;
    }
    Matrix k = flt.p * h.transpose() * sinv;
    const std::vector<double> dx = k * v;
    for (int i = 0; i < n; ++i) flt.x[i] += dx[i];

    Matrix iKh = Matrix::identity(n) - k * h;
    flt.p = iKh * flt.p * iKh.transpose() + k * r * k.transpose();
    for (int i = 0; i < n; ++i) {
        if (flt.p(i, i) < 1e-12) flt.p(i, i) = 1e-12;
        for (int j = i + 1; j < n; ++j) {
            const double a = 0.5 * (flt.p(i, j) + flt.p(j, i));
            flt.p(i, j) = flt.p(j, i) = a;
        }
    }
    return true;
}

// 根据 RINEX LLI 标志判断载波是否发生周跳。
static bool lliIndicatesSlip(std::map<std::string, int> &last, const std::string &key, int lli) {
    auto it = last.find(key);
    const bool known = it != last.end();
    const int prev = known ? it->second : 0;
    last[key] = lli;
    if (!known) return false;
    if ((lli & 1) != 0 && (prev & 1) == 0) return true;
    return false;
}

// LAMBDA 内部使用的列主序索引函数。
static int cm(int row, int col, int n) { return row + col * n; }

// LAMBDA 搜索中的最近整数取整。
static double lambdaRound(double x) { return std::floor(x + 0.5); }

// LAMBDA LD 分解步骤，把浮点模糊度协方差分解为 L 和 D。
static int lambdaLd(int n, const Matrix &q, std::vector<double> &l, std::vector<double> &d) {
    // LAMBDA step 1: LD factorization.
    // q is the covariance of float ambiguities. LD factorization prepares it for
    // integer decorrelation and search.
    std::vector<double> a(n * n, 0.0);
    for (int r = 0; r < n; ++r) for (int c = 0; c < n; ++c) a[cm(r, c, n)] = q(r, c);
    for (int i = n - 1; i >= 0; --i) {
        d[i] = a[cm(i, i, n)];
        if (d[i] <= 0.0) return -1;
        const double sq = std::sqrt(d[i]);
        for (int j = 0; j <= i; ++j) l[cm(i, j, n)] = a[cm(i, j, n)] / sq;
        for (int j = 0; j <= i - 1; ++j) {
            for (int k = 0; k <= j; ++k) a[cm(j, k, n)] -= l[cm(i, k, n)] * l[cm(i, j, n)];
        }
        for (int j = 0; j <= i; ++j) l[cm(i, j, n)] /= l[cm(i, i, n)];
    }
    return 0;
}

// LAMBDA Gauss 变换，用于降低模糊度之间的相关性。
static void lambdaGauss(int n, std::vector<double> &l, std::vector<double> &z, int i, int j) {
    const int mu = static_cast<int>(lambdaRound(l[cm(i, j, n)]));
    if (mu == 0) return;
    for (int k = i; k < n; ++k) l[cm(k, j, n)] -= static_cast<double>(mu) * l[cm(k, i, n)];
    for (int k = 0; k < n; ++k) z[cm(k, j, n)] -= static_cast<double>(mu) * z[cm(k, i, n)];
}

// LAMBDA 置换步骤，交换相邻模糊度维度以改善搜索空间形状。
static void lambdaPerm(int n, std::vector<double> &l, std::vector<double> &d, int j, double del,
                       std::vector<double> &z) {
    const double eta = d[j] / del;
    const double lam = d[j + 1] * l[cm(j + 1, j, n)] / del;
    d[j] = eta * d[j + 1];
    d[j + 1] = del;
    for (int k = 0; k <= j - 1; ++k) {
        const double a0 = l[cm(j, k, n)];
        const double a1 = l[cm(j + 1, k, n)];
        l[cm(j, k, n)] = -l[cm(j + 1, j, n)] * a0 + a1;
        l[cm(j + 1, k, n)] = eta * a0 + lam * a1;
    }
    l[cm(j + 1, j, n)] = lam;
    for (int k = j + 2; k < n; ++k) std::swap(l[cm(k, j, n)], l[cm(k, j + 1, n)]);
    for (int k = 0; k < n; ++k) std::swap(z[cm(k, j, n)], z[cm(k, j + 1, n)]);
}

// LAMBDA 降相关步骤，通过 Gauss 变换和置换减小模糊度相关性。
static void lambdaReduction(int n, std::vector<double> &l, std::vector<double> &d, std::vector<double> &z) {
    // LAMBDA step 2: integer decorrelation.
    // Gauss transforms and permutations reduce correlation between ambiguities,
    // making the later integer search more stable.
    int j = n - 2;
    int k = n - 2;
    while (j >= 0) {
        if (j <= k) {
            for (int i = j + 1; i < n; ++i) lambdaGauss(n, l, z, i, j);
        }
        const double del = d[j] + l[cm(j + 1, j, n)] * l[cm(j + 1, j, n)] * d[j + 1];
        if (del + 1e-6 < d[j + 1]) {
            lambdaPerm(n, l, d, j, del, z);
            k = j;
            j = n - 2;
        } else {
            --j;
        }
    }
}

static int lambdaSearch(int n, int m, const std::vector<double> &l, const std::vector<double> &d,
                        const std::vector<double> &zs, std::vector<double> &zn, std::vector<double> &s) {
    // LAMBDA step 3: search integer candidates.
    // The function keeps the best m candidates and their squared residuals.
    // The two best residuals are used for the ratio test.
    constexpr int kLoopMax = 10000;
    int k = n - 1;
    int nn = 0;
    int imax = 0;
    double maxdist = 1e99;
    std::vector<double> ss(n * n, 0.0), dist(n, 0.0), zb(n, 0.0), z(n, 0.0), step(n, 0.0);
    dist[k] = 0.0;
    zb[k] = zs[k];
    z[k] = lambdaRound(zb[k]);
    double y = zb[k] - z[k];
    step[k] = y <= 0.0 ? -1.0 : 1.0;
    int c = 0;
    for (; c < kLoopMax; ++c) {
        const double newdist = dist[k] + y * y / d[k];
        if (newdist < maxdist) {
            if (k != 0) {
                dist[--k] = newdist;
                for (int i = 0; i <= k; ++i) {
                    ss[cm(k, i, n)] = ss[cm(k + 1, i, n)] + (z[k + 1] - zb[k + 1]) * l[cm(k + 1, i, n)];
                }
                zb[k] = zs[k] + ss[cm(k, k, n)];
                z[k] = lambdaRound(zb[k]);
                y = zb[k] - z[k];
                step[k] = y <= 0.0 ? -1.0 : 1.0;
            } else {
                if (nn < m) {
                    if (nn == 0 || newdist > s[imax]) imax = nn;
                    for (int i = 0; i < n; ++i) zn[cm(i, nn, n)] = z[i];
                    s[nn++] = newdist;
                } else if (newdist < s[imax]) {
                    for (int i = 0; i < n; ++i) zn[cm(i, imax, n)] = z[i];
                    s[imax] = newdist;
                    imax = 0;
                    for (int i = 1; i < m; ++i) {
                        if (s[imax] < s[i]) imax = i;
                    }
                }
                if (nn >= m) maxdist = s[imax];
                z[0] += step[0];
                y = zb[0] - z[0];
                step[0] = -step[0] - (step[0] <= 0.0 ? -1.0 : 1.0);
            }
        } else {
            if (k == n - 1) break;
            ++k;
            z[k] += step[k];
            y = zb[k] - z[k];
            step[k] = -step[k] - (step[k] <= 0.0 ? -1.0 : 1.0);
        }
    }
    for (int i = 0; i < m - 1; ++i) {
        for (int j = i + 1; j < m; ++j) {
            if (s[i] <= s[j]) continue;
            std::swap(s[i], s[j]);
            for (int k2 = 0; k2 < n; ++k2) std::swap(zn[cm(k2, i, n)], zn[cm(k2, j, n)]);
        }
    }
    return (c >= kLoopMax || nn < m) ? -1 : 0;
}

static bool lambdaIntegerSearch(const std::vector<double> &a, Matrix q, std::vector<int> &best,
                                double &ratio, double *s0Out = nullptr, double *s1Out = nullptr) {
    // Public wrapper for ambiguity fixing.
    // Input:
    //   a: float ambiguity vector
    //   q: covariance matrix of the float ambiguities
    // Output:
    //   best: selected integer ambiguity vector
    //   ratio: second-best residual / best residual
    // If ratio is high enough, the fixed solution is considered reliable.
    const int n = static_cast<int>(a.size());
    if (n < 4 || n > 64) return false;
    for (int i = 0; i < n; ++i) {
        q(i, i) += 1e-10;
        for (int j = i + 1; j < n; ++j) {
            const double s = 0.5 * (q(i, j) + q(j, i));
            q(i, j) = q(j, i) = s;
        }
    }

    std::vector<double> l(n * n, 0.0), d(n, 0.0), z(n * n, 0.0);
    for (int i = 0; i < n; ++i) z[cm(i, i, n)] = 1.0;
    if (lambdaLd(n, q, l, d) != 0) return false;
    lambdaReduction(n, l, d, z);

    std::vector<double> zs(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) zs[i] += z[cm(j, i, n)] * a[j];
    }
    std::vector<double> zn(n * 2, 0.0), s(2, 1e99);
    if (lambdaSearch(n, 2, l, d, zs, zn, s) != 0 || s[0] <= 0.0 || s[1] >= 1e98) return false;

    Matrix zt(n, n);
    for (int r = 0; r < n; ++r) for (int c = 0; c < n; ++c) zt(r, c) = z[cm(c, r, n)];
    std::vector<double> e(n, 0.0);
    for (int i = 0; i < n; ++i) e[i] = zn[cm(i, 0, n)];
    const std::vector<double> f = solveLinear(zt, e);
    if (f.empty()) return false;
    best.assign(n, 0);
    for (int i = 0; i < n; ++i) best[i] = static_cast<int>(std::llround(f[i]));
    ratio = s[1] / s[0];
    if (s0Out) *s0Out = s[0];
    if (s1Out) *s1Out = s[1];
    return true;
}

struct FixDiagnostics {
    int n = 0;
    double qMin = 0.0;
    double qMax = 0.0;
    double qDiagCond = 0.0;
    double fracRms = 0.0;
    double s0 = 0.0;
    double s1 = 0.0;
    std::string status = "not_run";
};

// 检查浮点模糊度协方差质量，输出给调试日志解释固定失败原因。
static FixDiagnostics inspectAmbiguityCovariance(const std::vector<double> &afloat, const Matrix &qnn) {
    FixDiagnostics d;
    d.n = static_cast<int>(afloat.size());
    if (d.n <= 0) return d;
    d.qMin = 1e99;
    d.qMax = 0.0;
    double frac2 = 0.0;
    for (int i = 0; i < d.n; ++i) {
        const double qi = qnn(i, i);
        d.qMin = std::min(d.qMin, qi);
        d.qMax = std::max(d.qMax, qi);
        const double frac = afloat[i] - std::round(afloat[i]);
        frac2 += frac * frac;
    }
    d.qDiagCond = d.qMin > 0.0 ? d.qMax / d.qMin : 1e99;
    d.fracRms = std::sqrt(frac2 / std::max(1, d.n));
    return d;
}

// 对选中的双差模糊度尝试 LAMBDA 固定，并用固定整数更新基线坐标。
static bool tryFixAmbiguities(const RtkFilterState &flt, const std::vector<std::pair<int, int>> &activeDd,
                              const RtkOptions &opt, Vec3 &fixedPos, double &ratio, int &nfix,
                              std::vector<int> *fixedAmb = nullptr, std::vector<std::pair<int, int>> *fixedDd = nullptr,
                              FixDiagnostics *diag = nullptr) {
    std::vector<std::pair<int, int>> dd;
    for (const auto &p : activeDd) {
        if (p.first < kRtkAmbStart || p.second < kRtkAmbStart) continue;
        if (p.first >= static_cast<int>(flt.x.size()) || p.second >= static_cast<int>(flt.x.size())) continue;
        const int lockS = flt.ambLock.count(p.first) ? flt.ambLock.at(p.first) : 0;
        const int lockR = flt.ambLock.count(p.second) ? flt.ambLock.at(p.second) : 0;
        if (lockS < opt.minAmbLock || lockR < opt.minAmbLock) continue;
        if (std::find(dd.begin(), dd.end(), p) == dd.end()) dd.push_back(p);
    }
    auto ddVar = [&](const std::pair<int, int> &p) {
        return flt.p(p.first, p.first) + flt.p(p.second, p.second) - 2.0 * flt.p(p.first, p.second);
    };
    std::sort(dd.begin(), dd.end(), [&](const auto &a, const auto &b) { return ddVar(a) < ddVar(b); });
    const int maxN = opt.allowPartialAr ? std::min(8, static_cast<int>(dd.size())) : static_cast<int>(dd.size());
    if (maxN < 4) return false;

    bool ok = false;
    double bestRatio = ratio;
    int bestN = 0;
    Vec3 bestPos;
    std::vector<int> bestAmb;
    std::vector<std::pair<int, int>> bestDd;
    for (int k = maxN; k >= 4; --k) {
        std::vector<double> afloat(k, 0.0);
        Matrix qnn(k, k), qbn(3, k);
        for (int i = 0; i < k; ++i) {
            const int si = dd[i].first;
            const int ri = dd[i].second;
            afloat[i] = flt.x[si] - flt.x[ri];
            for (int j = 0; j < k; ++j) {
                const int sj = dd[j].first;
                const int rj = dd[j].second;
                qnn(i, j) = flt.p(si, sj) - flt.p(si, rj) - flt.p(ri, sj) + flt.p(ri, rj);
            }
            for (int r = 0; r < 3; ++r) qbn(r, i) = flt.p(r, si) - flt.p(r, ri);
        }

        FixDiagnostics cur = inspectAmbiguityCovariance(afloat, qnn);
        std::vector<int> nint;
        double rr = 0.0;
        double s0 = 0.0, s1 = 0.0;
        if (!lambdaIntegerSearch(afloat, qnn, nint, rr, &s0, &s1)) {
            cur.status = "lambda_fail";
            if (diag && diag->status == "not_run") *diag = cur;
            if (!opt.allowPartialAr) break;
            continue;
        }
        cur.s0 = s0;
        cur.s1 = s1;
        cur.status = rr >= opt.lambdaRatioThreshold ? "ratio_ok" : "ratio_low";
        if (diag && (!ok || cur.n > diag->n)) *diag = cur;
        if (rr < opt.lambdaRatioThreshold) {
            if (!opt.allowPartialAr) break;
            continue;
        }

        Matrix qnnInv;
        try {
            qnnInv = inverse(qnn);
        } catch (...) {
            continue;
        }
        std::vector<double> diff(k, 0.0);
        for (int i = 0; i < k; ++i) diff[i] = afloat[i] - static_cast<double>(nint[i]);
        const std::vector<double> gainDiff = (qbn * qnnInv) * diff;
        Vec3 candPos{flt.x[0] - gainDiff[0], flt.x[1] - gainDiff[1], flt.x[2] - gainDiff[2]};
        if (norm(candPos - Vec3{flt.x[0], flt.x[1], flt.x[2]}) >= 2.0) continue;
        if (!ok || k > bestN || (k == bestN && rr > bestRatio)) {
            ok = true;
            bestRatio = rr;
            bestN = k;
            bestPos = candPos;
            bestAmb = nint;
            bestDd.assign(dd.begin(), dd.begin() + k);
        }
        if (!opt.allowPartialAr) break;
    }
    if (!ok) return false;
    ratio = bestRatio;
    nfix = bestN;
    fixedPos = bestPos;
    if (fixedAmb) *fixedAmb = bestAmb;
    if (fixedDd) *fixedDd = bestDd;
    return true;
}

// 按 GPS/BDS 系统分别尝试固定，用于系统间质量差异较大时的部分固定。
static bool tryFixSystemAmbiguities(const RtkFilterState &flt, const std::vector<std::pair<int, int>> &activeDd,
                                    const std::vector<bool> &isBdsAmb, bool wantBds, const RtkOptions &opt,
                                    Vec3 &fixedPos, double &ratio, int &nfix) {
    std::vector<std::pair<int, int>> sub;
    for (const auto &p : activeDd) {
        if (p.first < 0 || p.second < 0 || p.first >= static_cast<int>(isBdsAmb.size()) ||
            p.second >= static_cast<int>(isBdsAmb.size())) {
            continue;
        }
        if (isBdsAmb[p.first] == wantBds && isBdsAmb[p.second] == wantBds) sub.push_back(p);
    }
    return tryFixAmbiguities(flt, sub, opt, fixedPos, ratio, nfix);
}

struct ActiveDdMeas {
    std::pair<int, int> idx;
    int satIndex = -1;
    int refIndex = -1;
    int freq = 0;
};

// 为每个系统/频点选择高度角最高的参考星，构造站星双差。
static std::map<std::pair<System, int>, int> selectDdRefs(const std::vector<DdSat> &sats) {
    std::map<std::pair<System, int>, int> ref;
    for (int i = 0; i < static_cast<int>(sats.size()); ++i) {
        for (int f = 0; f < 2; ++f) {
            if (!sats[i].ok[f]) continue;
            auto key = std::make_pair(sats[i].sat.sys, f);
            auto it = ref.find(key);
            if (it == ref.end() || sats[i].elev > sats[it->second].elev) ref[key] = i;
        }
    }
    return ref;
}

// 尽量保持参考星连续，只有参考星不可用或高度太低时才切换。
static std::map<std::pair<System, int>, int> selectStableDdRefs(const std::vector<DdSat> &sats,
                                                                RtkFilterState &flt) {
    std::map<std::pair<System, int>, int> refs;
    for (int f = 0; f < 2; ++f) {
        for (System sys : {System::GPS, System::BDS, System::GLO}) {
            const auto key = std::make_pair(sys, f);
            int keep = -1;
            auto pit = flt.refPrn.find(key);
            if (pit != flt.refPrn.end()) {
                for (int i = 0; i < static_cast<int>(sats.size()); ++i) {
                    if (sats[i].ok[f] && sats[i].sat.sys == sys && sats[i].sat.prn == pit->second &&
                        sats[i].elev >= 20.0) {
                        keep = i;
                        break;
                    }
                }
            }
            if (keep >= 0) {
                refs[key] = keep;
                continue;
            }
            int best = -1;
            for (int i = 0; i < static_cast<int>(sats.size()); ++i) {
                if (!sats[i].ok[f] || sats[i].sat.sys != sys) continue;
                if (best < 0 || sats[i].elev > sats[best].elev) best = i;
            }
            if (best >= 0) {
                refs[key] = best;
                flt.refPrn[key] = sats[best].sat.prn;
            }
        }
    }
    return refs;
}

// 检查流动站和基准站在某频点上是否同时具有码和相位观测。
static bool usableDualFreq(System sys, const ObsValue &r, const ObsValue &b, int f) {
    if (sys == System::BDS && f > 1) return false;
    if (sys == System::GPS && f > 1) return false;
    if (sys == System::GLO && f > 1) return false;
    return r.hasCode[f] && r.hasPhase[f] && b.hasCode[f] && b.hasPhase[f];
}

// 收集当前历元可用于双差的共视卫星，并计算单差码/相位观测。
static std::vector<DdSat> collectDdSats(const ObsEpoch &re, const ObsEpoch &be, const NavData &nav,
                                        const Vec3 &roverApprox, const RtkOptions &opt) {
    std::vector<DdSat> sats;
    const Llh approx = ecefToLlh(roverApprox);
    for (const auto &kv : re.obs) {
        const SatId sat = kv.first;
        if (sat.sys != System::GPS && sat.sys != System::BDS && sat.sys != System::GLO) continue;
        if (sat.sys == System::GPS && !opt.useGps) continue;
        if (sat.sys == System::BDS && !opt.useBds) continue;
        if (sat.sys == System::GLO && !opt.useGlo) continue;
        if (sat.sys == System::BDS && !opt.useBdsGeo && sat.prn <= 5) continue;
        auto bit = be.obs.find(sat);
        if (bit == be.obs.end()) continue;
        SatState ss;
        double rho;
        Vec3 los;
        const ObsValue &rov = kv.second;
        const ObsValue &bas = bit->second;
        const double elCode = rov.hasCode[0] ? rov.code[0] : (rov.hasCode[1] ? rov.code[1] : 0.0);
        if (!observedGeometry(nav, sat, re.time, elCode, roverApprox, rho, los)) continue;
        const Vec3 nedLos = ecefToNed(los, approx);
        const double el = std::asin(-nedLos.z) * kR2D;
        if (el < opt.elevationMaskDeg) continue;
        DdSat d;
        d.sat = sat;
        d.los = los;
        d.elev = el;
        for (int f = 0; f < 2; ++f) {
            if (!usableDualFreq(sat.sys, rov, bas, f)) continue;
            const double lam = carrierWavelength(sat, f);
            d.roverCode[f] = rov.code[f];
            d.baseCode[f] = bas.code[f];
            d.roverPhaseM[f] = rov.phase[f] * lam;
            d.basePhaseM[f] = bas.phase[f] * lam;
            d.sdCode[f] = rov.code[f] - bas.code[f];            d.sdPhaseM[f] = (rov.phase[f] - bas.phase[f]) * lam;
            d.roverLli[f] = rov.lli[f];
            d.baseLli[f] = bas.lli[f];
            d.ok[f] = true;
        }
        if (d.ok[0] || d.ok[1]) sats.push_back(d);
    }
    return sats;
}

// 计算某颗目标星相对参考星的双差几何距离和视线方向。
static bool ddGeometry(const ObsEpoch &re, const ObsEpoch &be, const NavData &nav, const Vec3 &rover,
                       const Vec3 &basePos, const DdSat &sat, const DdSat &ref, int f,
                       double &geomDd, Vec3 &losSat, Vec3 &losRef) {
    double rhoRu = 0.0, rhoRb = 0.0, rhoRefU = 0.0, rhoRefB = 0.0;
    Vec3 tmp;
    if (!observedGeometry(nav, sat.sat, re.time, sat.roverCode[f], rover, rhoRu, losSat)) return false;
    if (!observedGeometry(nav, sat.sat, be.time, sat.baseCode[f], basePos, rhoRb, tmp)) return false;
    if (!observedGeometry(nav, ref.sat, re.time, ref.roverCode[f], rover, rhoRefU, losRef)) return false;
    if (!observedGeometry(nav, ref.sat, be.time, ref.baseCode[f], basePos, rhoRefB, tmp)) return false;
    geomDd = (rhoRu - rhoRb) - (rhoRefU - rhoRefB);
    return true;
}

// 计算单颗卫星的站间单差几何距离。
static bool sdGeometry(const ObsEpoch &re, const ObsEpoch &be, const NavData &nav, const Vec3 &rover,
                       const Vec3 &basePos, const DdSat &sat, int f, double &geomSd, Vec3 &losSat) {
    double rhoRu = 0.0, rhoRb = 0.0;
    Vec3 tmp;
    if (!observedGeometry(nav, sat.sat, re.time, sat.roverCode[f], rover, rhoRu, losSat)) return false;
    if (!observedGeometry(nav, sat.sat, be.time, sat.baseCode[f], basePos, rhoRb, tmp)) return false;
    geomSd = rhoRu - rhoRb;
    return true;
}

// 固定解质量检查：用固定整数重新计算相位残差，RMS/最大残差超限则拒绝固定。
static bool validateFixedPhaseResiduals(const ObsEpoch &re, const ObsEpoch &be, const NavData &nav,
                                        const Vec3 &basePos, const std::vector<DdSat> &sats,
                                        const Vec3 &fixedPos,
                                        const std::vector<ActiveDdMeas> &activeMeas,
                                        const std::vector<std::pair<int, int>> &fixedDd,
                                        const std::vector<int> &fixedAmb, double rmsMax, double absMax,
                                        double &rms, double &maxAbs, int &used) {
    rms = 0.0;
    maxAbs = 0.0;
    used = 0;
    if (fixedDd.size() != fixedAmb.size() || fixedDd.empty()) return false;
    double sum2 = 0.0;
    for (size_t i = 0; i < fixedDd.size(); ++i) {
        const auto it = std::find_if(activeMeas.begin(), activeMeas.end(), [&](const ActiveDdMeas &m) {
            return m.idx == fixedDd[i];
        });
        if (it == activeMeas.end()) return false;
        if (it->satIndex < 0 || it->satIndex >= static_cast<int>(sats.size()) || it->refIndex < 0 ||
            it->refIndex >= static_cast<int>(sats.size())) {
            return false;
        }
        const DdSat &s = sats[it->satIndex];
        const DdSat &r = sats[it->refIndex];
        const int f = it->freq;
        double geomDd = 0.0;
        Vec3 losS, losR;
        if (!ddGeometry(re, be, nav, fixedPos, basePos, s, r, f, geomDd, losS, losR)) return false;
        const double obsDd = s.sdPhaseM[f] - r.sdPhaseM[f];
        const double lam = carrierWavelength(s.sat, f);
        const double res = obsDd - (geomDd + lam * static_cast<double>(fixedAmb[i]));
        sum2 += res * res;
        maxAbs = std::max(maxAbs, std::abs(res));
        ++used;
    }
    if (used <= 0) return false;
    rms = std::sqrt(sum2 / used);
    return rms <= rmsMax && maxAbs <= absMax;
}

// 仅用双差伪距求一个粗略基线坐标，作为 RTK 滤波初值或弱约束。
static bool solveDdCodePosition(const ObsEpoch &re, const ObsEpoch &be, const NavData &nav, const Vec3 &basePos,
                                const std::vector<DdSat> &sats, const RtkOptions &opt, Vec3 &rover,
                                double &rms, int &rowsUsed) {
    const auto refs = selectDdRefs(sats);
    rowsUsed = 0;
    rms = 0.0;
    for (int iter = 0; iter < 8; ++iter) {
        std::vector<std::vector<double>> rows;
        std::vector<double> rhs;
        std::vector<double> weights;
        for (int i = 0; i < static_cast<int>(sats.size()); ++i) {
            for (int f = 0; f < 2; ++f) {
                if (!sats[i].ok[f]) continue;
                auto rit = refs.find(std::make_pair(sats[i].sat.sys, f));
                if (rit == refs.end() || rit->second == i) continue;
                const DdSat &s = sats[i];
                const DdSat &r = sats[rit->second];
                double geomDd = 0.0;
                Vec3 losS, losR;
                if (!ddGeometry(re, be, nav, rover, basePos, s, r, f, geomDd, losS, losR)) continue;
                const double obsDd = s.sdCode[f] - r.sdCode[f];
                const Vec3 a = losR - losS;
                const double sinS = std::max(std::sin(std::max(s.elev, 5.0) * kD2R), 0.1);
                const double sinR = std::max(std::sin(std::max(r.elev, 5.0) * kD2R), 0.1);
                const double ddWeight = 1.0 / sqr(sinS) + 1.0 / sqr(sinR);
                const double w = 1.0 / std::max(opt.codeSigma * std::sqrt(ddWeight), 1e-6);
                rows.push_back({a.x * w, a.y * w, a.z * w});
                rhs.push_back((obsDd - geomDd) * w);
                weights.push_back(w);
            }
        }
        if (rows.size() < 4) return false;
        Matrix h(static_cast<int>(rows.size()), 3);
        for (int r = 0; r < h.rows(); ++r) for (int c = 0; c < 3; ++c) h(r, c) = rows[r][c];
        std::vector<double> dx;
        if (!leastSquares(h, rhs, dx, nullptr)) return false;
        rover.x += dx[0]; rover.y += dx[1]; rover.z += dx[2];
        rowsUsed = static_cast<int>(rows.size());
        double sum = 0.0;
        for (int r = 0; r < rowsUsed; ++r) {
            const double fit = rows[r][0] * dx[0] + rows[r][1] * dx[1] + rows[r][2] * dx[2];
            const double unweighted = (rhs[r] - fit) / std::max(weights[r], 1e-12);
            sum += unweighted * unweighted;
        }
        rms = std::sqrt(sum / std::max(1, rowsUsed - 3));
        if (std::sqrt(sqr(dx[0]) + sqr(dx[1]) + sqr(dx[2])) < 1e-4) return true;
    }
    return true;
}

struct DdAmbInfo {
    int stateIdx = -1;
    int refStateIdx = -1;
    int satIndex = -1;
    int refIndex = -1;
    int freq = 0;
    double elev = 0.0;
};

struct DdUpdateObs {
    Vec3 hpos;
    double innov = 0.0;
    double variance = 0.0;
    double refVariance = 0.0;
    int ambIdx = -1;
    int refAmbIdx = -1;
    int satIndex = -1;
    int refIndex = -1;
    int freq = 0;
    bool phase = false;
    System sys = System::OTHER;
    double elev = 0.0;
};

// 更新并检测 LLI 周跳标志。
static bool updateLliSlip(std::map<std::string, int> &last, const std::string &key, int lli) {
    return lliIndicatesSlip(last, key, lli);
}

// 更新 GF/MW 组合历史值，并检测组合值是否突跳。
static bool updateCombinationJump(std::map<std::string, double> &last, const std::string &key,
                                  double value, double threshold) {
    if (!std::isfinite(value)) return false;
    auto it = last.find(key);
    const bool known = it != last.end();
    const double prev = known ? it->second : value;
    last[key] = value;
    return known && std::abs(value - prev) > threshold;
}

// 使用 GF 和 MW 组合联合判断单接收机端是否发生周跳。
static bool detectReceiverComboSlip(std::map<std::string, double> &lastGf,
                                    std::map<std::string, double> &lastMw,
                                    const DdSat &s, bool rover) {
    if (!s.ok[0] || !s.ok[1]) return false;
    const double l1 = carrierWavelength(s.sat, 0);
    const double l2 = carrierWavelength(s.sat, 1);
    const double f1 = kClight / l1;
    const double f2 = kClight / l2;
    const double p1 = rover ? s.roverCode[0] : s.baseCode[0];
    const double p2 = rover ? s.roverCode[1] : s.baseCode[1];
    const double l1m = rover ? s.roverPhaseM[0] : s.basePhaseM[0];
    const double l2m = rover ? s.roverPhaseM[1] : s.basePhaseM[1];
    if (p1 <= 0.0 || p2 <= 0.0) return false;

    const std::string key = satToString(s.sat);
    const double gf = l1m - l2m;
    const double mw = (l1m - l2m) - (f1 * p1 + f2 * p2) / (f1 + f2);
    const bool gfSlip = updateCombinationJump(lastGf, key, gf, 0.50);
    const bool mwSlip = updateCombinationJump(lastMw, key, mw, 20.0);
    return gfSlip && mwSlip;
}

// 在模糊度固定后，仅用固定相位双差重新精化基线位置。
static bool refineFixedDdPosition(const ObsEpoch &re, const ObsEpoch &be, const NavData &nav, const Vec3 &basePos,
                                  const std::vector<DdSat> &sats, const std::vector<DdAmbInfo> &ambs,
                                  const std::vector<int> &nfix, Vec3 &pos,
                                  double &rms, double &maxAbs, int &used) {
    if (ambs.size() != nfix.size() || ambs.size() < 3) return false;
    for (int iter = 0; iter < 4; ++iter) {
        std::vector<std::vector<double>> rows;
        std::vector<double> rhs;
        for (size_t i = 0; i < ambs.size(); ++i) {
            const DdAmbInfo &aInfo = ambs[i];
            if (aInfo.satIndex < 0 || aInfo.refIndex < 0) continue;
            const DdSat &s = sats[aInfo.satIndex];
            const DdSat &r = sats[aInfo.refIndex];
            double geomDd = 0.0;
            Vec3 losS, losR;
            if (!ddGeometry(re, be, nav, pos, basePos, s, r, aInfo.freq, geomDd, losS, losR)) continue;
            const double lam = carrierWavelength(s.sat, aInfo.freq);
            const double obsDd = s.sdPhaseM[aInfo.freq] - r.sdPhaseM[aInfo.freq];
            const Vec3 h = losR - losS;
            rows.push_back({h.x, h.y, h.z});
            rhs.push_back(obsDd - lam * static_cast<double>(nfix[i]) - geomDd);
        }
        if (rows.size() < 3) return false;
        Matrix h(static_cast<int>(rows.size()), 3);
        for (int r = 0; r < h.rows(); ++r) for (int c = 0; c < 3; ++c) h(r, c) = rows[r][c];
        std::vector<double> dx;
        if (!leastSquares(h, rhs, dx, nullptr)) return false;
        pos.x += dx[0]; pos.y += dx[1]; pos.z += dx[2];
        if (std::sqrt(sqr(dx[0]) + sqr(dx[1]) + sqr(dx[2])) < 1e-4) break;
    }

    double sum2 = 0.0;
    rms = 0.0;
    maxAbs = 0.0;
    used = 0;
    for (size_t i = 0; i < ambs.size(); ++i) {
        const DdAmbInfo &aInfo = ambs[i];
        const DdSat &s = sats[aInfo.satIndex];
        const DdSat &r = sats[aInfo.refIndex];
        double geomDd = 0.0;
        Vec3 losS, losR;
        if (!ddGeometry(re, be, nav, pos, basePos, s, r, aInfo.freq, geomDd, losS, losR)) continue;
        const double lam = carrierWavelength(s.sat, aInfo.freq);
        const double obsDd = s.sdPhaseM[aInfo.freq] - r.sdPhaseM[aInfo.freq];
        const double res = obsDd - geomDd - lam * static_cast<double>(nfix[i]);
        sum2 += res * res;
        maxAbs = std::max(maxAbs, std::abs(res));
        ++used;
    }
    if (used <= 0) return false;
    rms = std::sqrt(sum2 / used);
    return true;
}

// 计算一个双差模糊度组合的方差，用于排序选择质量更好的模糊度。
static double ddAmbVariance(const RtkFilterState &flt, const DdAmbInfo &a) {
    if (a.stateIdx < 0 || a.refStateIdx < 0 ||
        a.stateIdx >= static_cast<int>(flt.x.size()) || a.refStateIdx >= static_cast<int>(flt.x.size())) {
        return 1e99;
    }
    return flt.p(a.stateIdx, a.stateIdx) + flt.p(a.refStateIdx, a.refStateIdx) -
           2.0 * flt.p(a.stateIdx, a.refStateIdx);
}

// 联合双差模糊度固定：选择锁定充分的模糊度，LAMBDA 固定并通过残差质量检验。
static bool tryJointDdFix(const ObsEpoch &re, const ObsEpoch &be, const NavData &nav,
                          const RtkFilterState &flt, const RtkOptions &opt,
                          const std::vector<DdSat> &sats, const std::vector<DdAmbInfo> &activeAmbs,
                          Vec3 &fixedPos, double &ratio, int &nfixTotal,
                          double &fixRms, double &fixMax, int &fixRows,
                          std::vector<DdAmbInfo> &fixedAmbs, std::vector<int> &fixedInts,
                          std::ostream *debug) {
    std::vector<DdAmbInfo> ambs;
    for (const auto &a : activeAmbs) {
        if (a.stateIdx < 0 || a.refStateIdx < 0 || a.satIndex < 0 || a.refIndex < 0) continue;
        if (a.stateIdx >= static_cast<int>(flt.x.size()) || a.refStateIdx >= static_cast<int>(flt.x.size())) continue;
        if (sats[a.satIndex].sat.sys == System::GLO) continue;
        const int satLock = flt.ambLock.count(a.stateIdx) ? flt.ambLock.at(a.stateIdx) : 0;
        const int refLock = flt.ambLock.count(a.refStateIdx) ? flt.ambLock.at(a.refStateIdx) : 0;
        if (satLock < opt.minAmbLock || refLock < opt.minAmbLock) continue;
        if (ddAmbVariance(flt, a) >= 1e90) continue;
        ambs.push_back(a);
    }
    std::sort(ambs.begin(), ambs.end(), [&](const DdAmbInfo &a, const DdAmbInfo &b) {
        const double qa = ddAmbVariance(flt, a);
        const double qb = ddAmbVariance(flt, b);
        if (std::abs(qa - qb) > 1e-9) return qa < qb;
        return a.elev > b.elev;
    });

    const int maxN = opt.allowPartialAr ? std::min(24, static_cast<int>(ambs.size())) : static_cast<int>(ambs.size());
    const int minN = 8;
    if (maxN < minN) {
        if (debug) *debug << " ar_joint_skip n=" << ambs.size();
        return false;
    }

    bool ok = false;
    double bestRatio = 0.0;
    int bestN = 0;
    Vec3 bestPos{flt.x[0], flt.x[1], flt.x[2]};
    double bestRms = 0.0, bestMax = 0.0;
    int bestRows = 0;
    std::vector<DdAmbInfo> bestAmbs;
    std::vector<int> bestInts;

    for (int ndd = maxN; ndd >= minN; --ndd) {
        std::vector<DdAmbInfo> sub(ambs.begin(), ambs.begin() + ndd);
        std::vector<double> afloat(ndd, 0.0);
        Matrix qnn(ndd, ndd), qbn(3, ndd);
        bool valid = true;
        for (int i = 0; i < ndd; ++i) {
            const int si = sub[i].stateIdx;
            const int ri = sub[i].refStateIdx;
            afloat[i] = flt.x[si] - flt.x[ri];
            for (int j = 0; j < ndd; ++j) {
                const int sj = sub[j].stateIdx;
                const int rj = sub[j].refStateIdx;
                qnn(i, j) = flt.p(si, sj) - flt.p(si, rj) - flt.p(ri, sj) + flt.p(ri, rj);
            }
            if (qnn(i, i) <= 0.0 || !std::isfinite(qnn(i, i))) valid = false;
            for (int r = 0; r < 3; ++r) qbn(r, i) = flt.p(r, si) - flt.p(r, ri);
        }
        if (!valid) continue;

        std::vector<int> nint;
        double rr = 0.0;
        double s0 = 0.0, s1 = 0.0;
        if (!lambdaIntegerSearch(afloat, qnn, nint, rr, &s0, &s1)) {
            if (debug) *debug << " ar_joint ndd=" << ndd << "/lambda_fail";
            continue;
        }
        if (debug) *debug << " ar_joint ndd=" << ndd << " ratio=" << rr;
        if (rr < opt.lambdaRatioThreshold) {
            if (debug) *debug << "/ratio_low";
            continue;
        }

        Matrix qnnInv;
        try {
            qnnInv = inverse(qnn);
        } catch (...) {
            if (debug) *debug << "/qinv_fail";
            continue;
        }
        std::vector<double> diff(ndd, 0.0);
        for (int i = 0; i < ndd; ++i) diff[i] = afloat[i] - static_cast<double>(nint[i]);
        const std::vector<double> posCorr = (qbn * qnnInv) * diff;
        Vec3 cand{flt.x[0] - posCorr[0], flt.x[1] - posCorr[1], flt.x[2] - posCorr[2]};
        double frms = 0.0, fmax = 0.0;
        int frows = 0;
        if (!refineFixedDdPosition(re, be, nav, flt.basePosition, sats, sub, nint, cand, frms, fmax, frows)) {
            if (debug) *debug << "/refine_fail";
            continue;
        }
        if (frms > opt.fixedPhaseRmsMax || fmax > opt.fixedPhaseAbsMax) {
            if (debug) *debug << "/residual_bad(" << frms << "," << fmax << ")";
            continue;
        }
        const double step = norm(cand - Vec3{flt.x[0], flt.x[1], flt.x[2]});
        if (step > 30.0) {
            if (debug) *debug << "/step_bad(" << step << ")";
            continue;
        }

        if (!ok || ndd > bestN || (ndd == bestN && rr > bestRatio)) {
            ok = true;
            bestRatio = rr;
            bestN = ndd;
            bestPos = cand;
            bestRms = frms;
            bestMax = fmax;
            bestRows = frows;
            bestAmbs = sub;
            bestInts = nint;
        }
        if (!opt.allowPartialAr) break;
    }

    if (!ok) return false;
    fixedPos = bestPos;
    ratio = bestRatio;
    nfixTotal = bestN;
    fixRms = bestRms;
    fixMax = bestMax;
    fixRows = bestRows;
    fixedAmbs = bestAmbs;
    fixedInts = bestInts;
    return true;
}

// 固定保持策略：若上一小段固定质量可靠，则用接近整数的双差模糊度继续保持固定。
static bool tryHeldDdFix(const ObsEpoch &re, const ObsEpoch &be, const NavData &nav,
                         const RtkFilterState &flt, const RtkOptions &opt,
                         const std::vector<DdSat> &sats, const std::vector<DdAmbInfo> &activeAmbs,
                         const Vec3 *ddCodePos, double ddCodeRms,
                         Vec3 &heldPos, int &nfixTotal,
                         double &fixRms, double &fixMax, int &fixRows,
                         std::vector<DdAmbInfo> &heldAmbs, std::vector<int> &heldInts) {
    if (flt.fixedHoldAge > 3) return false;
    std::vector<DdAmbInfo> ambs;
    std::vector<int> ints;
    for (const auto &a : activeAmbs) {
        if (a.stateIdx < 0 || a.refStateIdx < 0 || a.satIndex < 0 || a.refIndex < 0) continue;
        if (a.stateIdx >= static_cast<int>(flt.x.size()) || a.refStateIdx >= static_cast<int>(flt.x.size())) continue;
        if (sats[a.satIndex].sat.sys == System::GLO) continue;
        const int satLock = flt.ambLock.count(a.stateIdx) ? flt.ambLock.at(a.stateIdx) : 0;
        const int refLock = flt.ambLock.count(a.refStateIdx) ? flt.ambLock.at(a.refStateIdx) : 0;
        if (satLock < opt.minAmbLock || refLock < opt.minAmbLock) continue;
        const double afloat = flt.x[a.stateIdx] - flt.x[a.refStateIdx];
        const int nint = static_cast<int>(std::llround(afloat));
        if (std::abs(afloat - static_cast<double>(nint)) > 0.04) continue;
        ambs.push_back(a);
        ints.push_back(nint);
    }
    std::sort(ambs.begin(), ambs.end(), [&](const DdAmbInfo &a, const DdAmbInfo &b) {
        const double qa = ddAmbVariance(flt, a);
        const double qb = ddAmbVariance(flt, b);
        if (std::abs(qa - qb) > 1e-9) return qa < qb;
        return a.elev > b.elev;
    });
    if (ambs.size() < 10) return false;
    if (ambs.size() > 18) ambs.resize(18);
    ints.resize(ambs.size());
    for (size_t i = 0; i < ambs.size(); ++i) {
        const double afloat = flt.x[ambs[i].stateIdx] - flt.x[ambs[i].refStateIdx];
        ints[i] = static_cast<int>(std::llround(afloat));
    }

    Vec3 cand{flt.x[0], flt.x[1], flt.x[2]};
    double frms = 0.0, fmax = 0.0;
    int frows = 0;
    if (!refineFixedDdPosition(re, be, nav, flt.basePosition, sats, ambs, ints, cand, frms, fmax, frows)) return false;
    if (frms > std::min(opt.fixedPhaseRmsMax, 0.04) || fmax > std::min(opt.fixedPhaseAbsMax, 0.12)) return false;
    if (norm(cand - Vec3{flt.x[0], flt.x[1], flt.x[2]}) > 5.0) return false;
    if (ddCodePos && ddCodeRms > 0.0 && norm(cand - *ddCodePos) > std::max(2.0, 4.0 * ddCodeRms)) return false;

    heldPos = cand;
    nfixTotal = static_cast<int>(ambs.size());
    fixRms = frms;
    fixMax = fmax;
    fixRows = frows;
    heldAmbs = ambs;
    heldInts = ints;
    return true;
}

// 单历元双差 RTK 解算：完成共视卫星筛选、周跳处理、Kalman 更新和模糊度固定。
static bool solveDdEpoch(const ObsEpoch &re, const ObsEpoch &be, const NavData &nav,
                         RtkFilterState &flt, const RtkOptions &opt, GnssSolution &sol,
                         const Vec3 &sppPos, std::ostream *debug = nullptr) {
    // Solve one rover/base common epoch using double-differenced RTK.
    //
    // Main steps:
    //   1. collectDdSats: collect satellites observed by both rover and base;
    //   2. selectStableDdRefs: choose reference satellites for each system/frequency;
    //   3. slip detection: use LLI and geometry-free/MW-like combinations to reset
    //      ambiguity states when a cycle slip is suspected;
    //   4. build code/phase double-difference residual rows;
    //   5. Kalman update float baseline and ambiguity states;
    //   6. attempt integer ambiguity fixing with LAMBDA and quality checks.
    const Vec3 roverApprox = flt.initialized ? Vec3{flt.x[0], flt.x[1], flt.x[2]} : sppPos;
    std::vector<DdSat> sats = collectDdSats(re, be, nav, roverApprox, opt);
    int gpsF[2] = {}, bdsF[2] = {}, gloF[2] = {};
    for (const auto &s : sats) {
        for (int f = 0; f < 2; ++f) {
            if (!s.ok[f]) continue;
            if (s.sat.sys == System::GPS) ++gpsF[f];
            if (s.sat.sys == System::BDS) ++bdsF[f];
            if (s.sat.sys == System::GLO) ++gloF[f];
        }
    }
    int usableCount = 0;
    for (const auto &s : sats) if (s.ok[0] || s.ok[1]) ++usableCount;

    if (debug) *debug << re.time.tow << " base_tow=" << be.time.tow << " dt=" << timediff(be.time, re.time)
                      << " common_usable=" << usableCount
                      << " gps=" << gpsF[0] << "/" << gpsF[1]
                      << " bds=" << bdsF[0] << "/" << bdsF[1]
                      << " glo=" << gloF[0] << "/" << gloF[1];
    if (usableCount < 4) {
        if (debug) *debug << " fail=not_enough_common\n";
        return false;
    }

    RtkFilterState prior = flt;
    if (!flt.initialized) {
        Vec3 initPos = sppPos;
        double crms = 0.0;
        int crows = 0;
        if (solveDdCodePosition(re, be, nav, flt.basePosition, sats, opt, initPos, crms, crows)) {
            if (debug) *debug << " init_ddcode_rms=" << crms << " rows=" << crows;
        }
        initRtkFilter(flt, re.time, initPos, flt.basePosition);
        if (debug) *debug << " init_pos=" << initPos.x << "," << initPos.y << "," << initPos.z;
    }
    predictRtkFilter(flt, re.time, opt);
    alignSdAmbiguityOffsets(flt, sats);

    const auto refs = selectStableDdRefs(sats, flt);
    int nCodeUpdate = 0, nPhaseUpdate = 0, nSlipReset = 0;
    int nComboSlip = 0, nLliSlip = 0, nPhaseJumpSlip = 0, nRefSlipReset = 0;
    double codeSum2 = 0.0, phaseSum2 = 0.0;
    int dopUsed = 0;
    double dopRms = 0.0;
    bool haveDopVel = false;
    Vec3 dopVel{};
    std::vector<DdAmbInfo> activeAmbs;
    std::vector<DdUpdateObs> obsRows;
    std::map<std::string, bool> slipBySatFreq;

    for (const auto &s : sats) {
        const bool comboSlip = opt.useGfMwSlip &&
                               (detectReceiverComboSlip(flt.lastRoverGf, flt.lastRoverMw, s, true) ||
                                detectReceiverComboSlip(flt.lastBaseGf, flt.lastBaseMw, s, false));
        for (int f = 0; f < 2; ++f) {
            if (!s.ok[f]) continue;
            const std::string key = ambStateKey(s.sat, f);
            const bool lliSlip = updateLliSlip(flt.lastRoverLli, key, s.roverLli[f]) ||
                                 updateLliSlip(flt.lastBaseLli, key, s.baseLli[f]);
            const bool slip = comboSlip || lliSlip;
            if (comboSlip) ++nComboSlip;
            if (lliSlip) ++nLliSlip;
            slipBySatFreq[key] = slip;
        }
    }

    const Vec3 linearizePos{flt.x[0], flt.x[1], flt.x[2]};
    for (const auto &rkv : refs) {
        const int refIdx = rkv.second;
        const int f = rkv.first.second;
        if (refIdx < 0 || refIdx >= static_cast<int>(sats.size()) || !sats[refIdx].ok[f]) continue;
        const DdSat &ref = sats[refIdx];
        const std::string refSlipKey = ambStateKey(ref.sat, f);
        double geomSdRef = 0.0;
        Vec3 losRefSd;
        if (!sdGeometry(re, be, nav, linearizePos, flt.basePosition, ref, f, geomSdRef, losRefSd)) continue;
        const double lamRef = carrierWavelength(ref.sat, f);
        const double initialRefAmb = (ref.sdPhaseM[f] - ref.sdCode[f]) / lamRef;
        const int refAmbIdx = ensureAmbState(flt, refSlipKey, initialRefAmb);
        if (refAmbIdx < 0) continue;
        auto refSlipIt = slipBySatFreq.find(refSlipKey);
        if (refSlipIt != slipBySatFreq.end() && refSlipIt->second) {
            resetAmbState(flt, refAmbIdx, initialRefAmb, 100.0);
            ++nSlipReset;
            ++nRefSlipReset;
        }
        for (int i = 0; i < static_cast<int>(sats.size()); ++i) {
            if (i == refIdx || !sats[i].ok[f]) continue;
            const DdSat &s = sats[i];
            if (s.sat.sys != ref.sat.sys) continue;
            double geomDd = 0.0;
            Vec3 losS, losR;
            if (!ddGeometry(re, be, nav, linearizePos, flt.basePosition, s, ref, f, geomDd, losS, losR)) continue;
            const Vec3 hpos = losR - losS;
            const double sinS = std::max(std::sin(std::max(s.elev, 5.0) * kD2R), 0.1);
            const double sinR = std::max(std::sin(std::max(ref.elev, 5.0) * kD2R), 0.1);
            const double ddWeight = (1.0 / sqr(sinS) + 1.0 / sqr(sinR));

            DdUpdateObs codeObs;
            codeObs.hpos = hpos;
            codeObs.innov = (s.sdCode[f] - ref.sdCode[f]) - geomDd;
            codeObs.variance = sqr(opt.codeSigma) * ddWeight;
            codeObs.refVariance = sqr(opt.codeSigma) / sqr(sinR);
            codeObs.satIndex = i;
            codeObs.refIndex = refIdx;
            codeObs.freq = f;
            codeObs.phase = false;
            codeObs.sys = s.sat.sys;
            codeObs.elev = std::min(s.elev, ref.elev);
            obsRows.push_back(codeObs);
            ++nCodeUpdate;
            codeSum2 += codeObs.innov * codeObs.innov;

            double geomSdSat = 0.0;
            Vec3 losSatSd;
            if (!sdGeometry(re, be, nav, linearizePos, flt.basePosition, s, f, geomSdSat, losSatSd)) continue;
            const std::string ambKey = ambStateKey(s.sat, f);
            const double lam = carrierWavelength(s.sat, f);
            const double obsDdPhase = s.sdPhaseM[f] - ref.sdPhaseM[f];
            const double initialAmb = (s.sdPhaseM[f] - s.sdCode[f]) / lam;
            const int ambIdx = ensureAmbState(flt, ambKey, initialAmb);
            if (ambIdx < 0) continue;

            int ambLock = 0;
            auto lockIt = flt.ambLock.find(ambIdx);
            if (lockIt != flt.ambLock.end()) ambLock = lockIt->second;

            const std::string satSlipKey = ambStateKey(s.sat, f);
            bool slip = false;
            auto sit = slipBySatFreq.find(satSlipKey);
            auto rit = slipBySatFreq.find(refSlipKey);
            slip = (sit != slipBySatFreq.end() && sit->second) || (rit != slipBySatFreq.end() && rit->second);
            double phaseInnov = obsDdPhase - (geomDd + lam * (flt.x[ambIdx] - flt.x[refAmbIdx]));
            if (!slip && ambLock >= opt.minAmbLock && std::abs(phaseInnov) > 5.0) {
                slip = true;
                ++nPhaseJumpSlip;
            }
            if (slip) {
                resetAmbState(flt, ambIdx, initialAmb, 100.0);
                ++nSlipReset;
                phaseInnov = obsDdPhase - (geomDd + lam * (flt.x[ambIdx] - flt.x[refAmbIdx]));
            }

            DdUpdateObs phaseObs;
            phaseObs.hpos = hpos;
            phaseObs.innov = phaseInnov;
            phaseObs.variance = sqr(opt.phaseSigma) * ddWeight;
            phaseObs.refVariance = sqr(opt.phaseSigma) / sqr(sinR);
            phaseObs.ambIdx = ambIdx;
            phaseObs.refAmbIdx = refAmbIdx;
            phaseObs.satIndex = i;
            phaseObs.refIndex = refIdx;
            phaseObs.freq = f;
            phaseObs.phase = true;
            phaseObs.sys = s.sat.sys;
            phaseObs.elev = std::min(s.elev, ref.elev);
            obsRows.push_back(phaseObs);
            ++nPhaseUpdate;
            phaseSum2 += phaseInnov * phaseInnov;
        }
    }

    if (!obsRows.empty()) {
        const int nState = static_cast<int>(flt.x.size());
        const int nObs = static_cast<int>(obsRows.size());
        Matrix h(nObs, nState);
        Matrix r(nObs, nObs);
        std::vector<double> v(nObs, 0.0);
        for (int i = 0; i < nObs; ++i) {
            const DdUpdateObs &row = obsRows[i];
            h(i, 0) = row.hpos.x;
            h(i, 1) = row.hpos.y;
            h(i, 2) = row.hpos.z;
            if (row.phase && row.ambIdx >= 0 && row.ambIdx < nState) {
                h(i, row.ambIdx) = carrierWavelength(sats[row.satIndex].sat, row.freq);
                if (row.refAmbIdx >= 0 && row.refAmbIdx < nState) {
                    h(i, row.refAmbIdx) = -carrierWavelength(sats[row.satIndex].sat, row.freq);
                }
            }
            v[i] = row.innov;
            r(i, i) = row.variance;
        }
        for (int i = 0; i < nObs; ++i) {
            for (int j = i + 1; j < nObs; ++j) {
                const DdUpdateObs &a = obsRows[i];
                const DdUpdateObs &b = obsRows[j];
                if (a.phase != b.phase || a.sys != b.sys || a.freq != b.freq || a.refIndex != b.refIndex) continue;
                const double cov = std::min(a.refVariance, b.refVariance);
                r(i, j) = r(j, i) = cov;
            }
        }
        if (!batchKalmanUpdate(flt, h, v, r)) {
            if (debug) *debug << " fail=batch_update\n";
            flt = prior;
            return false;
        }
        for (const auto &row : obsRows) {
            if (!row.phase || row.ambIdx < 0) continue;
            flt.ambLock[row.ambIdx] += 1;
            if (row.refAmbIdx >= 0) flt.ambLock[row.refAmbIdx] += 1;
            if (flt.ambLock[row.ambIdx] >= opt.minAmbLock) {
                activeAmbs.push_back({row.ambIdx, row.refAmbIdx, row.satIndex, row.refIndex, row.freq, row.elev});
            }
        }
        haveDopVel = dopplerVelocitySolve(re, nav, Vec3{flt.x[0], flt.x[1], flt.x[2]}, dopVel, dopUsed, dopRms);
        if (haveDopVel) {
            applyVelocityPseudoMeasurement(flt, dopVel, std::max(0.25, dopRms));
        }
    }

    if (debug) *debug << " codeUpd=" << nCodeUpdate << " phaseUpd=" << nPhaseUpdate
                      << " slipResets=" << nSlipReset << "(" << nComboSlip << ","
                      << nLliSlip << "," << nPhaseJumpSlip << "," << nRefSlipReset << ")";
    if (nCodeUpdate < 4 || nPhaseUpdate < 3) {
        if (debug) *debug << " fail=not_enough_updates\n";
        flt = prior;
        return false;
    }

    bool fixed = false;
    double ratio = 0.0;
    int nfixTotal = 0;
    Vec3 fixedPos{flt.x[0], flt.x[1], flt.x[2]};
    double bestFixRms = 0.0, bestFixMax = 0.0;
    int bestFixRows = 0;
    std::vector<DdAmbInfo> fixedAmbs;
    std::vector<int> fixedInts;
    if (tryJointDdFix(re, be, nav, flt, opt, sats, activeAmbs,
                      fixedPos, ratio, nfixTotal, bestFixRms, bestFixMax, bestFixRows,
                      fixedAmbs, fixedInts, debug)) {
        fixed = true;
        for (size_t i = 0; i < fixedAmbs.size(); ++i) {
            applyDdAmbiguityHold(flt, fixedAmbs[i].stateIdx, fixedAmbs[i].refStateIdx,
                                 static_cast<double>(fixedInts[i]), 0.01);
        }
    }

    Vec3 fltPos{flt.x[0], flt.x[1], flt.x[2]};
    Vec3 ddCodePos = fltPos;
    double ddCodeRms = 0.0;
    int ddCodeRows = 0;
    const bool haveDdCode = solveDdCodePosition(re, be, nav, flt.basePosition, sats, opt, ddCodePos, ddCodeRms, ddCodeRows);
    const bool goodDdCode = haveDdCode && ddCodeRows >= 8 && ddCodeRms < 0.6;
    const bool reliableDdCode = haveDdCode && ddCodeRows >= 12 && ddCodeRms < 1.5;
    bool heldFixed = false;
    if (!fixed && reliableDdCode &&
        tryHeldDdFix(re, be, nav, flt, opt, sats, activeAmbs, &ddCodePos, ddCodeRms,
                     fixedPos, nfixTotal, bestFixRms, bestFixMax, bestFixRows,
                     fixedAmbs, fixedInts)) {
        fixed = true;
        heldFixed = true;
        ratio = std::max(ratio, opt.lambdaRatioThreshold);
        for (size_t i = 0; i < fixedAmbs.size(); ++i) {
            applyDdAmbiguityHold(flt, fixedAmbs[i].stateIdx, fixedAmbs[i].refStateIdx,
                                 static_cast<double>(fixedInts[i]), 0.01);
        }
    }
    if (fixed && goodDdCode && norm(fixedPos - ddCodePos) > std::max(15.0, 8.0 * ddCodeRms)) {
        fixed = false;
        heldFixed = false;
        ratio = 0.0;
        nfixTotal = 0;
    }
    const double codeAnchorDist = haveDdCode ? norm(fltPos - ddCodePos) : 0.0;
    const bool anchorToCode = !fixed && reliableDdCode && codeAnchorDist > std::max(1.0, 2.0 * ddCodeRms);
    if (anchorToCode) {
        resetRtkAmbiguities(flt);
        flt.x[0] = ddCodePos.x;
        flt.x[1] = ddCodePos.y;
        flt.x[2] = ddCodePos.z;
        for (int i = 3; i < kClkGpsIdx; ++i) flt.x[i] = 0.0;
        for (int i = 0; i < 3; ++i) flt.p(i, i) = std::max(flt.p(i, i), sqr(std::max(1.0, ddCodeRms)));
        fltPos = ddCodePos;
    }
    sol.time = re.time;
    const bool useDdCodeFloat = reliableDdCode && !haveDopVel;
    sol.ecef = fixed ? fixedPos : (useDdCodeFloat ? ddCodePos : fltPos);
    if (!std::isfinite(sol.ecef.x) || !std::isfinite(sol.ecef.y) || !std::isfinite(sol.ecef.z) ||
        norm(sol.ecef - sppPos) > 500.0) {
        if (debug) *debug << " fail=rtk_diverged\n";
        flt = prior;
        return false;
    }
    sol.baseline = sol.ecef - flt.basePosition;
    sol.llh = ecefToLlh(sol.ecef);
    if (haveDopVel) sol.velNed = ecefToNed(dopVel, sol.llh);
    sol.ns = usableCount;
    sol.q = fixed ? 1 : 2;
    sol.ratio = ratio;
    sol.mode = fixed ? (heldFixed ? "RTK_FIX_HOLD" : "RTK_FIX") : "RTK_FLOAT";
    if (fixed || useDdCodeFloat) {
        flt.x[0] = sol.ecef.x;
        flt.x[1] = sol.ecef.y;
        flt.x[2] = sol.ecef.z;
    }
    flt.fixedHoldAge = fixed ? 0 : std::min(1000000, flt.fixedHoldAge + 1);

    const double codeRms = nCodeUpdate > 0 ? std::sqrt(codeSum2 / nCodeUpdate) : 0.0;
    const double phaseRms = nPhaseUpdate > 0 ? std::sqrt(phaseSum2 / nPhaseUpdate) : 0.0;
    if ((!fixed && (codeRms > 80.0 || phaseRms > 20.0)) ||
        (fixed && (codeRms > 120.0 || phaseRms > 20.0))) {
        if (debug) *debug << " fail=postfit_quality codeRms=" << codeRms
                          << " phaseRms=" << phaseRms << " usable=" << usableCount << "\n";
        flt = prior;
        return false;
    }
    if (flt.hasLastAccepted) {
        const double dt = std::max(1.0, std::abs(timediff(sol.time, flt.lastAcceptedTime)));
        const Vec3 dn = ecefToNed(sol.ecef - flt.lastAcceptedPos, ecefToLlh(flt.lastAcceptedPos));
        const double horiz = std::sqrt(dn.x * dn.x + dn.y * dn.y);
        const double vert = std::abs(dn.z);
        const double horizLimit = fixed ? 30.0 * dt : 8.0 * dt;
        const double vertLimit = fixed ? 30.0 * dt : 4.0 * dt;
        if (horiz > horizLimit || vert > vertLimit || sol.llh.h < -20.0 || sol.llh.h > 80.0) {
            if (debug) *debug << " fail=motion_quality horiz=" << horiz
                              << " vert=" << vert << " h=" << sol.llh.h << "\n";
            flt = prior;
            return false;
        }
    }

    sol.stdN = sol.stdE = fixed ? std::max(0.01, bestFixRms) : std::min(5.0, std::max(0.2, codeRms));
    sol.stdD = fixed ? std::max(0.02, bestFixRms * 2.0) : std::min(10.0, std::max(0.5, codeRms * 2.0));
    flt.acceptedSolution = true;
    flt.hasLastAccepted = true;
    flt.lastAcceptedTime = sol.time;
    flt.lastAcceptedPos = sol.ecef;

    if (debug) {
        *debug << " ok=" << sol.mode << " ratio=" << ratio
               << " bl=" << sol.baseline.x << "," << sol.baseline.y << "," << sol.baseline.z
               << " codeRms=" << codeRms << " phaseRms=" << phaseRms
               << " dopVel=" << (haveDopVel ? 1 : 0) << "/" << dopUsed << "/" << dopRms
               << " ddCode=" << (haveDdCode ? 1 : 0) << "/" << ddCodeRows << "/" << ddCodeRms
               << " fltH=" << ecefToLlh(fltPos).h << " ddH=" << (haveDdCode ? ecefToLlh(ddCodePos).h : 0.0)
               << " anchor=" << (anchorToCode ? 1 : 0) << "/" << codeAnchorDist
               << " nfix=" << nfixTotal << " fixRms=" << bestFixRms
               << " fixMax=" << bestFixMax << " fixRows=" << bestFixRows
               << " pos=" << sol.ecef.x << "," << sol.ecef.y << "," << sol.ecef.z << "\n";
    }
    return true;
}

// === Single-Difference Kalman Filter RTK ===
// State: pos(3), vel(3), clkGps(1), clkBds(1), amb_sd(N)
// SD observations: rover - base, with receiver clock difference estimation

static bool solveSdEpoch(const ObsEpoch &re, const ObsEpoch &be, const NavData &nav,
                         RtkFilterState &flt, const RtkOptions &opt, GnssSolution &sol,
                         const Vec3 &sppPos, std::ostream *debug = nullptr) {
    // Collect SD observations
    const Vec3 roverApprox = flt.initialized ? Vec3{flt.x[0], flt.x[1], flt.x[2]} : sppPos;
    std::vector<DdSat> sats = collectDdSats(re, be, nav, roverApprox, opt);
    int gpsF[2] = {}, bdsF[2] = {};
    for (const auto &s : sats) {
        for (int f = 0; f < 2; ++f) {
            if (!s.ok[f]) continue;
            if (s.sat.sys == System::GPS) ++gpsF[f];
            if (s.sat.sys == System::BDS) ++bdsF[f];
        }
    }
    // Count usable SD pairs (satellites with at least one freq valid)
    int usableCount = 0;
    for (const auto &s : sats) if (s.ok[0] || s.ok[1]) ++usableCount;

    if (debug) *debug << re.time.tow << " base_tow=" << be.time.tow << " dt=" << timediff(be.time, re.time)
                      << " common_usable=" << usableCount
                      << " gps=" << gpsF[0] << "/" << gpsF[1] << " bds=" << bdsF[0] << "/" << bdsF[1];

    if (usableCount < 4) {
        if (debug) *debug << " fail=not_enough_common\n";
        return false;
    }

    RtkFilterState prior = flt;

    // Initialize filter if first epoch
    if (!flt.initialized) {
        Vec3 codePos = sppPos;
        // Try to refine with DD code position
        double crms; int crows;
        Vec3 refinedPos = sppPos;
        if (solveDdCodePosition(re, be, nav, flt.basePosition, sats, opt, refinedPos, crms, crows)) {
            codePos = refinedPos;
        }
        initRtkFilter(flt, re.time, codePos, flt.basePosition);
        if (debug) *debug << " init_pos=" << codePos.x << "," << codePos.y << "," << codePos.z;
    }

    // Predict filter to current epoch
    predictRtkFilter(flt, re.time, opt);

    // Clear slip flags for this epoch
    flt.ambSlipFlag.clear();

    // === Sequential SD measurement updates ===
    int nCodeUpdate = 0, nPhaseUpdate = 0, nSlipReset = 0;
    double codeRmsSum = 0.0, phaseRmsSum = 0.0;

    for (int i = 0; i < static_cast<int>(sats.size()); ++i) {
        const DdSat &d = sats[i];
        for (int f = 0; f < 2; ++f) {
            if (!d.ok[f]) continue;

            const double lam = carrierWavelength(d.sat, f);
            const std::string ambKey = ambStateKey(d.sat, f);

            // Compute geometric SD range
            double rhoRov = 0.0, rhoBase = 0.0;
            Vec3 los;
            {
                Vec3 tmp;
                if (!observedGeometry(nav, d.sat, re.time, d.roverCode[f],
                                      {flt.x[0], flt.x[1], flt.x[2]}, rhoRov, los)) continue;
                if (!observedGeometry(nav, d.sat, be.time, d.baseCode[f],
                                      flt.basePosition, rhoBase, tmp)) continue;
            }
            const double rhoSd = rhoRov - rhoBase;

            // ---- SD Code Update ----
            {
                const double obsSdCode = d.sdCode[f];
                const double predicted = rhoSd + (d.sat.sys == System::GPS ? flt.x[kClkGpsIdx] : flt.x[kClkBdsIdx]);
                const double innovation = obsSdCode - predicted;

                // H: [los(3), 0(3), clk, 0...]
                std::vector<double> h(flt.x.size(), 0.0);
                h[0] = -los.x; h[1] = -los.y; h[2] = -los.z;
                if (d.sat.sys == System::GPS) h[kClkGpsIdx] = 1.0;
                else h[kClkBdsIdx] = 1.0;

                double codeVar = sqr(opt.codeSigma) / sqr(std::sin(std::max(d.elev, 10.0) * kD2R));
                if (scalarKalmanUpdate(flt, h, innovation, codeVar)) {
                    ++nCodeUpdate;
                    codeRmsSum += innovation * innovation;
                }
            }

            // ---- SD Phase Update ----
            {
                // Ensure SD ambiguity state
                const double clkSd0 = d.sat.sys == System::GPS ? flt.x[kClkGpsIdx] : flt.x[kClkBdsIdx];
                const int ambIdx = ensureAmbState(flt, ambKey, (d.sdPhaseM[f] - rhoSd - clkSd0) / lam);
                if (ambIdx < 0) continue;

                // Cycle slip detection via LLI
                const int roverLli = d.roverLli[f];
                const int baseLli = d.baseLli[f];
                bool slip = false;

                // LLI bit 0 rising edge ? cycle slip
                auto itR = flt.lastRoverLli.find(ambKey);
                if (itR != flt.lastRoverLli.end()) {
                    if ((roverLli & 1) && !(itR->second & 1)) slip = true;
                    if ((roverLli & 2) != (itR->second & 2)) slip = true;
                }
                flt.lastRoverLli[ambKey] = roverLli;

                // Also check base LLI for the same sat+freq
                auto itB = flt.lastBaseLli.find(ambKey);
                int prevBaseLli = (itB != flt.lastBaseLli.end()) ? itB->second : 0;
                if ((baseLli & 1) && !(prevBaseLli & 1)) slip = true;
                flt.lastBaseLli[ambKey] = baseLli;

                if (slip) {
                    resetAmbState(flt, ambIdx, 0.0, 100.0);
                    flt.ambSlipFlag[ambIdx] = 1;
                    ++nSlipReset;
                }

                const double obsSdPhaseM = d.sdPhaseM[f];
                const double predicted = rhoSd + (d.sat.sys == System::GPS ? flt.x[kClkGpsIdx] : flt.x[kClkBdsIdx])
                                        + lam * flt.x[ambIdx];
                const double innovation = obsSdPhaseM - predicted;

                std::vector<double> h(flt.x.size(), 0.0);
                h[0] = -los.x; h[1] = -los.y; h[2] = -los.z;
                if (d.sat.sys == System::GPS) h[kClkGpsIdx] = 1.0;
                else h[kClkBdsIdx] = 1.0;
                h[ambIdx] = lam;

                double phaseVar = sqr(opt.phaseSigma) / sqr(std::sin(std::max(d.elev, 10.0) * kD2R));
                if (scalarKalmanUpdate(flt, h, innovation, phaseVar)) {
                    ++nPhaseUpdate;
                    phaseRmsSum += innovation * innovation;

                    // Increment lock count
                    auto lit = flt.ambLock.find(ambIdx);
                    if (lit != flt.ambLock.end()) lit->second++;
                    else flt.ambLock[ambIdx] = 1;
                }
            }
        }
    }

    if (debug) *debug << " codeUpd=" << nCodeUpdate << " phaseUpd=" << nPhaseUpdate << " slipResets=" << nSlipReset;

    if (nPhaseUpdate < 3) {
        if (debug) *debug << " fail=not_enough_phase_updates\n";
        return false;
    }

    // === LAMBDA AR on DD formed from SD filter states ===
    // Collect SD ambiguity states with sufficient lock
    struct SdAmbInfo {
        int stateIdx;
        SatId sat;
        int freq;
        double floatVal;
        double elev;
    };
    std::vector<SdAmbInfo> sdList;
    for (int i = 0; i < static_cast<int>(sats.size()); ++i) {
        for (int f = 0; f < 2; ++f) {
            if (!sats[i].ok[f]) continue;
            const std::string ambKey = ambStateKey(sats[i].sat, f);
            auto ait = flt.ambIndex.find(ambKey);
            if (ait == flt.ambIndex.end()) continue;
            const int idx = ait->second;
            if (idx >= static_cast<int>(flt.x.size())) continue;
            auto lit = flt.ambLock.find(idx);
            const int lock = (lit != flt.ambLock.end()) ? lit->second : 0;
            if (lock < opt.minAmbLock) continue;

            SdAmbInfo info;
            info.stateIdx = idx;
            info.sat = sats[i].sat;
            info.freq = f;
            info.floatVal = flt.x[idx];
            info.elev = sats[i].elev;
            sdList.push_back(info);
        }
    }

    // Group by (sys, freq), pick reference (highest elevation)
    std::map<std::pair<System, int>, std::vector<int>> groups;
    for (int i = 0; i < static_cast<int>(sdList.size()); ++i) {
        groups[std::make_pair(sdList[i].sat.sys, sdList[i].freq)].push_back(i);
    }

    bool fixed = false;
    double ratio = 0.0;
    int nfixTotal = 0;
    Vec3 fixedPos = {flt.x[0], flt.x[1], flt.x[2]};

    // For each group, form DD and try LAMBDA
    for (auto &gkv : groups) {
        auto &indices = gkv.second;
        if (indices.size() < 2) continue;

        // Pick reference: highest elevation
        int refIdx = indices[0];
        for (int idx : indices) if (sdList[idx].elev > sdList[refIdx].elev) refIdx = idx;

        // Form DD
        std::vector<double> ddFloat;
        int ndd = 0;
        std::vector<int> ddSatGlobal, ddRefGlobal;
        int ddFreq = sdList[refIdx].freq;

        for (int idx : indices) {
            if (idx == refIdx) continue;
            // DD = N_sat - N_ref
            double val = sdList[idx].floatVal - sdList[refIdx].floatVal;
            ddFloat.push_back(val);
            ddSatGlobal.push_back(sdList[idx].stateIdx);
            ddRefGlobal.push_back(sdList[refIdx].stateIdx);
            ++ndd;
        }

        if (ndd < 4) continue;

        // Build DD covariance Qnn from filter P
        Matrix qnn(ndd, ndd);
        for (int r = 0; r < ndd; ++r) {
            const int si = ddSatGlobal[r];
            const int ri = ddRefGlobal[r];
            for (int c = 0; c < ndd; ++c) {
                const int sj = ddSatGlobal[c];
                const int rj = ddRefGlobal[c];
                qnn(r, c) = flt.p(si, sj) - flt.p(si, rj) - flt.p(ri, sj) + flt.p(ri, rj);
            }
        }

        // LAMBDA
        std::vector<int> nint;
        double rr = 0.0;
        if (!lambdaIntegerSearch(ddFloat, qnn, nint, rr)) continue;

        if (debug) *debug << " ar_" << (gkv.first.first == System::GPS ? "G" : "C")
                          << gkv.first.second << " ndd=" << ndd << " ratio=" << rr;

        if (rr < opt.lambdaRatioThreshold) continue;

        // Compute fixed position from fixed DD phase measurements
        std::vector<std::vector<double>> posRows;
        std::vector<double> posRhs;
        Vec3 curPos = {flt.x[0], flt.x[1], flt.x[2]};

        for (int k = 0; k < ndd; ++k) {
            // Find the sat and ref in sats
            int satIdxInSats = -1, refIdxInSats = -1;
            for (int i2 = 0; i2 < static_cast<int>(sats.size()); ++i2) {
                for (int f2 = 0; f2 < 2; ++f2) {
                    if (!sats[i2].ok[f2]) continue;
                    const std::string k2 = ambStateKey(sats[i2].sat, f2);
                    auto it2 = flt.ambIndex.find(k2);
                    if (it2 == flt.ambIndex.end()) continue;
                    if (it2->second == ddSatGlobal[k]) satIdxInSats = i2;
                    if (it2->second == ddRefGlobal[k]) refIdxInSats = i2;
                }
            }
            if (satIdxInSats < 0 || refIdxInSats < 0) continue;

            const DdSat &s = sats[satIdxInSats];
            const DdSat &r = sats[refIdxInSats];
            double Nfix = static_cast<double>(nint[k]);
            double lam = carrierWavelength(s.sat, ddFreq);

            double geomDd = 0.0;
            Vec3 losS, losR;
            if (!ddGeometry(re, be, nav, curPos, flt.basePosition, s, r, ddFreq, geomDd, losS, losR)) continue;

            double obsDd = s.sdPhaseM[ddFreq] - r.sdPhaseM[ddFreq];
            double phaseRange = obsDd - lam * Nfix;
            Vec3 a = losR - losS;
            posRows.push_back({a.x, a.y, a.z});
            posRhs.push_back(phaseRange - geomDd);
        }

        if (static_cast<int>(posRows.size()) >= 3) {
            Matrix h(static_cast<int>(posRows.size()), 3);
            for (int r = 0; r < h.rows(); ++r)
                for (int c = 0; c < 3; ++c) h(r, c) = posRows[r][c];
            std::vector<double> dx;
            if (leastSquares(h, posRhs, dx, nullptr)) {
                Vec3 candPos = {curPos.x + dx[0], curPos.y + dx[1], curPos.z + dx[2]};
                if (norm(candPos - curPos) < 2.0) {
                    fixedPos = candPos;
                    fixed = true;
                    nfixTotal = ndd;
                    if (rr > ratio) ratio = rr;
                }
            }
        }
    }

    // === Output ===
    sol.time = re.time;
    Vec3 fltPos = {flt.x[0], flt.x[1], flt.x[2]};
    Vec3 ddCodePos = fltPos;
    double ddCodeRms = 0.0;
    int ddCodeRows = 0;
    const bool haveDdCode = solveDdCodePosition(re, be, nav, flt.basePosition, sats, opt, ddCodePos, ddCodeRms, ddCodeRows);
    const bool goodDdCode = haveDdCode && ddCodeRows >= 8 && ddCodeRms < 5.0;
    if (fixed && goodDdCode && norm(fixedPos - ddCodePos) > 2.0) {
        fixed = false;
        ratio = 0.0;
        nfixTotal = 0;
    }
    const bool useDdCodePos = !fixed && goodDdCode;
    sol.ecef = fixed ? fixedPos : (useDdCodePos ? ddCodePos : fltPos);
    // Guard against filter divergence; keep fallback independent and reset RTK state.
    if (!std::isfinite(sol.ecef.x) || !std::isfinite(sol.ecef.y) || !std::isfinite(sol.ecef.z) ||
        norm(sol.ecef - sppPos) > 500.0) {
        if (debug) *debug << " fail=rtk_diverged\n";
        if (prior.acceptedSolution) flt = prior;
        return false;
    }
    sol.baseline = sol.ecef - flt.basePosition;
    sol.llh = ecefToLlh(sol.ecef);
    sol.ns = usableCount;
    sol.q = fixed ? 1 : 2;
    sol.ratio = ratio;
    sol.mode = fixed ? "RTK_FIX" : "RTK_FLOAT";

    // Update filter position state if fixed
    if (fixed || useDdCodePos) {
        flt.x[0] = sol.ecef.x; flt.x[1] = sol.ecef.y; flt.x[2] = sol.ecef.z;
    }

    double codeRms = nCodeUpdate > 0 ? std::sqrt(codeRmsSum / nCodeUpdate) : 0.0;
    double phaseRms = nPhaseUpdate > 0 ? std::sqrt(phaseRmsSum / nPhaseUpdate) : 0.0;
    if ((!fixed && (usableCount < 8 || codeRms > 80.0 || phaseRms > 20.0)) ||
        (fixed && (codeRms > 120.0 || phaseRms > 30.0))) {
        if (debug) {
            *debug << " fail=postfit_quality codeRms=" << codeRms
                   << " phaseRms=" << phaseRms << " usable=" << usableCount << "\n";
        }
        if (prior.acceptedSolution) flt = prior;
        return false;
    }
    if (flt.hasLastAccepted) {
        const double dt = std::max(1.0, std::abs(timediff(sol.time, flt.lastAcceptedTime)));
        const Vec3 dn = ecefToNed(sol.ecef - flt.lastAcceptedPos, ecefToLlh(flt.lastAcceptedPos));
        const double horiz = std::sqrt(dn.x * dn.x + dn.y * dn.y);
        const double vert = std::abs(dn.z);
        if (horiz > 8.0 * dt || vert > 4.0 * dt || sol.llh.h < 0.0 || sol.llh.h > 40.0) {
            if (debug) {
                *debug << " fail=motion_quality horiz=" << horiz
                       << " vert=" << vert << " h=" << sol.llh.h << "\n";
            }
            if (prior.acceptedSolution) flt = prior;
            return false;
        }
    }
    sol.stdN = sol.stdE = std::min(5.0, std::max(0.02, codeRms));
    sol.stdD = std::min(10.0, std::max(0.04, codeRms * 2.0));
    flt.acceptedSolution = true;
    flt.hasLastAccepted = true;
    flt.lastAcceptedTime = sol.time;
    flt.lastAcceptedPos = sol.ecef;

    if (debug) {
        *debug << " ok=" << sol.mode << " ratio=" << ratio
               << " bl=" << sol.baseline.x << "," << sol.baseline.y << "," << sol.baseline.z
               << " codeRms=" << codeRms << " phaseRms=" << phaseRms
               << " ddCode=" << (haveDdCode ? 1 : 0) << "/" << ddCodeRows << "/" << ddCodeRms
               << " ddPos=" << ddCodePos.x << "," << ddCodePos.y << "," << ddCodePos.z
               << " nfix=" << nfixTotal << " pos=" << sol.ecef.x << "," << sol.ecef.y << "," << sol.ecef.z << "\n";
    }
    return true;
}

// 通过毫秒级周内秒索引查找同名历元，加快流动站和基准站配对。
static const ObsEpoch *findEpochByTow(const std::map<int, const ObsEpoch *> &epochs, const GTime &t) {
    const int key = static_cast<int>(std::llround(t.tow * 1000.0));
    auto it = epochs.find(key);
    return it == epochs.end() ? nullptr : it->second;
}

// 执行一次 RTK Kalman 滤波 pass，可按正向或反向历元顺序处理。
static std::vector<GnssSolution> runRtkPass(const RinexObsFile &rover, const RinexObsFile &base,
                                            const NavData &nav, const RtkOptions &opt,
                                            bool reverse, const std::string &debugPath) {
    // Run one RTK filtering pass.
    // reverse=false: process epochs from start to end.
    // reverse=true : process epochs from end to start, then reverse the result outside.
    //
    // Each epoch first gets an SPP fallback solution. If a matching base epoch exists,
    // solveDdEpoch tries to replace it by float/fixed RTK.
    std::vector<GnssSolution> out;
    Vec3 last = rover.approxPos;
    bool hasLastVelRef = false;
    GnssSolution lastVelRef;

    // Single filter for GPS+BDS combined SD EKF
    RtkFilterState flt;
    flt.basePosition = base.approxPos;

    std::map<int, const ObsEpoch *> baseByTow;
    for (const auto &be : base.epochs) {
        baseByTow[static_cast<int>(std::llround(be.time.tow * 1000.0))] = &be;
    }

    std::ofstream debug(debugPath);
    if (debug) debug << "# rover_tow base_tow dt common_usable codeUpd phaseUpd slipResets ar_sys ndd ratio status pos\n";

    const int nEpoch = static_cast<int>(rover.epochs.size());
    for (int ep = 0; ep < nEpoch; ++ep) {
        const int idx = reverse ? (nEpoch - 1 - ep) : ep;
        const ObsEpoch &re = rover.epochs[idx];
        const ObsEpoch *be = findEpochByTow(baseByTow, re.time);

        // SPP for fallback and initial position
        GnssSolution spp;
        Vec3 initPos = last;
        if (!sppSolve(re, nav, initPos, spp)) continue;

        GnssSolution sol = spp;
        sol.mode = "PPP_FALLBACK";
        sol.q = 5;
        sol.baseline = sol.ecef - base.approxPos;

        if (be) {
            GnssSolution rtk;
            if (solveDdEpoch(re, *be, nav, flt, opt, rtk, spp.ecef, debug ? &debug : nullptr)) {
                sol = rtk;
            }
        }

            if (norm(sol.velNed) < 1e-6 && hasLastVelRef) {
                const double dt = timediff(sol.time, lastVelRef.time);
                if (dt > 0.0 && dt < 5.0) {
                    const Vec3 dn = ecefToNed(sol.ecef - lastVelRef.ecef, lastVelRef.llh);
                sol.velNed = dn / dt;
            }
        }
        if (sol.q != 5) {
            lastVelRef = sol;
            hasLastVelRef = true;
        }
        out.push_back(sol);
        last = sol.ecef;
    }
    return out;
}

// 合并前向和后向 RTK 结果：质量好的一侧优先，两侧都好时按方差加权平均。
static GnssSolution combineForwardBackward(const GnssSolution &f, const GnssSolution &b, const Vec3 &basePos) {
    if (f.time.week != b.time.week || std::abs(f.time.tow - b.time.tow) > 1e-4) return f;
    if (f.q == 5 && b.q != 5) return b;
    if (b.q == 5 && f.q != 5) return f;
    if (f.q == 5 && b.q == 5) return f;

    const Vec3 fbDn = ecefToNed(f.ecef - b.ecef, f.llh);
    const double fbSep = norm(fbDn);
    auto demoteUnconfirmedFix = [](GnssSolution s) {
        s.q = 2;
        s.ratio = 0.0;
        s.stdN = std::max(s.stdN, 0.5);
        s.stdE = std::max(s.stdE, 0.5);
        s.stdD = std::max(s.stdD, 1.0);
        s.mode = "RTK_FLOAT_UNCONFIRMED_FIX";
        return s;
    };
    auto reliableOneSidedFix = [](const GnssSolution &s) {
        return s.ratio >= 30.0 && s.stdN <= 0.03 && s.stdE <= 0.03 && s.stdD <= 0.06;
    };
    if (f.q == 1 && b.q != 1) return (fbSep < 1.0 || (fbSep < 1.2 && reliableOneSidedFix(f))) ? f : demoteUnconfirmedFix(f);
    if (b.q == 1 && f.q != 1) return (fbSep < 1.0 || (fbSep < 1.2 && reliableOneSidedFix(b))) ? b : demoteUnconfirmedFix(b);

    const double sf = std::max(0.02, (f.stdN + f.stdE + f.stdD) / 3.0);
    const double sb = std::max(0.02, (b.stdN + b.stdE + b.stdD) / 3.0);
    const double wf = 1.0 / (sf * sf);
    const double wb = 1.0 / (sb * sb);
    GnssSolution c = f;
    c.ecef = (f.ecef * wf + b.ecef * wb) / (wf + wb);
    c.baseline = c.ecef - basePos;
    c.llh = ecefToLlh(c.ecef);
    c.q = (f.q == 1 && b.q == 1 && fbSep < 0.5) ? 1 : 2;
    c.mode = c.q == 1 ? "RTK_FIX_COMBINED" : "RTK_FLOAT_COMBINED";
    c.ratio = std::max(f.ratio, b.ratio);
    c.ns = std::max(f.ns, b.ns);
    c.stdN = c.stdE = std::max(0.02, std::sqrt(1.0 / (wf + wb)));
    c.stdD = std::max(0.04, c.stdN * 2.0);
    return c;
}

// SD Kalman Filter RTK: processes rover+base epochs with forward/backward combined filtering.
std::vector<GnssSolution> runRtkFloat(const RinexObsFile &rover, const RinexObsFile &base, const NavData &nav,
                                      const RtkOptions &opt) {
    // RTK 主入口：
    // 内部先分别做前向/后向 Kalman RTK 处理，再组合得到较稳定的 GNSS 历元解。
    // 输出的 GnssSolution 作为松组合的 1 Hz 位置/速度观测。
    // 注意：这里是 GNSS 预处理，pseudo-online 的“不用未来数据”约束针对组合导航滤波循环，
    // 组合导航自身只按当前到达的 GNSS 历元进行 update。
    std::vector<GnssSolution> fwd = runRtkPass(rover, base, nav, opt, false,
                                              "output/rtk_debug.txt");
    std::vector<GnssSolution> bwd = runRtkPass(rover, base, nav, opt, true,
                                              "output/rtk_debug_bwd.txt");
    std::reverse(bwd.begin(), bwd.end());
    if (fwd.size() != bwd.size()) return fwd;
    std::vector<GnssSolution> out;
    out.reserve(fwd.size());
    bool hasLastVelRef = false;
    GnssSolution lastVelRef;
    for (size_t i = 0; i < fwd.size(); ++i) {
        GnssSolution sol = combineForwardBackward(fwd[i], bwd[i], base.approxPos);
        if (norm(sol.velNed) < 1e-6 && hasLastVelRef) {
            const double dt = timediff(sol.time, lastVelRef.time);
            if (dt > 0.0 && dt < 5.0) {
                const Vec3 dn = ecefToNed(sol.ecef - lastVelRef.ecef, lastVelRef.llh);
                sol.velNed = dn / dt;
            }
        }
        if (sol.q != 5) {
            lastVelRef = sol;
            hasLastVelRef = true;
        }
        out.push_back(sol);
    }
    bool hasLastFix = false;
    GnssSolution lastFix;
    for (auto &s : out) {
        if (s.q != 1) continue;
        if (s.stdD > 0.06) {
            s.q = 2;
            s.ratio = 0.0;
            s.stdN = std::max(s.stdN, 0.5);
            s.stdE = std::max(s.stdE, 0.5);
            s.stdD = std::max(s.stdD, 1.0);
            s.mode = "RTK_FLOAT_FIX_STD_REJECT";
            continue;
        }
        if (hasLastFix) {
            const double dt = timediff(s.time, lastFix.time);
            const double maxDh = 0.5 + 0.02 * std::max(0.0, dt);
            if (dt > 0.0 && std::abs(s.llh.h - lastFix.llh.h) > maxDh) {
                s.q = 2;
                s.ratio = 0.0;
                s.stdN = std::max(s.stdN, 0.5);
                s.stdE = std::max(s.stdE, 0.5);
                s.stdD = std::max(s.stdD, 1.0);
                s.mode = "RTK_FLOAT_FIX_HEIGHT_REJECT";
                continue;
            }
        }
        lastFix = s;
        hasLastFix = true;
    }
    return out;
}

// 将 RTK 解算结果写入文本文件，供组合导航和可视化程序继续读取。
void writeGnssSolutions(const std::string &path, const std::vector<GnssSolution> &sols) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write gnss solution: " + path);
    out << "# week tow lat_deg lon_deg h bx by bz vn ve vd q ns ratio stdN stdE stdD mode\n";
    out << std::fixed << std::setprecision(9);
    for (const auto &s : sols) {
        out << s.time.week << ' ' << std::setprecision(3) << s.time.tow << ' '
            << std::setprecision(12) << s.llh.lat * kR2D << ' ' << s.llh.lon * kR2D << ' '
            << std::setprecision(4) << s.llh.h << ' ' << s.baseline.x << ' ' << s.baseline.y << ' '
            << s.baseline.z << ' ' << s.velNed.x << ' ' << s.velNed.y << ' ' << s.velNed.z
            << ' ' << s.q << ' ' << s.ns << ' ' << std::setprecision(3) << s.ratio << ' '
            << s.stdN << ' ' << s.stdE << ' ' << s.stdD << ' ' << s.mode << '\n';
    }
}

} // namespace gnssins


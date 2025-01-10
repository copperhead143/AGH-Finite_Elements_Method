#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <numeric>


#define PLIK 2
#define NPC 2
#define NPCBC 2


#if PLIK == 1
std::string file_path = "C:/Users/ssiko/CLionProjects/MES_refactor/Test1_4_4.txt";
#endif

#if PLIK == 2
std::string file_path = "C:/Users/ssiko/CLionProjects/MES_refactor/Test2_4_4_MixGrid.txt";
#endif

#if PLIK == 3
std::string file_path = "C:/Users/ssiko/CLionProjects/MES_refactor/Test3_31_31_kwadrat.txt";
#endif

#if NPC == 2
std::vector<double> ksi = { -1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0), -1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0) };
std::vector<double> eta = { -1.0 / std::sqrt(3.0), -1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0) };
std::vector<double> weights{ 1.0, 1.0 };
#endif
#if NPC == 3
std::vector<double> ksi = { -std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0),
                           -std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0),
                           -std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0) };
std::vector<double> eta = { -std::sqrt(3.0 / 5.0), -std::sqrt(3.0 / 5.0), -std::sqrt(3.0 / 5.0),
                           0.0, 0.0, 0.0,
                           std::sqrt(3.0 / 5.0), std::sqrt(3.0 / 5.0), std::sqrt(3.0 / 5.0) };
std::vector<double> weights = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
#endif
#if NPC == 4
std::vector<double> ksi = { -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                           -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                           -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                           -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                           -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                           -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                           -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                           -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
};
std::vector<double> eta = { -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
                            std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
};
std::vector<double> weights = { (18.0 - std::sqrt(30.0)) / 36.0,
                               (18.0 + std::sqrt(30.0)) / 36.0,
                               (18.0 + std::sqrt(30.0)) / 36.0,
                               (18.0 - std::sqrt(30.0)) / 36.0 };
#endif

#if NPCBC == 2
std::vector<double> ksibc = { -1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0), 1.0, 1.0, 1.0 / std::sqrt(3.0), -1.0 / std::sqrt(3.0), -1.0, -1.0 };
std::vector<double> etabc = { -1.0, -1.0, -1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0), 1.0, 1.0, 1.0 / std::sqrt(3.0), -1.0 / std::sqrt(3.0) };
std::vector<double> weightsbc = { 1.0, 1.0 };
#endif
#if NPCBC == 3
std::vector<double> ksibc = { -std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0),
                              1.0, 1.0, 1.0,
                              std::sqrt(3.0 / 5.0), 0.0, -std::sqrt(3.0 / 5.0),
                              -1.0, -1.0, -1.0 };
std::vector<double> etabc = { -1.0, -1.0, -1.0,
                              -std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0),
                               1.0, 1.0, 1.0,
                               std::sqrt(3.0 / 5.0), 0.0, -std::sqrt(3.0 / 5.0) };
std::vector<double> weightsbc = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
#endif
#if NPCBC == 4
std::vector<double> ksibc = {
    -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
     std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
     std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
     1.0, 1.0, 1.0, 1.0,
     std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
     std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    -1.0, -1.0, -1.0, -1.0
};

std::vector<double> etabc = {
    -1.0, -1.0, -1.0, -1.0,
    -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
     std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
     std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
     1.0, 1.0, 1.0, 1.0,
     std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
     std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0))
};

std::vector<double> weightsbc = {
    (18.0 - std::sqrt(30.0)) / 36.0,
    (18.0 + std::sqrt(30.0)) / 36.0,
    (18.0 + std::sqrt(30.0)) / 36.0,
    (18.0 - std::sqrt(30.0)) / 36.0
};
#endif


struct GlobalData {
    int SimulationTime;
    int SimulationStepTime;
    int Conductivity;
    int Alfa;
    int Tot;
    int InitialTemp;
    int Density;
    int SpecificHeat;
    int NodesNumber;
    int ElementsNumber;
};

struct Node {
  double x;
  double y;
  bool BC = false;
};

struct UnivElement {
  std::vector<std::vector<double>> dNdKsi;
  std::vector<std::vector<double>> dNdEta;
  std::vector<std::vector<double>> dNdBC;
  std::vector<std::vector<double>> N;

 UnivElement();

  void calcedN_values();

  void calcSurface();

  void calcN();
};

struct Element {
  std::vector<int> nodeID;
  using Vector = std::vector<double>;
  using Matrix = std::vector<std::vector<double>>;
  Matrix H_mat;
  Matrix Hbc;
  Vector P_vec;
  Matrix C_mat;
  std::vector<Matrix> Jac;
  std::vector<Matrix> Jac1;
  Vector detJ;
  Vector detJbc;

  void calcJakobian();
  void calcDetJ();
  void invJakobian();
  void calcH();

  void calcDetJBc();
  void calcHbc();

  void calcP();

  void calcC();

  void setup();
};

struct Solution{
  using Matrix = std::vector<std::vector<double>>;
  using Vector = std::vector<double>;
  Matrix Hglobal;
  Vector Pglobal;
  Matrix Cglobal;
  Matrix Global;

  void setup();
  void agregate();
  void printP();
  void printH();
  void printC();
  void solve();
  Vector Gauss(const Matrix& A, Vector& b);
};

std::vector<Node> loadNodes(const std::string& filename) {
    std::vector<Node> nodes;
    std::ifstream file(filename);
    std::string line;
    bool nodeSection = false;
    bool bcSection = false;

    while (getline(file, line)) {
        if (line.find("*Node") != std::string::npos) {
            nodeSection = true;
            continue;
        }
        if (line.find("*BC") != std::string::npos) {
            bcSection = true;
            nodeSection = false;
            continue;
        }
        if (line.find('*') != std::string::npos) {
            nodeSection = false;
            bcSection = false;
        }
        if (nodeSection && !line.empty()) {
            std::istringstream iss(line);
            std::string temp;
            double x, y;

            getline(iss, temp, ',');
            iss >> x;
            iss.ignore(1, ',');
            iss >> y;
            if (std::isnan(x) || std::isnan(y)) {
                std::cerr << "Error: NaN in node coordinates. Line: " << line << std::endl;
            }
            nodes.push_back({x, y});
        }
        if (bcSection && !line.empty()) {
            std::istringstream iss(line);
            std::string temp;

            while (getline(iss, temp, ',')) {
                int nodeId = std::stoi(temp) - 1;
                if (nodeId >= 0 && nodeId < nodes.size()) {
                    nodes[nodeId].BC = true;
                } else {
                    std::cerr << "Error: Invalid node ID in BC section. Line: " << line << std::endl;
                }
            }
        }
    }

    return nodes;
}
GlobalData loadData(const std::string& filename) {
    GlobalData data;
    std::ifstream file (filename);
    std::string line;
    while (getline(file, line)) {
        std::istringstream iss(line);
        std::string key, extra;
        int value;
        if (!(iss >> key)) continue;
        if ((key == "Nodes" || key == "Elements") && (iss >> extra) && (extra == "number")) {
            if (iss >> value) {
                if (key == "Nodes") {
                    data.NodesNumber = value;
                } else if (key == "Elements") {
                    data.ElementsNumber = value;
                }
            }
        }
        else if (iss >> value) {
            if (key == "SimulationTime") data.SimulationTime = value;
            else if (key == "SimulationStepTime") data.SimulationStepTime = value;
            else if (key == "Conductivity") data.Conductivity = value;
            else if (key == "Alfa") data.Alfa = value;
            else if (key == "Tot") data.Tot = value;
            else if (key == "InitialTemp") data.InitialTemp = value;
            else if (key == "Density") data.Density = value;
            else if (key == "SpecificHeat") data.SpecificHeat = value;
        }
    }
    return data;
}
std::vector<Element> loadElements(const std::string& filename) {
    std::vector<Element> elements;
    std::ifstream file(filename);
    std::string line;
    bool elementSection = false;

    while (getline(file, line)) {
        if (line.find("Element, type=DC2D4") != std::string::npos) {
            elementSection = true;
            continue;
        }

        if (elementSection && line.find("*BC") != std::string::npos) {
            break;
        }

        if (elementSection && !line.empty()) {
            std::istringstream iss(line);
            std::string temp;
            std::vector<int> nodeIds;

            getline(iss, temp, ',');

            while (getline(iss, temp, ',')) {
                nodeIds.push_back(stoi(temp));
            }

            elements.push_back({nodeIds});
        }
    }

    return elements;
}

std::vector<Node> nodes = loadNodes(file_path);
GlobalData data = loadData(file_path);
std::vector<Element> elements = loadElements(file_path);

UnivElement::UnivElement() {
    calcN();
    calcSurface();
    calcedN_values();
};

void UnivElement::calcedN_values() {
    for (int i = 0; i < NPC * NPC; ++i) {
        dNdKsi.push_back({-0.25 * (1 - eta[i]),
                          0.25 * (1 - eta[i]),
                          0.25 * (1 + eta[i]),
                          -0.25 * (1 + eta[i])});
        dNdEta.push_back({-0.25 * (1 - ksi[i]),
                          -0.25 * (1 + ksi[i]),
                          0.25 * (1 + ksi[i]),
                          0.25 * (1 - ksi[i])});
    }
}

void UnivElement::calcSurface() {
    for (int i = 0; i< NPCBC*4;i++) {
        dNdBC.push_back({0.25 * (1 - ksibc[i]) * (1 - etabc[i]),
                    0.25 * (1 + ksibc[i]) * (1 - etabc[i]),
                    0.25 * (1 + ksibc[i]) * (1 + etabc[i]),
                    0.25 * (1 - ksibc[i]) * (1 + etabc[i])});
    }
}

void UnivElement::calcN() {
    for (int i = 0; i < NPC * NPC; ++i) {
        N.push_back({0.25 * (1 - ksi[i]) * (1 - eta[i]),
                     0.25 * (1 + ksi[i]) * (1 - eta[i]),
                     0.25 * (1 + ksi[i]) * (1 + eta[i]),
                     0.25 * (1 - ksi[i]) * (1 + eta[i])});
    }
}

UnivElement univElement; //globala zmienna dla elementu uniwersalnego

void Element::setup() {
    calcJakobian();
    calcDetJ();
    invJakobian();
    calcH();
    calcDetJBc();
    calcHbc();
    calcP();
    calcC();
}

//obliczanie macierzy H
void Element::calcJakobian() {
    Jac.resize(NPC * NPC, Matrix(2,Vector(2, 0.0)));
    for (int pc = 0; pc < NPC * NPC; ++pc) {
        for (int j = 0; j < 4; ++j) {
            Jac[pc][0][0] += univElement.dNdKsi[pc][j] * nodes[nodeID[j] - 1].x;
            Jac[pc][0][1] += univElement.dNdKsi[pc][j] * nodes[nodeID[j] - 1].y;
            Jac[pc][1][0] += univElement.dNdEta[pc][j] * nodes[nodeID[j] - 1].x;
            Jac[pc][1][1] += univElement.dNdEta[pc][j] * nodes[nodeID[j] - 1].y;
        }
    }
}

void Element::calcDetJ() {
    for (int pc = 0; pc < NPC * NPC; pc++) {
        double determinant = Jac[pc][0][0] * Jac[pc][1][1] - Jac[pc][0][1] * Jac[pc][1][0];
        if (std::abs(determinant) < 1e-12) { // Jeśli wyznacznik jest zbyt mały lub zerowy
            std::cerr << "Warning: Small or zero determinant at integration point " << pc << std::endl;
        }
        detJ.push_back(determinant);
    }
}

void Element::invJakobian() {
    Jac1 = Jac;
    for (int pc = 0; pc < NPC * NPC; pc++) {
        Jac1[pc][0][0] = Jac[pc][1][1] / detJ[pc];
        Jac1[pc][0][1] = -Jac[pc][0][1] / detJ[pc];
        Jac1[pc][1][0] = -Jac[pc][1][0] / detJ[pc];
        Jac1[pc][1][1] = Jac[pc][0][0] / detJ[pc];
    }
}

void Element::calcH() {
    Matrix matrix_x(NPC * NPC, Vector(4, 0.0));
    Matrix matrix_y(NPC * NPC, Vector(4, 0.0));
    std::vector<Matrix> matrix_x_pc(NPC * NPC, Matrix(4, Vector(4, 0.0)));
    std::vector<Matrix> matrix_y_pc(NPC * NPC, Matrix(4, Vector(4, 0.0)));
    std::vector<Matrix> matrix_h_pc(NPC * NPC, Matrix(4, Vector(4, 0.0)));
    H_mat.resize(4, Vector(4, 0.0));

    for (int pc = 0; pc < NPC * NPC; ++pc) {
        for (int j = 0; j < 4; ++j) {
            matrix_x[pc][j] = Jac1[pc][0][0] * univElement.dNdKsi[pc][j] + Jac1[pc][0][1] * univElement.dNdEta[pc][j];
            matrix_y[pc][j] = Jac1[pc][1][0] * univElement.dNdKsi[pc][j] + Jac1[pc][1][1] * univElement.dNdEta[pc][j];
        }
    }

    for (int pc = 0; pc < NPC * NPC; ++pc) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                matrix_x_pc[pc][i][j] = matrix_x[pc][i] * matrix_x[pc][j];
                matrix_y_pc[pc][i][j] = matrix_y[pc][i] * matrix_y[pc][j];
            }
        }
    }

    for (int pc = 0; pc < NPC * NPC; ++pc) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                matrix_h_pc[pc][i][j] = data.Conductivity * (matrix_x_pc[pc][i][j] + matrix_y_pc[pc][i][j]) * detJ[pc];
            }
        }
    }

    for (int k = 0; k < NPC; ++k) {
        for (int l = 0; l < NPC; ++l) {
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    H_mat[i][j] += matrix_h_pc[k * NPC + l][i][j] * weights[k] * weights[l];
                }
            }
        }
    }
}


void Element::calcDetJBc() {
    Vector L;
    Vector x;
    Vector y;
    double dx;
    double dy;
    for (int i = 0; i<4;i++) {
        x.push_back(nodes[nodeID[i] - 1].x);
        y.push_back(nodes[nodeID[i]-1].y);
    }

    for (int i = 0; i < 4; ++i) {
        dx = nodes[nodeID[i]-1].x - nodes[nodeID[(i+1)%4] - 1].x;
        dy = nodes[nodeID[i]-1].y - nodes[nodeID[(i+1)%4] - 1].y;
        detJbc.push_back(sqrt(dx*dx + dy*dy)/2);
    }
}

void Element::calcHbc() {
    Matrix HbcSuma(4, Vector(16, 0.0));

    for (int edge = 0; edge < 4; ++edge) {
        for (int point = 0; point < NPCBC; ++point) {
            int index = edge * NPCBC + point;

            for (int row = 0; row < 4; ++row) {
                for (int col = 0; col < 4; ++col) {
                    HbcSuma[edge][row * 4 + col] += univElement.dNdBC[index][row] * univElement.dNdBC[index][col] * weightsbc[point];
                }
            }
        }

        for (int k = 0; k < 16; ++k) {
            HbcSuma[edge][k] *= data.Alfa * detJbc[edge];
        }
    }

    Hbc.clear();
    Hbc.resize(4, Vector(4, 0.0));

    for (int edge = 0; edge < 4; ++edge) {
        if (nodes[nodeID[edge] - 1].BC && nodes[nodeID[(edge + 1) % 4] - 1].BC) {
            for (int row = 0; row < 4; ++row) {
                for (int col = 0; col < 4; ++col) {
                    Hbc[row][col] += HbcSuma[edge][row * 4 + col];
                }
            }
        }
    }
}

void Element::calcP() {
    Matrix temp(NPCBC * 4, Vector(4, 0.0));

    for (int i = 0; i < NPCBC * 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            temp[i][j] = univElement.dNdBC[i][j] * data.Tot;
        }
    }

    Matrix Psum(4, Vector(4, 0.0));

    for (int edge = 0; edge < 4; ++edge) {
        for (int point = 0; point < NPCBC; ++point) {
            int index = edge * NPCBC + point;
            for (int i = 0; i < 4; ++i) {
                Psum[edge][i] += temp[index][i] * weightsbc[point];
            }
        }
    }

    P_vec.clear();
    P_vec.resize(4, 0.0);

    for (int edge = 0; edge < 4; ++edge) {
        if (nodes[nodeID[edge] - 1].BC && nodes[nodeID[(edge + 1) % 4] - 1].BC) {
            for (int i = 0; i < 4; ++i) {
                P_vec[i] += Psum[edge][i] * data.Alfa * detJbc[edge];
            }
        }
    }
}

void Element::calcC() {
    C_mat.resize(4, Vector(4, 0.0));

//zlozonosc obliczeniowa? Tak
    for (int i = 0; i < NPC; i++) {
        for (int j = 0; j < NPC; j++) {
            for (int k = 0; k < 4; k++) {
                for (int l = 0; l < 4; l++) {
                    C_mat[k][l] += (data.Density * data.SpecificHeat * detJ[i*NPC+j] * univElement.N[i*NPC+j][k] * univElement.N[i*NPC+j][l] * weights[i] * weights[j]);
                }
            }
        }
    }
}

void Solution::setup() {
    agregate();
}

void Solution::agregate() {
    Hglobal.resize(data.NodesNumber, Vector(data.NodesNumber, 0.0));
    Pglobal.resize(data.NodesNumber, 0.0);
    Cglobal.resize(data.NodesNumber, Vector(data.NodesNumber, 0.0));

    int globalJ, globalK;

    for (int i = 0; i < data.ElementsNumber; ++i) {
        Element element = elements[i];
        element.setup();
        for (int j = 0; j < 4; ++j) {
            globalJ = element.nodeID[j] - 1;
            for (int k = 0; k < 4; ++k) {
                globalK = element.nodeID[k] - 1;
                Hglobal[globalJ][globalK] += element.H_mat[j][k] + element.Hbc[j][k];
                Cglobal[globalJ][globalK] += element.C_mat[j][k];
            }
            Pglobal[globalJ] += element.P_vec[j];
        }
    }
}

void Solution::printH() {
    for (int i = 0; i < data.NodesNumber; ++i) {
        for (int j = 0; j < data.NodesNumber; ++j) {
            std::cout << std::setw(8) << std::setprecision(5) << Hglobal[i][j] << ", ";
        }
        std::cout << std::endl;
    }
    std::cout<< std::endl;
}

void Solution::printP() {
    for (int i = 0; i < data.NodesNumber; ++i) {
        std::cout << std::setw(8) << std::setprecision(6) << Pglobal[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
}

void Solution::printC() {
    for (int i = 0; i < data.NodesNumber; ++i) {
        for (int j = 0; j < data.NodesNumber; ++j) {
            std::cout << std::setw(8) << std::setprecision(6) << Cglobal[i][j] << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

std::vector<double> Solution::Gauss(const std::vector<std::vector<double>>& A, std::vector<double>& b) {
    int n = A.size();
    std::vector<std::vector<double>> augmentedMatrix = A;
    for (int i = 0; i < n; ++i) {
        augmentedMatrix[i].push_back(b[i]);
    }
    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(augmentedMatrix[k][i]) > std::abs(augmentedMatrix[maxRow][i])) {
                maxRow = k;
            }
        }
        std::swap(augmentedMatrix[i], augmentedMatrix[maxRow]);
        double pivot = augmentedMatrix[i][i];
        for (int j = i; j <= n; ++j) {
            augmentedMatrix[i][j] /= pivot;
        }
        for (int k = i + 1; k < n; ++k) {
            double factor = augmentedMatrix[k][i];
            for (int j = i; j <= n; ++j) {
                augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
            }
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        b[i] = augmentedMatrix[i][n];
        for (int j = i + 1; j < n; ++j) {
            b[i] -= augmentedMatrix[i][j] * b[j];
        }
    }

    return b;
}

void Solution::solve() {
    Global.resize(data.NodesNumber, Vector(data.NodesNumber, 0.0));
    for (int i = 0; i < data.NodesNumber; i++) {
        for (int j = 0; j < data.NodesNumber; j++) {
            if (std::isnan(Hglobal[i][j]) || std::isnan(Cglobal[i][j])) {
                std::cerr << "Error: NaN in global matrices at (" << i << ", " << j << ")" << std::endl;
            }
            Global[i][j] = Hglobal[i][j] + Cglobal[i][j] / data.SimulationStepTime;
        }
    }

    Vector t0(data.NodesNumber, data.InitialTemp);
    Vector t = t0;

    for (double i = data.SimulationStepTime; i <= data.SimulationTime; i += data.SimulationStepTime) {
        for (int row = 0; row < data.NodesNumber; row++) {
            t[row] = 0;
            for (int col = 0; col < data.NodesNumber; col++) {
                t[row] += Cglobal[row][col] / data.SimulationStepTime * t0[col];
            }
            t[row] += Pglobal[row];
        }

        t = Gauss(Global, t);

        std::cout << "Temperature at time " << i << "s:\n";
        std::cout << "MIN: " << *std::min_element(t.begin(), t.end()) << " MAX: " << *std::max_element(t.begin(), t.end()) << std::endl;

        t0 = t;
    }

}

void printNodes(const std::vector<Node>& nodes) {
    std::cout << "Nodes:" << std::endl;
    for (const auto& node : nodes) {
        std::cout << "x: " << node.x << ", y: " << node.y << ", BC: " << node.BC << std::endl;
    }
    std::cout << std::endl;
}

void printElements(const std::vector<Element>& elements) {
    std::cout << "Elements:" << std::endl;
    for (const auto& element : elements) {
        std::cout << "Node IDs: ";
        for (const auto& id : element.nodeID) {
            std::cout << id << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void printGlobalData(const GlobalData& data) {
    std::cout << "Global Data:" << std::endl;
    std::cout << "SimulationTime: " << data.SimulationTime << std::endl;
    std::cout << "SimulationStepTime: " << data.SimulationStepTime << std::endl;
    std::cout << "Conductivity: " << data.Conductivity << std::endl;
    std::cout << "Alfa: " << data.Alfa << std::endl;
    std::cout << "Tot: " << data.Tot << std::endl;
    std::cout << "InitialTemp: " << data.InitialTemp << std::endl;
    std::cout << "Density: " << data.Density << std::endl;
    std::cout << "SpecificHeat: " << data.SpecificHeat << std::endl;
    std::cout << "NodesNumber: " << data.NodesNumber << std::endl;
    std::cout << "ElementsNumber: " << data.ElementsNumber << std::endl;
    std::cout << std::endl;
}

int main() {

    printNodes(nodes);
    printElements(elements);
    printGlobalData(data);
    Solution sol;
    sol.setup();

    std::cout<<"Macierz Hglobal: \n";
    sol.printH();
    std::cout<<"Macierz Cglobal: \n";
    sol.printC();
    std::cout<<"Wektor Pglobal: \n";
    sol.printP();
    sol.solve();

    return 0;
}

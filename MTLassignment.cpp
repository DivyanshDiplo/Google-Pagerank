#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
//2022EE11782
//Divyansh Kumar
using namespace std;

class Graph {
private:
    unordered_map<int, vector<int>> adjList;
    vector<double> weights;

public:
    vector<vector<double>> adjacencyMatrix;
    int n = -1;
    void addEdge(int u, int v) {
        adjList[u].push_back(v);
        n = max(n, u);
        n = max(n, v);
    }

    void genweights() {
        weights.resize(n + 1);
        for (const auto& pair : adjList) {
            int outDegree = pair.second.size();
            weights[pair.first] = (outDegree != 0) ? 1.0 / outDegree : 0.0;
        }
    }

    void generateS() {
        adjacencyMatrix.resize(n + 1, vector<double>(n + 1, 0.0));
        for (const auto& pair : adjList) {
            int u = pair.first;
            for (int v : pair.second) {
                adjacencyMatrix[u][v] = weights[u];
            }
        }

        double danglingWeight = 1.0 / (n + 1);
        for (int i = 0; i <= n; ++i) {
            bool isDangling = true; 
            for (int j = 0; j <= n; ++j) {
                if (adjacencyMatrix[i][j] != 0.0) {
                    isDangling = false; 
                    break;
                }
            }
            if (isDangling) {
                for (int j = 0; j <= n; ++j) {
                    adjacencyMatrix[i][j] = danglingWeight;
                }
            }
        }
    }

    void clear() {
        adjList.clear();
        weights.clear();
        adjacencyMatrix.clear();
    }

};

//matrix operations
vector<vector<double>> matrixAddition(const vector<vector<double>>& matrix1, const vector<vector<double>>& matrix2) {
    vector<vector<double>> result(matrix1.size(), vector<double>(matrix1[0].size()));
    for (size_t i = 0; i < matrix1.size(); ++i) {
        for (size_t j = 0; j < matrix1[0].size(); ++j) {
            result[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
    return result;
}

vector<vector<double>> matrixSubtraction(const vector<vector<double>>& matrix1, const vector<vector<double>>& matrix2) {
    vector<vector<double>> result(matrix1.size(), vector<double>(matrix1[0].size()));
    for (size_t i = 0; i < matrix1.size(); ++i) {
        for (size_t j = 0; j < matrix1[0].size(); ++j) {
            result[i][j] = matrix1[i][j] - matrix2[i][j];
        }
    }
    return result;
}


vector<vector<double>> scalarMatrixMultiplication(const vector<vector<double>>& matrix, double scalar) {
    vector<vector<double>> result(matrix.size(), vector<double>(matrix[0].size()));
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[0].size(); ++j) {
            result[i][j] = matrix[i][j] * scalar;
        }
    }
    return result;
}

vector<vector<double>> matrixMatrixMultiplication(const vector<vector<double>>& matrix1, const vector<vector<double>>& matrix2) {
    vector<vector<double>> result(matrix1.size(), vector<double>(matrix2[0].size(), 0.0));
    for (size_t i = 0; i < matrix1.size(); ++i) {
        for (size_t j = 0; j < matrix2[0].size(); ++j) {
            for (size_t k = 0; k < matrix2.size(); ++k) {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
    return result;
}


// Function to print a matrix
void printMatrix(const vector<vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (double element : row) {
            cout << element << " ";
        }
        cout << endl;
    }
}

double computeDistance(const vector<double>& vec1, const vector<double>& vec2) {
    double distance = 0.0;
    for (int i = 0; i < vec1.size(); ++i) {
        distance += pow(vec1[i] - vec2[i], 2);
    }
    return sqrt(distance);
}

int main() {
    Graph graph;

    return 0;
}
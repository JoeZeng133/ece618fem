#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "windows.h"

#define dirac(i, j) ((i == j)? 1 : 0)
using namespace std;

template<typename T>
class mat {
private:
	vector<T> databuffer;
	int dimx, dimy;

public:
	mat() = default;

	mat(int x, int y) : dimx(x), dimy(y) {
		zeros(x, y);
	}

	void zeros(const int x, const int y) {
		databuffer.resize(x * y, 0);
		dimx = x;
		dimy = y;
	}

	T& operator()(const int x, const int y) {
		return databuffer[(x - 1) * dimy + y - 1];
	}

	T& operator()(const int x) {
		return databuffer[x - 1];
	}

	T* data() {
		return databuffer.data();
	}

	int size(int dim) {
		if (dim == 1)
			return dimx;
		else
			return dimy;
	}
};

class sparse {
public:
	vector<int> i, j;
	vector<double> v;

	void reserve(int n) {
		i.reserve(n);
		j.reserve(n);
		v.reserve(n);
	}

	void push_back(const int i, const int j, const double v) {
		this->i.push_back(i);
		this->j.push_back(j);
		this->v.push_back(v);
	}
};

void assembly(mat<double>& nodeArr, mat<int>& emArr, sparse& A, sparse& B);
void read(mat<double>& nodeArr, mat<int>& emArr);
void output(sparse& A, sparse& B);

int main() {
	mat<double> nodeArr;
	mat<int> emArr;
	sparse A, B;

#if _DEBUG
	char buffer[MAX_PATH];
	GetModuleFileName(NULL, buffer, MAX_PATH);
	cout << "Current Directory is " << buffer << endl;  
#endif // _DEBUG

	read(nodeArr, emArr);
	assembly(nodeArr, emArr, A, B);
	output(A, B);
}

void read(mat<double>& nodeArr, mat<int>& emArr) {
	string filename = "output/mesh.bin";
	fstream file;
	file.open(filename, ios::in | ios::binary);

	if (!file.is_open()) {
		cout << "mesh file open failure" << endl;
		exit(1);
	}

	int numNode, numEm;
	file.read((char*)&numNode, sizeof(int));
	nodeArr.zeros(numNode, 2);
	file.read((char*)nodeArr.data(), numNode * 2 * sizeof(double));

	file.read((char*)&numEm, sizeof(int));
	emArr.zeros(numEm, 3);
	file.read((char*)emArr.data(), numEm * 3 * sizeof(int));

	file.close();
}

void output(sparse& A, sparse& B) {
	string filename = "output/equation.bin";
	fstream file;
	file.open(filename, ios::out | ios::binary);

	if (!file.is_open()) {
		cout << "equation file open failure" << endl;
		exit(1);
	}

	int numEqs = A.i.size();
	file.write((char*)&numEqs, sizeof(int));
	
	file.write((char*)A.i.data(), numEqs * sizeof(int));
	file.write((char*)A.j.data(), numEqs * sizeof(int));
	file.write((char*)A.v.data(), numEqs * sizeof(double));

	file.write((char*)B.i.data(), numEqs * sizeof(int));
	file.write((char*)B.j.data(), numEqs * sizeof(int));
	file.write((char*)B.v.data(), numEqs * sizeof(double));

	file.close();
}

void assembly(mat<double>& nodeArr, mat<int>& emArr, sparse& A, sparse& B) {
	int numNode = nodeArr.size(1);
	int numEm = emArr.size(1);

	mat<int> node(3, 1);
	mat<double> a(numEm, 3), b(numEm, 3), c(numEm, 3);
	mat<double> area(numEm, 1);

	for (int e = 1; e <= numEm; ++e) {
		node(1) = emArr(e, 1);
		node(2) = emArr(e, 2);
		node(3) = emArr(e, 3);

		a(e, 1) = nodeArr(node(2), 1) * nodeArr(node(3), 2) - nodeArr(node(2), 2) * nodeArr(node(3), 1);
		a(e, 2) = nodeArr(node(3), 1) * nodeArr(node(1), 2) - nodeArr(node(3), 2) * nodeArr(node(1), 1);
		a(e, 3) = nodeArr(node(1), 1) * nodeArr(node(2), 2) - nodeArr(node(1), 2) * nodeArr(node(2), 1);

		b(e, 1) = nodeArr(node(2), 2) - nodeArr(node(3), 2);
		b(e, 2) = nodeArr(node(3), 2) - nodeArr(node(1), 2);
		b(e, 3) = nodeArr(node(1), 2) - nodeArr(node(2), 2);

		c(e, 1) = nodeArr(node(3), 1) - nodeArr(node(2), 1);
		c(e, 2) = nodeArr(node(1), 1) - nodeArr(node(3), 1);
		c(e, 3) = nodeArr(node(2), 1) - nodeArr(node(1), 1);

		area(e) = 0.5 * (b(e, 1) * c(e, 2) - b(e, 2) * c(e, 1));

		A.reserve(numEm * 9);
		B.reserve(numEm * 9);

		for (int i = 1; i <= 3; ++i)
			for (int j = i; j <= 3; ++j) {
				double tmp = 0.25 / area(e) * (b(e, i) * b(e, j) + c(e, i) * c(e, j));
				A.push_back(node(i), node(j), tmp);
				if (i != j) A.push_back(node(j), node(i), tmp);

				tmp = area(e) / 12 * (1 + dirac(i, j));
				B.push_back(node(i), node(j), tmp);
				if (i != j) B.push_back(node(j), node(i), tmp);
			}
	}
}




#pragma once
#include "Vector.h"

template <class T>
class MatrixT {
	VectorT<Vector> data; 
	int n;
	int m;
public:
	int getN() const { return n; };
	int getM() const { return m; };
	MatrixT(int n = 0, int m = 0) : n(n), m(m) {
		data = VectorT<Vector>(n);
		for (int i = 0; i < n; i++) this->data[i] = Vector(m);
	}
	MatrixT(const Vector & v) :n(v.getSize()), m(1){
		data = VectorT<Vector>(n);
		for (int i = 0; i < n; i++) this->data[i] = Vector(1);
		for (int i = 0; i < n; i++) this->data[i][0] = v[i];
	}
	double & operator()(int ind1, int ind2) {
		return this->data[ind1][ind2];
	}
	double & operator()(int ind1, int ind2) const{
		return this->data[ind1][ind2];
	}
	Vector & operator[](int index) const {
		return data[index];
	}
	friend std::ostream& operator<<(std::ostream& ostr, const MatrixT & M) {
		for (int i = 0; i < M.getN(); i++) ostr << M[i];
		return ostr;
	}
	Vector  operator*(const Vector & a){
		if (this->n == a.getSize()) {
			VectorT<T> tmp(a), res(n);
			for (int i = 0; i < n; i++) res[i]= (int)(data[i] * tmp);
			return res;
		}
		else throw "Inccorect operator*=\n";
	}
	/*void Transposition() {
		VectorT<Vector> tmp=data;
		this->n = tmp.m;
		this->m = tmp.n;
		data = VectorT<Vector>(n);
		for (int i = 0; i < n; i++) data[i] = Vector(m);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				data[i][j] = tmp[j][i];
			}
		}
	}*/
};
using Matrix = MatrixT<double>;



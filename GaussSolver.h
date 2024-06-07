#pragma once
#include "Matrix.h"
#include <vector>


class GaussSolver {
public:
	bool comp(double& b) const {
		double epsilone=0.;
		return abs(b)==epsilone;
	}
	void swap(const Matrix & A, int i, int j) {
		for (int t = 0; t < A.getM(); t++) {
			double  tmp = A[i][t];
			A[i][t] = A[j][t];
			A[j][t] = tmp;
		}
	}
	std::vector<Vector> solve(const Matrix& A, const Vector& v){
		Matrix M(A.getN(), A.getM() + 1);
		for (int i = 0; i < M.getN(); i++) {
			M[i][M.getM() - 1] = v[i];
			for (int j = 0; j < A.getM(); j++) M[i][j] = A[i][j];
		}
		//std::cout << M << std::endl;
		int n = M.getN(), m = M.getM();
		int rank = 0;
		VectorT<int>  indexVed(n);
		for (int i = 0; i < n; i++) indexVed[i] = 100000;

		for (int i = 0; i < n; i++) {
			double Ved = -100000000;
			int indexIV = -1;
			for (int j = 0; j < m-1; j++) {
				for (int t = i; t < n; t++) {
					if ((M[t][j] != 0) && ((abs(M[t][j]) > abs(Ved))||(Ved== -100000000))) {
						Ved = M[t][j];
						indexIV = t;
						indexVed[i] = j;
					}
				}
				if (Ved != -100000000) {
					rank++;
					swap(M, indexIV, i);
					break;
				}
			}
			//std::cout << Ved << std::endl;
			//std::cout << M << std::endl;
			if(Ved!= -100000000){
				for (int j = 0; j < m; j++) M[i][j] /= Ved;
				Ved = 1;
				for (int t = 0; t < n; t++) {
					if (t != i) {
						double tmp = (double)M[t][indexVed[i]] / Ved;
						for (int j = 0; j < m; j++) {
							M[t][j] -= (double)M[i][j] * tmp;

						}
					}
				}
				//std::cout << M << std::endl;
			}
		}
		for (int i = 1; i < n; i++) {
			for (int j = 1; j < n; j++) {
				if (indexVed[j] < indexVed[j - 1]) {
					swap(M, j, j - 1);
				}
			}
		}
		bool Flag = true;
		for (int i = 0; i < n; i++) {
			int j = 0;
			bool flag = false;
			while (j < m - 1) {
				if (!comp(M[i][j])) {
					flag = true;
					break;
				}
				j++;
			}
			if ((flag == false)&&(!comp(M[i][m-1]))) Flag = false;
		}
		//std::cout << M << std::endl;
		//std::cout << indexVed << rank << std::endl;
		if (Flag == true) {
			std::vector<Vector> ans(m - rank);
			for (int i = 0; i < m - rank; i++) ans[i] = Vector(m - 1);
			int p = 1;
			for (int i = 0; i < m - 1; i++) {
				bool flag = true;
				for (int j = 0; j < n; j++) { // проверка: текущий столбец отвечает за базисную переменную или нет
					if (indexVed[j] == i) {
						flag = false;
						ans[0][j] = M[j][m - 1];
						break;
					}
				}
				if (flag == true) {
					int a = 0;
					if (m-1< n) a = m-1;
					else a = n;
					for (int j = 0; j < a; j++) ans[p][j] = -M[j][i]; //
					ans[p][i] = 1;
					p++;
				}
			}
			return ans;
		}
		std::vector<Vector> ans;
		return ans;
	}
};

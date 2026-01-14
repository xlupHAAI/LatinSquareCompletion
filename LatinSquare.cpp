#include "LatinSquare.h"

#include <random>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <numeric>

using namespace std;

class Domain {
private:
	vector<int> idx;
public:
	int size;
	vector<int> elems;

	Domain(const vector<int>& v) :size(v.size()), elems(v), idx(vector<int>(MAXN)) {
		for (int i = 0; i < size; ++i)	idx[v[i]] = i;
	}
	Domain():size(0) {}

	void erase(int val) {
		if (idx[val] < 0)	return;
		idx[elems[--size]] = idx[val];
		elems[idx[val]] = elems[size];
		idx[val] = NONE;
	}

	Domain& operator=(const Domain& y) {
		size = y.size;
		idx = y.idx;
		elems = y.elems;
		return *this;
	}

	inline bool has(int val) { return idx[val] > NONE; }
};

class bi_set {
private:
	vector<int> idx;
public:
	int size;
	vector<int> elems;

	bi_set() :idx(vector<int>(MAXN, -1)), elems(vector<int>(MAXN)), size(0) {  }

	void insert(int val) {
		if (idx[val] == NONE) {
			elems[size] = val;
			idx[val] = size++;
		}
	}

	void erase(int val) {
		if (idx[val] > NONE) {
			idx[elems[--size]] = idx[val];
			elems[idx[val]] = elems[size];
			idx[val] = NONE;
		}
	}

	bi_set& operator=(const bi_set& y) {
		size = y.size;
		idx = y.idx;
		elems = y.elems;
		return *this;
	}
	inline bool has(int val) { return idx[val] > NONE; }
};

namespace szx {

class Solver {
	// random number generator.
	mt19937 pseudoRandNumGen;
	void initRand(int seed) { pseudoRandNumGen = mt19937(seed); srand(seed); }
	int fastRand(int ub) { return pseudoRandNumGen() % ub; }

private:
	struct Move {
		Num row;
		Num u, v;
		int benifit;
		Move(Num _r, Num _u, Num _v, int _f) :row(_r), u(_u), v(_v), benifit(_f) {  }
		Move() :benifit(INF) {  }
	};

	Num n;					// dimension of latin-square
	uint nbiter;			// iteration count. reset when each time restart happens
	int accu = 0, rt = 0;	// params in adaptive restart procedure

	vector<vector<bool>> is_fixed;	// if fixed after the constraints propagate
	vector<NumArr> not_fixed;		// cells whose color is not fixed

	Table S;			// current solution of the input square
	Table nbcc;			// nbcc[col][color]: cells number of each color on each column 
	Table nbcc_fixed;	// fixed cells number of each color on each column

	Table S_best;		// best solution found so far
	Table nbcc_best;	// nbcc under S_best

	int conflictEdgeNum;	// total amount of conflict edges at current
	int conflictNumBest;	// conflictEdgeNum of S_best

	struct tabu_item {
		Num i, u, Cu;
		unsigned int tenure;
		tabu_item(Num row, Num a, Num Ca, unsigned int t) :i(row), u(a), Cu(Ca), tenure(t) {  }
		bool operator<(const tabu_item& y)const {
			return tenure > y.tenure;
		}
	};
	priority_queue<tabu_item> tabu_queue;
	vector<vector<bool>> tabu_list[MAXN];

	vector<bi_set> conflict_cells;
	vector<bi_set> conflict_cells_best;
	
	void initLSC(LatinSquare&);
	inline void initParams(LatinSquare&);
	inline void initMemory();
	inline void initfS();
	inline void fixCell(vector<vector<Domain>>&, int, int, int);

	void tabuExpired();
	Move findMove();
	inline void makeMove(Move&);
	void makeMoveImpl(Num, Num, Num, Num, Num);
	void AdaptiveRestart();

	inline void resetTabuInfo();
	inline void replaceBestSolution();
	inline void restartFromTheBest();

	// debug
	bool checkErrSolution();

public:
	void solve(Table& output, LatinSquare& input, function<long long()> restMilliSec, int seed) {
		
		initRand(seed);

		initLSC(input);

		if (conflictEdgeNum == 0) {

			//if (checkErrSolution()) {
			//	cout << "\n\nerror answer.";
			//	exit(-1);
			//}
			//printf("duration: %lld.\n\n", 100000000 - restMilliSec());
			output = S; 
			return;
		}

		while (restMilliSec() > 0) {

			tabuExpired();

			Move best_move = findMove();

			makeMove(best_move);

			if (conflictEdgeNum == 0) {
				
				//if (checkErrSolution()) {
				//	cout << "\n\nerror answer.";
				//	exit(-1);
				//}
				//printf("\nduration: %lld.\n", 100000 - restMilliSec());
				output = S;
				return;
			}
			else if (conflictEdgeNum <= conflictNumBest) {

				replaceBestSolution();
			}
			else {

				AdaptiveRestart();
			}
		}

		//cout<<conflictNumBest;
		output = S_best;
	}
};

void Solver::tabuExpired() {
	while (!tabu_queue.empty()) {
		const tabu_item &ti = tabu_queue.top();
		if (ti.tenure < nbiter) {
			tabu_list[ti.i][ti.u][ti.Cu] = false;
			tabu_queue.pop();
		}
		else break;
	}
}

inline void Solver::resetTabuInfo() {
	nbiter = 1;
	priority_queue<tabu_item>().swap(tabu_queue);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			fill(tabu_list[i][j].begin(), tabu_list[i][j].end(), false);
}

inline void Solver::initParams(LatinSquare& ls) {
	n = ls.n;
	accu = 0;
	rt = rt0;
	nbiter = 1;
}

inline void Solver::initMemory() {
	is_fixed = vector<vector<bool>>(n, vector<bool>(n));
	not_fixed = vector<NumArr>(n);

	S = Table(n, NumArr(n));
	//nbcc = Table(n, NumArr(n));
	nbcc_fixed = Table(n, NumArr(n));

	for (int i = 0; i < n; ++i)	tabu_list[i] = vector<vector<bool>>(n, vector<bool>(n));

	conflict_cells = vector<bi_set>(n, bi_set());
}

inline void Solver::initfS() {
	conflictEdgeNum = 0;
	for (Num i = 0; i < n; i++)
		for (Num j = 0; j < n; j++) {
			conflictEdgeNum += (nbcc[i][j] - 1) * nbcc[i][j] >> 1;
			if (!is_fixed[i][j] && nbcc[j][S[i][j]] > 1)
				conflict_cells[i].insert(j);
		}
}

inline void Solver::fixCell(vector<vector<Domain>>& D, int i, int j, int c) {
	S[i][j] = c;
	nbcc_fixed[j][c]++;
	is_fixed[i][j] = true;

	for (int k = 0; k < n; ++k) {
		if (!is_fixed[k][j]) {
			D[k][j].erase(c);
		}
		if (!is_fixed[i][k]) {
			D[i][k].erase(c);
		}
	}
}

void Solver::initLSC(LatinSquare& ls) {

	initParams(ls);

	initMemory();

	vector<int> color_domain(n);		
	for (int i = 1; i < n; i++)
		color_domain[i] = i;

	vector<vector<Domain>> D(n, vector<Domain>(n, Domain(color_domain)));
	for (Assignment& a : ls.fixedNums) {
		fixCell(D, a.row, a.col, a.num);
	}

_propagate:
	
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			if (!is_fixed[i][j] && D[i][j].size == 1) {
				fixCell(D, i, j, D[i][j].elems[0]);
				goto _propagate;
			}

	vector<int> nbcolor(n);
	for (int i = 0; i < n; ++i) {
		fill(nbcolor.begin(), nbcolor.end(), 0);
		for (int j = 0; j < n; ++j)
			if (!is_fixed[i][j])
				for (int ei = 0; ei < D[i][j].size; ++ei)
					nbcolor[D[i][j].elems[ei]]++;

		for (int c = 0; c < n; c++)
			if (nbcolor[c] == 1) {
				int k;
				for (k = 0; k < n; ++k)
					if (!is_fixed[i][k] && D[i][k].has(c))
						break;
				fixCell(D, i, k, c);
				goto _propagate;
			}
	}

	for (int j = 0; j < n; ++j) {
		fill(nbcolor.begin(), nbcolor.end(), 0);
		for (int i = 0; i < n; ++i)
			if (!is_fixed[i][j])
				for (int ei = 0; ei < D[i][j].size; ++ei)
					nbcolor[D[i][j].elems[ei]]++;

		for (int c = 0; c < n; c++)
			if (nbcolor[c] == 1) {
				int k;
				for (k = 0; k < n; ++k)
					if (!is_fixed[k][j] && D[k][j].has(c))
						break;
				fixCell(D, k, j, c);
				goto _propagate;
			}
	}

	nbcc = nbcc_fixed;
		
	for (Num i = 0; i < n; ++i) {
		Domain row_domain(color_domain);
		for (Num j = 0; j < n; ++j)
			if (!is_fixed[i][j])
				not_fixed[i].push_back(j);
			else row_domain.erase(S[i][j]);

		random_shuffle(row_domain.elems.begin(), row_domain.elems.begin() + row_domain.size);

		for (Num j = 0; j < row_domain.size; ++j) {
			Num u = not_fixed[i][j];
			S[i][u] = row_domain.elems[j];
			nbcc[u][S[i][u]]++;
		}
	}

	initfS();

	replaceBestSolution();
}

void Solver::AdaptiveRestart() {
	if (conflictEdgeNum - conflictNumBest > rt) {
		restartFromTheBest();
		resetTabuInfo();
		if (rt < rtub && ++accu == accuub) {
			accu = 0;
			++rt;
		}
	}
}

inline void Solver::replaceBestSolution() {
	conflictNumBest = conflictEdgeNum;
	S_best = S;
	nbcc_best = nbcc;
	conflict_cells_best = conflict_cells;
}

inline void Solver::restartFromTheBest() {
	conflictEdgeNum = conflictNumBest;
	S = S_best;
	nbcc = nbcc_best;
	conflict_cells = conflict_cells_best;
}

void Solver::makeMoveImpl(Num i, Num u, Num v, Num Cu, Num Cv) {
	int tenure = nbiter++ + alpha * conflictEdgeNum + fastRand(10);
	
	if (nbcc[u][Cu] > 1) {
		tabu_queue.push(tabu_item(i, u, Cu, tenure));
		tabu_list[i][u][Cu] = true;
	}
	if (nbcc[v][Cv] > 1) {
		tabu_queue.push(tabu_item(i, v, Cv, tenure));
		tabu_list[i][v][Cv] = true;
	}

	nbcc[u][Cu]--;	nbcc[v][Cv]--;
	nbcc[u][Cv]++;	nbcc[v][Cu]++;	

	swap(S[i][u], S[i][v]);

	if (nbcc[u][Cu] == 1 && nbcc_fixed[u][Cu] == 0) {
		for (int k = 0; k < n; ++k)
			if (S[k][u] == Cu) {
				conflict_cells[k].erase(u);
				break;
			}
	}
	if (nbcc[u][Cv] > 1)	conflict_cells[i].insert(u);
	else conflict_cells[i].erase(u);
	if (nbcc[u][Cv] == 2 && nbcc_fixed[u][Cv] == 0) {
		for (int k = 0; k < n; ++k)
			if (S[k][u] == Cv && k != i) {
				conflict_cells[k].insert(u);
				break;
			}
	}

	if (nbcc[v][Cv] == 1 && nbcc_fixed[v][Cv] == 0) {
		for (int k = 0; k < n; ++k)
			if (S[k][v] == Cv) {
				conflict_cells[k].erase(v);
				break;
			}
	}
	if (nbcc[v][Cu] > 1)	conflict_cells[i].insert(v);
	else conflict_cells[i].erase(v);
	if (nbcc[v][Cu] == 2 && nbcc_fixed[v][Cu] == 0) {
		for (int k = 0; k < n; ++k)
			if (S[k][v] == Cu && k != i) {
				conflict_cells[k].insert(v);
				break;
			}
	}
}

inline void Solver::makeMove(Move& m) {
	conflictEdgeNum += m.benifit;
	makeMoveImpl(m.row, m.u, m.v, S[m.row][m.u], S[m.row][m.v]);
}


Solver::Move Solver::findMove() {
	Move best_move_T, best_move_NT;
	int countT = 1, countNT = 1;
	int delta_fixed_T = INF, delta_fixed_NT = INF;

	for (Num i = 0; i < n; ++i) {
		Domain not_fixed_cells(not_fixed[i]);

		for (int idx = 0; idx < conflict_cells[i].size; ++idx) {
			Num u = conflict_cells[i].elems[idx];
			Num Cu = S[i][u];
			not_fixed_cells.erase(u);
			for (int vidx = 0; vidx < not_fixed_cells.size; ++vidx) {
				Num v = not_fixed_cells.elems[vidx];
				Num Cv = S[i][v];

				int benifit = nbcc[u][Cv] + nbcc[v][Cu] - nbcc[u][Cu] - nbcc[v][Cv];
				
				if (!tabu_list[i][u][Cv] && !tabu_list[i][v][Cu]) {

					if (benifit == best_move_NT.benifit) {

						int delta = nbcc_fixed[u][Cv] + nbcc_fixed[v][Cu] - nbcc_fixed[u][Cu] - nbcc_fixed[v][Cv];

						if (delta < delta_fixed_NT) {
							countNT = 1;
							best_move_NT = Move(i, u, v, benifit);
							delta_fixed_NT = delta;
						}
						else if (delta == delta_fixed_NT && fastRand(++countNT) == 0) {
							best_move_NT = Move(i, u, v, benifit);
						}

					}
					else if (benifit < best_move_NT.benifit) {
						countNT = 1;
						best_move_NT = Move(i, u, v, benifit);
						delta_fixed_NT = nbcc_fixed[u][Cv] + nbcc_fixed[v][Cu] - nbcc_fixed[u][Cu] - nbcc_fixed[v][Cv];
					}				
				}

				else {
					if (benifit == best_move_T.benifit) {

						int delta = nbcc_fixed[u][Cv] + nbcc_fixed[v][Cu] - nbcc_fixed[u][Cu] - nbcc_fixed[v][Cv];

						if (delta < delta_fixed_T) {
							countT = 1;
							best_move_T = Move(i, u, v, benifit);
							delta_fixed_T = delta;
						}
						else if (delta == delta_fixed_T && fastRand(++countT) == 0) {
							best_move_T = Move(i, u, v, benifit);
						}
					}
					else if (benifit < best_move_T.benifit) {
						countT = 1;
						best_move_T = Move(i, u, v, benifit);
						delta_fixed_T = nbcc_fixed[u][Cv] + nbcc_fixed[v][Cu] - nbcc_fixed[u][Cu] - nbcc_fixed[v][Cv];
					}
				}

			}
		}
	}

	best_move_T.benifit += 2;
	best_move_NT.benifit += 2;

	return (best_move_T.benifit < conflictNumBest - conflictEdgeNum)
		&& (best_move_T.benifit < best_move_NT.benifit
			|| (best_move_T.benifit == best_move_NT.benifit && delta_fixed_T < delta_fixed_NT))
		? best_move_T
		: best_move_NT;
}





bool Solver::checkErrSolution() {
	NumArr row_exists(n), col_exists(n);
	for (Num i = 0; i < n; ++i) {
		fill(row_exists.begin(), row_exists.end(), 0);
		fill(col_exists.begin(), col_exists.end(), 0);
		for (Num j = 0; j < n; ++j) {
			if (row_exists[S[i][j]]++ > 0)	return true;
			if (col_exists[S[j][i]]++ > 0)	return true;
		}
	}
	return false;
} 








// solver.
void solveLatinSquare(Table& output, LatinSquare& input, function<long long()> restMilliSec, int seed) {
	Solver().solve(output, input, restMilliSec, seed);
}

}

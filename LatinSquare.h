////////////////////////////////
/// usage : 1.	SDK for graph coloring solver.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_SZX_NPBENCHMARK_LATIN_SQUARE_H
#define CN_HUST_SZX_NPBENCHMARK_LATIN_SQUARE_H


#include <array>
#include <vector>
#include <functional>
#include <unordered_set>

#define NONE -1
#define INF 0x7ffff000
#define MAXN 70

namespace szx {

using Num = int;
struct Assignment {
	Num row;
	Num col;
	Num num;
};

struct LatinSquare {
	Num n;
	std::vector<Assignment> fixedNums; // fixed numbers.
};

using uint = unsigned int;
using NumArr = std::vector<Num>;
using Table = std::vector<std::vector<Num>>; // a 2D array of numbers.

constexpr double alpha = 0.4;
constexpr int rt0 = 10;
constexpr int rtub = 15;
constexpr int accuub = 1000;








void solveLatinSquare(Table& output, LatinSquare& input, std::function<long long()> restMilliSec, int seed);

}


#endif // CN_HUST_SZX_NPBENCHMARK_LATIN_SQUARE_H

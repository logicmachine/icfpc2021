#ifndef __RANDOM_HPP__
#define __RANDOM_HPP__
using ll = long long int;

// [lb, ub] の閉区間内の値をランダムに返す構造体
// #include <random> しよう

struct Rand {
public:
    Rand() = default;
    Rand(std::mt19937::result_type seed) : eng(seed) {}
    int NextInt(int lb, int ub) {
        return std::uniform_int_distribution<int>{lb, ub}(eng);
    }
    ll NextLong(ll lb, ll ub) {
        return std::uniform_int_distribution<ll>{lb, ub}(eng);
    }
    double NextDouble(double lb = 0.0, double ub = 1.0) {
        return std::uniform_real_distribution<double>{lb, ub}(eng);
    }
private:
    std::mt19937 eng{std::random_device{}()};
};

#endif // !__RANDOM_HPP__

/* 
// example.
int main() {
    Rand rnd(114514);
    int l, r; scanf("%d%d", &l, &r);
    printf("l = %d, r = %d, value = %d\n", l, r, rnd.NextInt(l, r));

    double L, R; scanf("%lf%lf", &L, &R);
    printf("L = %.12f, R = %.12f, value = %.12f\n", L, R, rnd.NextDouble(L, R));

    long long int a, b; scanf("%lld%lld", &a, &b);
    printf("a = %lld, b = %lld, value = %lld\n", a, b, rnd.NextLong(a, b));
    return 0;
}
*/

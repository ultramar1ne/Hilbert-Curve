#include <vector>
#include <array>
#include <algorithm>
#include <limits>
#include <utility>
#include <iostream>
#include <iomanip>

#include <parlay/parallel.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>

const int dim_d = 5;

template<typename T, size_t d>
class point_d {
public:
    std::array<T, d> coordinates_;
    point_d(std::initializer_list<T> v) {
        std::copy(v.begin(), v.end(), coordinates_.begin());
    }
    point_d()= default;
    explicit point_d(const std::array<T, d>& arr) : coordinates_{arr} {}
    auto operator[](size_t dimension) { return &coordinates_[dimension]; }  // can be edited

    //const T operator[](size_t dimension) const { return coordinates_[dimension]; }
};

using points = std::vector<point_d<int,dim_d>>;
using point = point_d<int,dim_d>;

class HilbertCurve{
public:
    HilbertCurve(int p, int dim);
    int p;
    int dim;
    std::vector<int> distance;
    auto Points2Dis(points& ps);
    auto Point2Dis( point& p);

//private:
    long transpose2HInt(point& p);
    std::string  binRepr(int num, int width);
};

HilbertCurve::HilbertCurve(int p, int dim) {
    this->dim = dim;
    this->p = p;
}

auto HilbertCurve::Point2Dis( point& p) {
    auto m = 1<<(this->p -1);
    auto q = m;
    while(q>1){
        auto pp = q-1;
        for(int i = 0; i< this->dim; i++){
            if(*p[i]&q){
                *p[0] ^= pp;
            }else{
                auto t = (*p[0] ^ *p[i]) & pp;
                *p[0] ^= t;
                *p[i] ^= t;
            }
        }
        q >>=1;
    }

    //gray code
    for(int i = 1; i< this->dim; i++){
        *p[i] ^= *p[i-1];
    }
    auto t = 0;
    q = m;
    while(q>1){
        if(*p[this->dim-1] &q){
            t ^= q-1;
        }
        q>>=1;
    }
    for(auto i = 0; i< this->dim; i++){
        *p[i] ^= t;
    }
    return this->transpose2HInt(p);
}

auto HilbertCurve::Points2Dis(points &ps) {
    auto n = ps.size();
    long* res = new long[n];
    auto f = [&](size_t i){
        res[i] = this->Point2Dis(ps[i]);
    };
    parlay::parallel_for(0,n,f);
    return res;
}

long HilbertCurve::transpose2HInt(point &p) {
    auto x = std::vector<std::string>(dim_d);
    for(int i = 0; i<dim_d; i++){
        x[i] = binRepr(*p[i],this->p);
    }
    std::string h;
    for(int j = 0; j<this->p; j++){
        for(int i = 0; i<dim_d; i++){
            h+= x[i][j];
        }
    }
    return stol(h, 0, 2);;
}

std::string HilbertCurve::binRepr(int num, const int width) {
    std::string binnum = std::bitset<25>(num).to_string(); //todo: hardcoded
    return  binnum;
}

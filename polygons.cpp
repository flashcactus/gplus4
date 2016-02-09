#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

/*

^ Y
|
|
|      X
*----->

*/

template <typename T> struct MyVec2 {
    T x;
    T y;

    MyVec2(T x=0, T y=0) {
        this->x = x;
        this->y = y;
    }

    template<typename H> friend MyVec2<H> operator+ (const MyVec2<H> v1, const MyVec2<H> v2);
    template<typename H> friend MyVec2<H> operator- (const MyVec2<H> v1, const MyVec2<H> v2);
    MyVec2 operator- () const {return MyVec2(-x, -y);}
    
    MyVec2 operator/ (const T n) const { return MyVec2(x/n, y/n); }
    template<typename H> friend H operator* (const MyVec2<H> v1, const MyVec2<H> v2);
    template<typename H> friend MyVec2 operator* (const H n, const MyVec2<H> v);
    template<typename H> friend MyVec2 operator* (const MyVec2<H> v, const H n);

    double arg() const { return atan2(y, x); }

    double abs() const { return sqrt(x*x + y*y); }
    double sqabs() const { return x*x + y*y; }

    MyVec2 normalized() const { return (*this)/this->abs(); }

    MyVec2 rotcw() const { return MyVec2(y, -x); }
    MyVec2 rotccw() const { return MyVec2(-y, x); }
};

template<typename H> 
MyVec2<H> operator+ (const MyVec2<H> v1, const MyVec2<H> v2) { 
    return MyVec2<H>(v1.x+v2.x, v1.y+v2.y); 
}

template<typename H> 
MyVec2<H> operator- (const MyVec2<H> v1, const MyVec2<H> v2) { 
    return MyVec2<H>(v1.x-v2.x, v1.y-v2.y); 
}

template<typename H> 
H operator* (const MyVec2<H> v1, const MyVec2<H> v2) { 
    return v1.x*v2.x + v1.y*v2.y; 
}

template<typename H> 
MyVec2<H> operator* (const H n, const MyVec2<H> v) {
    return MyVec2<H>(v.x*n, v.y*n);
}

template<typename H> 
MyVec2<H> operator* (const MyVec2<H> v, const H n) {
    return MyVec2<H>(v.x*n, v.y*n);
}


typedef MyVec2<float> fvec ;
typedef MyVec2<double> dvec ;
typedef fvec pt;

//-----------------------------------------------

template <typename T> struct Line {
    private:
        MyVec2<T> norm;
    public:
        T h;
    
    T setnrm(const MyVec2<T> n) {
        norm = n.normalized(); 
    }
    T setvec(const MyVec2<T> v) {
        norm = v.normalized().rotcw();
    }

    Line(const MyVec2<T> n, const T height) { 
        setnrm(n);
        h = height;
    }
    Line(const MyVec2<T> p1, const MyVec2<T> p2) {
        setvec(p2-p1);
        h = p1*norm;
    }
    
    MyVec2<T> getnrm() { return norm; }
    MyVec2<T> getvec() { return norm.rotccw(); } 

    void fliptoward(MyVec2<T> pnt) {
        if ( (pnt*norm - h) < 0 ) {
            h = -h;
            norm = -norm;
        }
    }

    T sgdist(MyVec2<T> pnt) const {
       return pnt*norm - h;
    }
    
    MyVec2<T> operator& (const Line<T> oth) { //find intersection
        MyVec2<T> r = (oth.norm.y*h - norm.y*oth.h , - oth.norm.x*h + norm.x*oth.h );
        T c = norm.x * oth.norm.y - norm.y * oth.norm.x;
        return r/c;
    }
};

//-------------------------------------------------

template <typename T> struct ConvPolygon {
    T EPS = 0.001;
    vector<Line<T> > lines;
    
    ConvPolygon(){};
   
    ConvPolygon(const MyVec2<T>* pts, int size) {
        ConvPolygon(vector<MyVec2<T>>(pts, pts+size));
    }

    ConvPolygon(const vector<MyVec2<T>> pts) {
        vector<MyVec2<T>> points = pts;
        sort(
            points.begin(), points.end(), 
            [] (const pt &p1, const pt &p2) { return (p1.arg() < p2.arg()); }
        );
        mpv(points);
        MyVec2<T> CoM = *(points.begin());
        for(auto i=points.begin()+1; i<points.end(); ++i) {
            CoM = CoM + *i;
            lines.push_back(Line<T>(*(i-1), *i));
        }
        lines.push_back(Line<T>(*(points.end()-1), *(points.begin())));
        mpl(lines);
        CoM = CoM/points.size();
        pv(CoM);
        for(auto i=lines.begin()+1; i<lines.end(); ++i) {
            i->fliptoward(CoM);
        }
        mpl(lines);
        cout<<lines.size()<<endl;
    }
   
    
    bool is_in(MyVec2<T> pnt, T eps = 0) {
        for(auto k=lines.begin(); k!= lines.end(); ++k) {
            if( k->sgdist(pnt) < -eps ) {return false;}
        }
        return true;
    }

    vector<MyVec2<T>> getpts() const {
        vector<MyVec2<T>> points;
        MyVec2<T> X;
        bool ok;
        for ( auto i=lines.begin(); i!= lines.end(); ++i ) {
           for (auto j=i+1; j!= lines.end(); ++j) {
                X = *i & *j;
                if(is_in(X,EPS)) { points.push_back(X); }
            }
        }
        sort(
            points.begin(), points.end(), 
            [] (const pt &p1, const pt &p2) { return (p1.arg() < p2.arg()); }
        );
    }
};

template <typename T> void pv(MyVec2<T> v) {
    cout<<"( "<<v.x<<", "<<v.y<<" )"<<endl;
}
template <typename T> void mpv(vector< MyVec2<T> > vv) {
    for(auto i = vv.begin(); i<vv.end(); ++i) { pv(*i); }
}

template <typename T> void pl(Line<T> l) {
    cout<<l.h<<' ';
    pv(l.getnrm());
}
template <typename T> void mpl(vector< Line<T> > vl) {
    for(auto i = vl.begin(); i<vl.end(); ++i) { pl(*i); }
}


int main() {
    fvec x(1,0);
    fvec y(0,1);
    
    pt points[] = {pt(3, -1), pt(-3, -1), pt(0, 3)};
    ConvPolygon<float> poly(points, 3);
    cout<<poly.lines.size();
    mpl(poly.lines);
}


#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <GLFW/glfw3.h>
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
    
    void setnrm(const MyVec2<T> n) {
        norm = n.normalized(); 
    }
    void setvec(const MyVec2<T> v) {
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
        if ( sgdist(pnt) <= 0 ) {
            h = -h;
            norm = -norm;
        }
    }

    T sgdist(MyVec2<T> pnt) const {
       return pnt*norm - h;
    }
    
    MyVec2<T> operator& (const Line<T> oth) const { //find intersection
        MyVec2<T> r(oth.norm.y*h - norm.y*oth.h, - oth.norm.x*h + norm.x*oth.h );
        T c = norm.x * oth.norm.y - norm.y * oth.norm.x;
        return r/c;
    }
};

//-------------------------------------------------

template <typename T> struct ConvPolygon {
    public:
        T EPS = 0.001;
        vector<Line<T> > lines;
        
        ConvPolygon(){};
       
        ConvPolygon(const MyVec2<T>* pts, int size) {
            fill(pts, size);
        }

        ConvPolygon(const vector<MyVec2<T>> pts) {
            fill(pts);
        }
        
        vector<Line<T>> fill(const MyVec2<T>* pts, int size) {
            return fill(vector<MyVec2<T>>(pts, pts+size));
        }

        vector<Line<T>> fill(const vector<MyVec2<T>> pts) {
            vector<MyVec2<T>> points = pts;
            sort(
                points.begin(), points.end(), 
                [] (const pt &p1, const pt &p2) { return (p1.arg() < p2.arg()); }
            );
            MyVec2<T> CoM = *(points.begin());
            for(auto i=points.begin()+1; i<points.end(); ++i) {
                CoM = CoM + *i;
                lines.push_back(Line<T>(*(i-1), *i));
            }
            lines.push_back(Line<T>(*(points.end()-1), *(points.begin())));
            CoM = CoM/points.size();
            pv(CoM);
            for(auto i=lines.begin(); i<lines.end(); ++i) {
                i->fliptoward(CoM);
            }
            return lines;
        }

        ConvPolygon(const ConvPolygon& p) {
            EPS = p.EPS;
            lines = p.lines;
        }
        
        bool is_in(MyVec2<T> pnt, T eps = 0) const {
            for(auto k=lines.begin(); k!= lines.end(); ++k) {
                if( !(k->sgdist(pnt) >= -eps) ) { return false;}
            }
            return true;
        }

        vector<MyVec2<T>> getpts() const {
            vector<MyVec2<T>> points;
            MyVec2<T> X;
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
            return points;
        }

        ConvPolygon<T> intersect(ConvPolygon<T> &other) const {
            ConvPolygon<T> poly;
            poly.lines = lines;
            poly.lines.insert( poly.lines.end(), other.lines.begin(), other.lines.end() );
            return poly;
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


void error_callback(int error, const char* description)
{
        fputs(description, stderr);
}

const int steps = 500;
int stepsleft;
static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
        if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
            glfwSetWindowShouldClose(window, GL_TRUE);
        
        if (key == GLFW_KEY_SPACE && action == GLFW_PRESS){
            stepsleft = steps;
        }
}


int main() {
    fvec x(1,0);
    fvec y(0,1);
    
    pt points[4] = {pt(3, -1), pt(-3, -1), pt(1, 3), pt(-1,3)};
    pt antipoints[4];
    for(auto i =0; i<4; ++i) { antipoints[i] = -points[i]; }

    ConvPolygon<float> poly(points, 4);
    vector<pt> pts1 = poly.getpts();
    ConvPolygon<float> poly2(antipoints, 4);
    vector<pt> pts2 = poly2.getpts();
    
    //get values for morphing
    vector<float> h1start, h1step, h2start, h2step;
    float h;
    
    //fit p1 over p2
    for (auto l = poly.lines.begin(); l!=poly.lines.end(); ++l){
        h = l->h;
        for(auto p = pts2.begin(); p!=pts2.end(); ++p) {
            if(*p * l->getnrm() < h) {h = *p * l->getnrm();}
        }
        h1start.push_back(l->h);//start small
        h1step.push_back((h - l->h)/steps);//grow big
    }

    //fit p2 over p1 now
    for (auto l = poly2.lines.begin(); l!=poly2.lines.end(); ++l){
        h = l->h;
        for(auto p = pts1.begin(); p!=pts1.end(); ++p) {
            if(*p * l->getnrm() < h) {h = *p * l->getnrm();}
        }
        h2start.push_back(h);//start big
        h2step.push_back((l->h - h)/steps);//shrink down
        l->h = h; //move the starting point now
    }
    ConvPolygon<float> ipoly = poly.intersect(poly2);
    mpl(ipoly.lines);
    vector<pt> ppts = ipoly.getpts();
    mpv(ppts);
     
    //start up GL
    if (!glfwInit())
            exit(EXIT_FAILURE);
    glfwSetErrorCallback(error_callback);
    
    //create window
    GLFWwindow* window = glfwCreateWindow(800, 800, "polygons", NULL, NULL); 
    if (!window)
    {
            glfwTerminate();
                exit(EXIT_FAILURE);
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    int stepsleft=steps;//crutch
    //main rendering loop
    while (!glfwWindowShouldClose(window))
    {
        //setup view
        float ratio;
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        ratio = width / (float) height;
        glViewport(0, 0, width, height);
        glClear(GL_COLOR_BUFFER_BIT);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(-ratio, ratio, -1.f, 1.f, 1.f, -1.f);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        
        //morph the polygons
        if(stepsleft > 0) {
            for (unsigned i = 0; i<poly.lines.size(); ++i){
               poly.lines[i].h += h1step[i];
            }
            for (unsigned i = 0; i<poly2.lines.size(); ++i){
               poly2.lines[i].h += h2step[i];
            }
            ipoly = poly.intersect(poly2);
            ppts = ipoly.getpts();
            stepsleft--; 
        } else { 
            glfwWaitEvents();
        }
        
        //draw them
        glBegin(GL_TRIANGLE_FAN);
        for(auto i=ppts.begin(); i!= ppts.end(); ++i) {
            glColor3f(1.f,0.f,0.f);
            glVertex3f(i->x/5,i->y/5,0.f);
        }
        glEnd();
        
        //output
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    
    //exit
    glfwDestroyWindow(window);
    glfwTerminate();
}


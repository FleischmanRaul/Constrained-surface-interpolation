#include <QApplication>
#include <iostream>
#include <ctime>
#include "GUI/MainWindow.h"

using namespace cagd;

using namespace std;

class Base
{
protected:
    int _b;
public:
    Base(int b):_b(b){}
    virtual int g() const
    {
        return _b;
    }
    virtual Base* clone() const
    {
        return new Base(*this);
    }
};


class Derived1:public Base
{
protected:
    int _d;
public:
    Derived1(int b, int d):Base(b),_d(d){}
    int g() const
    {
        return _b+_d;
    }
    Derived1* clone() const
    {
        return new Derived1(*this);
    }
};

class Derived2:public Base
{
protected:
    int _d;
public:
    Derived2(int b, int d):Base(b),_d(d){}
    int g() const
    {
        return _b*_d;
    }
    Derived2* clone() const
    {
        return new Derived2(*this);
    }

};

int main(int argc, char **argv)
{
    // creating an application object and setting one of its attributes
    QApplication app(argc, argv);
    srand(time(NULL));
//    Base *p = new Derived2(1,2);
//    cout << p->g() << endl;
//    Base *q = p->clone();
//    cout << q->g() << endl;
    // if you have installed a different version of Qt, it may happen that
    // the application attribute Qt::AA_UseDesktopOpenGL is not recognized
    // on Windows its existence is critical for our applications
    // on Linux or Mac you can uncomment this line
    app.setAttribute(Qt::AA_UseDesktopOpenGL, true);

    // creating a main window object
    MainWindow mwnd;
    mwnd.showMaximized();
    // running the application
    cout << "This is the end for you my friend" << endl;
    return app.exec();
}

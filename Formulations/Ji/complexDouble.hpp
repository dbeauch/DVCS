#ifndef COMPLEXDOUBLE_HPP
#define COMPLEXDOUBLE_HPP


class complexDouble {

private:
    double realPart;    
    double imaginaryPart;

public:
    complexDouble();
    complexDouble(double a, double b);
    ~complexDouble();
    double real();
    double imaginary();
    void setReal(double newReal);
    void setImaginary(double newImaginary);

};

complexDouble cdstar(complexDouble c, complexDouble d);
complexDouble conjugate(complexDouble z);


#endif
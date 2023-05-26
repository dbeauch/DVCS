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
    double getReal();
    double getImaginary();
    void subtract(complexDouble otherTerm);
    void add(complexDouble otherTerm);
    void multiplyByConstant(double c);
    complexDouble operator + (const complexDouble& obj){
        complexDouble temp;
        temp.realPart = realPart + obj.realPart;
        temp.imaginaryPart = imaginaryPart + obj.imaginaryPart;
        return temp;
    }
    complexDouble operator - (const complexDouble& obj){
        complexDouble temp;
        temp.realPart = realPart - obj.realPart;
        temp.imaginaryPart = imaginaryPart - obj.imaginaryPart;
        return temp;
    }
    complexDouble operator * (double myConst){
        complexDouble temp;
        temp.realPart = realPart * myConst;
        temp.imaginaryPart = imaginaryPart * myConst;
        return temp;
    }
    void operator = (const complexDouble& z) {
        realPart = z.realPart;
        imaginaryPart = z.imaginaryPart;
    }

};

complexDouble cdstar(complexDouble c, complexDouble d);
complexDouble conjugate(complexDouble z);


#endif
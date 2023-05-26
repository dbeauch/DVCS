#include "complexDouble.hpp"


complexDouble::complexDouble() {
    realPart = 0;
    imaginaryPart = 0;
}

complexDouble::complexDouble(double a, double b) {
    realPart = a;
    imaginaryPart = b;
}

complexDouble::~complexDouble(){

}

double complexDouble::real(){
    return realPart;
}

double complexDouble::imaginary(){
    return imaginaryPart;
}

void complexDouble::subtract(complexDouble otherTerm){
    realPart -= otherTerm.getReal();
    imaginaryPart -= otherTerm.getReal();
}

void complexDouble::add(complexDouble otherTerm){
    realPart += otherTerm.getReal();
    imaginaryPart += otherTerm.getReal();
}

void complexDouble::multiplyByConstant(double c){
    realPart *= c;
    imaginaryPart *= c;
}

void complexDouble::setReal(double newReal){
    realPart = newReal;
    return;
}

void complexDouble::setImaginary(double newImaginary){
    imaginaryPart = newImaginary;
    return;
}

double complexDouble::getReal(){
    return realPart;
}

double complexDouble::getImaginary(){
    return imaginaryPart;
}

complexDouble conjugate(complexDouble z){
    return complexDouble(z.real(), -1.0 * z.imaginary());
}

complexDouble cdstar(complexDouble c, complexDouble d){
    complexDouble dstar = conjugate(d);
    return complexDouble(c.real() * dstar.real() - c.imaginary() * dstar.imaginary(), c.real() * dstar.imaginary() + c.imaginary() * dstar.real());
}
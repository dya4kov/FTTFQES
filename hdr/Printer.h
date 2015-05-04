#pragma once
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

enum Pos {left, right};
/**
* @brief Simple formatted void printing interface.
*/
class Printer {
public:
	/**
	* @brief Print formatted string into stringstream.
	*/
	void printString(std::stringstream &ss,
				const char* data, 
				const int width = 20,
				const Pos pos = left, 
				const char fill = ' ');
	/**
	* @brief Print formatted string into stringstream.
	*/
	void printString(std::stringstream &ss,
				const std::string data, 
				const int width = 20,
				const Pos pos = left, 
				const char fill = ' ');
	/**
	* @brief Print formatted double into stringstream.
	*/
	void printDouble(std::stringstream &ss,
				const double value, 
				const int digitsToPrint = 3,
				const int width = 20,
				const Pos pos = left,
				const char fill = ' ');
	/**
	* @brief Print formatted scientific double into stringstream.
	*/
	void printSciDouble(std::stringstream &ss,
				const double value, 
				const int digitsToPrint = 3,
				const int width = 20,
				const Pos pos = left,
				const char fill = ' ');
	/**
	* @brief Print formatted integer into stringstream.
	*/
	void printInt(std::stringstream &ss,
				const int value,
				const int width = 20,
				const Pos pos = left,
				const char fill = ' ');
	/**
	* @brief Print formatted string into cout.
	*/
	void printString(const char* data, 
				const int width = 20,
				const Pos pos = left, 
				const char fill = ' ');
	/**
	* @brief Print formatted string into cout.
	*/
	void printString(const std::string data, 
				const int width = 20,
				const Pos pos = left, 
				const char fill = ' ');
	/**
	* @brief Print formatted double into cout.
	*/
	void printDouble(const double value, 
				const int digitsToPrint = 3,
				const int width = 20,
				const Pos pos = left,
				const char fill = ' ');
	/**
	* @brief Print formatted scientific double into cout.
	*/
	void printSciDouble(const double value, 
				const int digitsToPrint = 3,
				const int width = 20,
				const Pos pos = left,
				const char fill = ' ');
	/**
	* @brief Print formatted integer into cout.
	*/
	void printInt(const int value,
				const int width = 20,
				const Pos pos = left,
				const char fill = ' ');
	/**
	* @brief Print formatted string into filestream out.
	*/
	void printString(std::ofstream &out, 
				const char* data, 
				const int width = 20,
				const Pos pos = left, 
				const char fill = ' ');
	/**
	* @brief Print formatted string into filestream out.
	*/
	void printString(std::ofstream &out, 
				const std::string data, 
				const int width = 20,
				const Pos pos = left, 
				const char fill = ' ');
	/**
	* @brief Print formatted double into filestream out.
	*/
	void printDouble(std::ofstream &out, 
				const double value, 
				const int digitsToPrint = 3,
				const int width = 20,
				const Pos pos = left,
				const char fill = ' ');
	/**
	* @brief Print formatted scientific double into filestream out.
	*/
	void printSciDouble(std::ofstream &out, 
				const double value, 
				const int digitsToPrint = 3,
				const int width = 20,
				const Pos pos = left,
				const char fill = ' ');
	/**
	* @brief Print formatted integer into filestream out.
	*/
	void printInt(std::ofstream &out, 
				const int value,
				const int width = 20,
				const Pos pos = left,
				const char fill = ' ');
};

void Printer::printString(std::stringstream& ss,
			const char* data, 
			const int width,
			const Pos pos, 
			const char fill) {
	if (pos == left) ss.setf(std::ios::left);
	if (pos == right) ss.setf(std::ios::right);
    ss.width(width);
    ss.fill(fill);
    ss << data;
}

void Printer::printString(std::stringstream& ss,
			const std::string data, 
			const int width,
			const Pos pos, 
			const char fill) {
	if (pos == left)  ss.setf(std::ios::left);
	if (pos == right) ss.setf(std::ios::right);
    ss.width(width);
    ss.fill(fill);
    ss << data;
}

void Printer::printDouble(std::stringstream& ss,
			const double value, 
			const int digitsToPrint,
			const int width,
			const Pos pos,
			const char fill) {
	if (pos == left)  ss.setf(std::ios::left );
	if (pos == right) ss.setf(std::ios::right);
    ss.width(width);
    ss.precision(digitsToPrint);
    ss.fill(fill);
    ss << value;
}

void Printer::printSciDouble(std::stringstream& ss,
			const double value, 
			const int digitsToPrint,
			const int width,
			const Pos pos,
			const char fill) {
	if (pos == left)  ss.setf(std::ios::scientific | std::ios::left  | std::ios::showpos);
	if (pos == right) ss.setf(std::ios::scientific | std::ios::right | std::ios::showpos);
    ss.width(width);
    ss.precision(digitsToPrint);
    ss.fill(fill);
    ss << value;
}

void Printer::printInt(std::stringstream& ss,
			const int value,
			const int width,
			const Pos pos,
			const char fill) {
	if (pos == left)  ss.setf(std::ios::left );
	if (pos == right) ss.setf(std::ios::right);
    ss.width(width);
    ss.fill(fill);
    ss << value;
}

void Printer::printString(const char* data, 
			const int width,
			const Pos pos, 
			const char fill) {
	if (pos == left) std::cout.setf(std::ios::left);
	if (pos == right) std::cout.setf(std::ios::right);
    std::cout.width(width);
    std::cout.fill(fill);
    std::cout << data;
}

void Printer::printString(const std::string data, 
			const int width,
			const Pos pos, 
			const char fill) {
	if (pos == left) std::cout.setf(std::ios::left);
	if (pos == right) std::cout.setf(std::ios::right);
    std::cout.width(width);
    std::cout.fill(fill);
    std::cout << data;
}

void Printer::printDouble(const double value, 
			const int digitsToPrint,
			const int width,
			const Pos pos,
			const char fill) {
	if (pos == left)  std::cout.setf(std::ios::left );
	if (pos == right) std::cout.setf(std::ios::right);
    std::cout.width(width);
    std::cout.precision(digitsToPrint);
    std::cout.fill(fill);
    std::cout << value;
}

void Printer::printSciDouble(const double value, 
			const int digitsToPrint,
			const int width,
			const Pos pos,
			const char fill) {
	if (pos == left)  std::cout.setf(std::ios::scientific | std::ios::left  | std::ios::showpos);
	if (pos == right) std::cout.setf(std::ios::scientific | std::ios::right | std::ios::showpos);
    std::cout.width(width);
    std::cout.precision(digitsToPrint);
    std::cout.fill(fill);
    std::cout << value;
}

void Printer::printInt(const int value,
			const int width,
			const Pos pos,
			const char fill) {
	if (pos == left)  std::cout.setf(std::ios::left );
	if (pos == right) std::cout.setf(std::ios::right);
    std::cout.width(width);
    std::cout.fill(fill);
    std::cout << value;
}

void Printer::printString(std::ofstream &out, 
			const char* data, 
			const int width,
			const Pos pos, 
			const char fill) {
	if (pos == left)  out.setf(std::ios::left);
	if (pos == right) out.setf(std::ios::right);
    out.width(width);
    out.fill(fill);
    out << data;
}

void Printer::printString(std::ofstream &out, 
			const std::string data, 
			const int width,
			const Pos pos, 
			const char fill) {
	if (pos == left)  out.setf(std::ios::left);
	if (pos == right) out.setf(std::ios::right);
    out.width(width);
    out.fill(fill);
    out << data;
}

void Printer::printDouble(std::ofstream &out, 
			const double value, 
			const int digitsToPrint,
			const int width,
			const Pos pos,
			const char fill) {
	if (pos == left)  out.setf(std::ios::left );
	if (pos == right) out.setf(std::ios::right);
    out.width(width);
    out.precision(digitsToPrint);
    out.fill(fill);
    out << value;
}

void Printer::printSciDouble(std::ofstream &out, 
			const double value, 
			const int digitsToPrint,
			const int width,
			const Pos pos,
			const char fill) {
	if (pos == left)  out.setf(std::ios::scientific | std::ios::left  | std::ios::showpos);
	if (pos == right) out.setf(std::ios::scientific | std::ios::right | std::ios::showpos);
    out.width(width);
    out.precision(digitsToPrint);
    out.fill(fill);
    out << value;
}

void Printer::printInt(std::ofstream &out, 
			const int value,
			const int width,
			const Pos pos,
			const char fill) {
	if (pos == left)  out.setf(std::ios::left );
	if (pos == right) out.setf(std::ios::right);
    out.width(width);
    out.fill(fill);
    out << value;
}
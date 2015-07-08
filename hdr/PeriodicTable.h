#pragma once
#include <iostream>
#include <fstream>
#include <map>
#include "stringUtils.h"
#include "Printer.h"

#define PERIODIC_TABLE_FILE "res/PeriodicTableSimple.dat"
#define PERIODIC_TABLE_TEST "out/PeriodicTableTest.dat"

class PeriodicTable {
public:
	struct Data {
		double mass;
		double Z;
	};
	void print();
	PeriodicTable();
	Data& operator[] (std::string element) {
		return dataMap[element];
	}
private:
	void load();
	std::map<std::string, Data> dataMap;
	Printer printer;
};

PeriodicTable::PeriodicTable() {
	load();
}

void PeriodicTable::load() {
	std::ifstream datafile(PERIODIC_TABLE_FILE, std::ios::in);
	std::string elName;
	Data elem;
	if (!datafile.good()) {
		std::cerr << "cannot read periodic table from" 
		<< PERIODIC_TABLE_FILE 
		<< std::endl;
	}
	datafile >> elName;
	datafile >> elName;
	datafile >> elName;
    while (!datafile.eof()) {
    	datafile >> elName;
    	datafile >> elem.Z;
    	datafile >> elem.mass;
    	dataMap[elName] = elem;
    }
    datafile.close();
}

void PeriodicTable::print() {
	std::ofstream outfile(PERIODIC_TABLE_TEST, std::ios::out);
	if (!outfile.good()) {
		std::cerr << "cannot write periodic table in " 
		<< PERIODIC_TABLE_FILE 
		<< std::endl;
	}
	std::map<std::string, Data>::iterator iter;
	printer.printString(outfile, "Name", 15, left);
	printer.printString(outfile, "Charge", 15, left);
	printer.printString(outfile, "Mass", 15, left);
	outfile << std::endl;
	for (iter = dataMap.begin(); iter != dataMap.end(); iter++) {
		printer.printString(outfile, iter->first, 15, left);
		printer.printDouble(outfile, iter->second.Z, 5, 15, left);
		printer.printDouble(outfile, iter->second.mass, 5, 15, left);
		outfile << std::endl;
	}
    outfile.close();
}
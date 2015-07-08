#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>

#define BUF_SIZE 100

void readInput(std::string filename,
			   std::string& element,
			   std::string& volume, 
	           std::string& temperature,
	           std::string& output,
	           double& tolerance) {
	std::ifstream inputfile(("./in/" + filename).c_str(), std::ios::in);
	char* buffer = new char[BUF_SIZE];
	std::string data;
	while (!isalnum(buffer[0]))
		inputfile.getline(buffer, BUF_SIZE);
	data = buffer;
	element = split(data, '#')[0];
	element = split(data, ' ')[0];
	std::cout << "element: " << element << std::endl;
	inputfile.getline(buffer, BUF_SIZE);
	while (!isalnum(buffer[0]))
		inputfile.getline(buffer, BUF_SIZE);
	data = buffer;
	volume = split(data, '#')[0];
	std::cout << "volume: " << volume << std::endl;
	inputfile.getline(buffer, BUF_SIZE);
	while (!isalnum(buffer[0]))
		inputfile.getline(buffer, BUF_SIZE);
	data = buffer;
	temperature = split(data, '#')[0];
	std::cout << "temperature: " << temperature << std::endl;
	inputfile.getline(buffer, BUF_SIZE);
	while (!isalnum(buffer[0]))
		inputfile.getline(buffer, BUF_SIZE);
	data = buffer;
	output = split(data, '#')[0];
	std::cout << "output data: " << output << std::endl;
	inputfile.getline(buffer, BUF_SIZE);
	while (!isalnum(buffer[0]))
		inputfile.getline(buffer, BUF_SIZE);
	data = buffer;
	data = split(data, '#')[0];
	data = split(data, ' ')[0];
	tolerance = atof(data.c_str());
	std::cout << "tolerance: " << tolerance << std::endl;
}
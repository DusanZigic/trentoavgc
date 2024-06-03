#include <iostream>
#include <string>
#include <cstring>
#include <cctype>
#include <sstream>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
#include <fstream>
#include <iomanip>

static std::string srcDir;
static size_t eventIDLow, eventIDHigh;
static double xRangeLow, xRangeHigh, xBinW;
static double yRangeLow, yRangeHigh, yBinW;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int getSDSize(size_t &xsize, size_t &ysize)
//function that determines size of entropy density file
{
	std::string pathIn = srcDir + std::to_string(eventIDLow) + ".dat";
	std::ifstream fileIn(pathIn, std::ios_base::in);
	if (!fileIn.is_open()) {
		std::cerr << "Error: unable to open entropy density file when determining grid sizes. Aborting..." << std::endl;
		return -1;
	}

	xsize = 0; std::vector<size_t> ysizes;

	std::string line; double buff;
	while (std::getline(fileIn, line)) {
		std::stringstream ss(line);
		ysize = 0; while (!ss.eof()) {ss >> buff; ++ysize;}
		ysizes.push_back(ysize);
		++xsize;
	}

	std::sort(ysizes.begin(), ysizes.end()); ysizes.erase(unique(ysizes.begin(), ysizes.end()), ysizes.end());
	if (ysizes.size() != 1) {
		std::cerr << "Error: not all entropy density sizes are the same. Aborting..." << std::endl;
		return -2;
	}

	ysize = ysizes[0];

	fileIn.close();

	return 1;
}

int addSD(std::vector<std::vector<double>> &sd, size_t event_id)
//function that adds entropy density
{
	std::string pathIn = srcDir + std::to_string(event_id) + ".dat";
	std::ifstream fileIn(pathIn);
	if (!fileIn.is_open()) {
		std::cerr << "Error: unable to open sd input file for event " + std::to_string(event_id) + ". Aborting..." << std::endl;
		return -1;
	}

	std::string line; double buff; size_t x1Cnt = 0, x2Cnt = 0;

	while (std::getline(fileIn, line)) {
		std::stringstream ss(line);
		while (!ss.eof()) {ss >> buff; sd[x1Cnt][x2Cnt] += buff; ++x2Cnt;}
		++x1Cnt;
		x2Cnt=0;
	}

	fileIn.close();

	return 1;
}

int normSD(std::vector<std::vector<double>> &sd)
//function that norms entropy density
{
	double eventN = static_cast<double>(eventIDHigh - eventIDLow + 1);

	for (size_t ix1=0; ix1<sd.size(); ++ix1)
		for (size_t ix2=0; ix2<sd[ix1].size(); ++ix2)
			sd[ix1][ix2] /= eventN;

	return 1;
}

int exportSD(std::vector<std::vector<double>> &sd)
//function that exports averaged entropy density:
{
	std::string pathOut = "./sdavg.dat";
	std::ofstream fileOut(pathOut, std::ios_base::out);
	if (!fileOut.is_open()) {
		std::cerr << "Error: unable to open output file. Aborting..." << std::endl;
		return -1;
	}

	for (size_t ix1=0; ix1<sd.size(); ++ix1) {
		for (size_t ix2=0; ix2<sd[ix1].size(); ++ix2) {
			fileOut << sd[ix1][ix2] << " ";
		}
		fileOut << "\n";
	}

	fileOut.close();

	return 1;	
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<double> xBins, yBins, xGrid, yGrid;

int generateGridsBins()
//function that generates bins and grids:
{
	{
		double x_0 = xRangeLow - abs(xBinW)/2.0, x_step = xBinW, x_max = xRangeHigh + abs(xBinW)/2.0;
		double x   = x_0; while (x <= x_max) {xBins.push_back(x); x += x_step;}
		
		double y_0 = yRangeLow - abs(yBinW)/2.0, y_step = yBinW, y_max = yRangeHigh + abs(yBinW)/2.0;
		double y   = y_0; while (y <= y_max) {yBins.push_back(y); y += y_step;}
	}
	
	{
		double x_0 = xRangeLow, x_step = xBinW, x_max = xRangeHigh;
		double x   = x_0; while (x <= x_max) {xGrid.push_back(x); x += x_step;}
		
		double y_0 = yRangeLow, y_step = yBinW, y_max = yRangeHigh;
		double y   = y_0; while (y <= y_max) {yGrid.push_back(y); y += y_step;}
	}

	return 1;
}

int loadBCP(std::vector<std::vector<double>> &bcp, size_t event_id)
//function that loads binary collision points:
{
	std::string pathIn = srcDir + "bcp" + std::to_string(event_id) + ".dat";
	std::ifstream fileIn(pathIn, std::ios_base::in);
	if (!fileIn.is_open()) {
		std::cerr << "Error: unable to open bcp file for event " + std::to_string(event_id) + ". Aborting..." << std::endl;
		return -1;
	}

	std::string line; double buffX, buffY;

	while (std::getline(fileIn, line)) {
		if (line.rfind("#", 0) == 0) continue;
		std::stringstream ss(line);
		ss >> buffX; ss >> buffY; bcp.push_back(std::vector<double>{buffX, buffY});
	}

	fileIn.close();

	return 1;
}

int countBCP(std::vector<std::vector<double>> &hist, const std::vector<std::vector<double>> &bcp)
//function that counts binary collision points into a histogram:
{
	for (const auto& xy : bcp) {
		double x = xy[0];
		size_t i=0; while(x > xBins[i]) ++i; --i;
		if ((i<0) || (i>xGrid.size()-1)) {
			std::cerr << "Error: x point value outside histogram's range. Aborting..." << std::endl;
			return -1;
		}

		double y = xy[1];
		size_t j=0; while(y > yBins[j]) ++j; --j;
		if ((j<0) || (j>yGrid.size()-1)) {
			std::cerr << "Error: y point value outside histogram's range. Aborting..." << std::endl;
			return -2;
		}

		hist[i][j] += 1.0;
	}

	return 1;
}

int normHIST(std::vector<std::vector<double>> &hist)
//function that norms histogram:
{
	double histMax = 0.0;

	for (const auto& h: hist) {
		double rowMax = *max_element(begin(h), end(h));
		histMax = std::max(histMax, rowMax);
	}
	
	double tenCnt = 1.0;
	while(histMax >= 10) {
		tenCnt += 1;
		histMax /= 10.0;
	}
	double histNorm = pow(10.0, tenCnt);

	for (size_t iH=0; iH<hist.size(); ++iH) {
		std::for_each(hist[iH].begin(), hist[iH].end(), [histNorm](double &el){el /= histNorm;});
	}

	return 1;
}

int exportHIST(const std::vector<std::vector<double>> &hist)
//function that exports histogram
{
	std::string pathOut = "./bcdensity.dat";

	std::ofstream fileOut(pathOut, std::ios_base::out);
	if (!fileOut.is_open()) {
		std::cerr << "Error: unable to open output file. Aborting..." << std::endl;
		return -1;
	}

	for (size_t iX=0; iX<xGrid.size(); ++iX)
		for (size_t iY=0; iY<yGrid.size(); ++iY)
			fileOut << std::fixed << std::setw(5) << std::setprecision(2) <<	xGrid[iX] << " "
					<< std::fixed << std::setw(5) << std::setprecision(2) <<    yGrid[iY] << " "
					<< std::fixed << std::setw(7) << std::setprecision(6) << hist[iX][iY] << std::endl;

	fileOut.close();

	return 1;
}

int main(int argc, char const *argv[])
{
	if (argc != 10) {
		std::cerr << "usage: ./trentoavgc src_dir eventIDLow eventIDHigh xRangeLow xRangeHigh xBinW yRangeLow yRangeHigh yBinW" << std::endl;
		return -1;
	}

	srcDir = argv[1]; if (srcDir.back() != '/') srcDir += "/";
	eventIDLow = std::atoi(argv[2]); eventIDHigh = std::atoi(argv[3]);
	xRangeLow  = std::atof(argv[4]); xRangeHigh  = std::atof(argv[5]); xBinW = std::atof(argv[6]);
	yRangeLow  = std::atof(argv[7]); yRangeHigh  = std::atof(argv[8]); yBinW = std::atof(argv[9]);

	size_t xSDsize, ySDsize; if (getSDSize(xSDsize, ySDsize) != 1) return -2;
	std::vector<std::vector<double>> SD(xSDsize, std::vector<double>(ySDsize, 0.0));

	if (generateGridsBins() != 1) return -3;
	std::vector<std::vector<double>> HIST(xGrid.size(), std::vector<double>(yGrid.size(), 0.0));

	for (size_t eventID=eventIDLow; eventID<=eventIDHigh; eventID++) {
		if (addSD(SD, eventID)    != 1) return -4;
		std::vector<std::vector<double>> BCP;
		if (loadBCP(BCP, eventID) != 1) return -5;
		if (countBCP(HIST, BCP)   != 1) return -6;
	}

	if (normSD(SD)       != 1) return -7;
	if (exportSD(SD)     != 1) return -8;
	if (normHIST(HIST)   != 1) return -9;
	if (exportHIST(HIST) != 1) return -10;

	return 0;
}
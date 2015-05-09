/**
Course: CS536 - Computer Graphics
Student: Zhichao Cao
Email: zc77@drexel.edu
Title: CS536 - Extra credit
*/

#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
using namespace std ;

string readFile(string strPath);
double calSlope(int qx, int qy, int rx, int ry);
void bresenham(int qx, int qy, int rx, int ry, int sign, int index, double scal, int rota, int xTran, int yTran);
void writePixel(int x, int y, int sign, int index, int length);
void transform(double& x, double& y, int sign, int index, double scal, int rota, int xTran, int yTran, int umin, int vmin, int umax, int vmax, int xl, int yl, int xu, int yu);
int calXpoint(int qx, int rx);
void drawPic(string & strOutput, int countIndex, int xLowerBound, int xUpperBound, int yLowerBound, int yUpperBound, int viewXLowerBound, int viewYLowerBound, int viewXUpperBound, int viewYUpperBound);
void scanLineFill(int viewXLowerBound, int viewYLowerBound, int viewXUpperBound, int viewYUpperBound, int xyLnIndex);
void CalculateBezierPoint(int &x, int &y, double t, double px1, double py1, double px2, double py2, double px3, double py3, double px4, double py4);

double xyOp[4000][4];  //Set of pairs of points
double xyOpSlope[4000];  //Set of slope of lines
double xOpPath[4000][4000];  //Set x-coop of points on line
double yOpPath[4000][4000];  //Set y-coop of points on line
double typeOpPath[4000];  //Type of slope of lines
double lenOpPath[4000];  //Number of points on line
double lenOpPathBack[4000];
int realWorldCoor[501][501];  //For cooperation of real world
int viewWorldCoor[501][501];  //For viewport
//int scanLine[501][501];
//int scanCount[501];
//int scanLineCount[501];
//int scanLineSet[501][501];
//int sortedInter[501][501];
//int yMax = -1;
//int yMin = 2001;  //For scan-line filling algorithm
int xCurveOpPath[4000][4000];
int yCurveOpPath[4000][4000];
int lenCurveOpPath[4000];
int xyCurLnOpPath[4000][4];
int xyCurLnOpPathIndex;
//int xCurLnOpPath[8000][5000];
//int yCurLnOpPath[5000][5000];

const double PI = 3.14159265358979323846;

int main(int argc, char* argv[])
{
	string strOutput = "";  //String for output xpm file
	//string strOutputFile = "out.xpm";  //Output file name
	string strInputFile = "ExtraCredit.ps";  //-f
	double scalingFactor = 1.0;  //-s
	int cClockwiseRotation = 0;  //-r
	int xTranslation = 0;  //-m
	int yTranslation = 0;  //-n
	int xLowerBound = 0;  //-a
	int yLowerBound = 0;  //-b
	int xUpperBound = 250;  //-c
	int yUpperBound = 250;  //-d
	int viewXLowerBound = 0;  //-j
	int viewYLowerBound = 0;  //-k
	int viewXUpperBound = 200;  //-o
	int viewYUpperBound = 200;  //-p
	string strPS = "";  //Content string in the .ps file
	string strSetPS[500];  //Split strPS string
	int stIndex = 0;  
	int edIndex = 0;
	int stIndexOp = 0;
	string tempStrOp = "";
	int xyLnIndex = 0;  //Number of lines
	int xyOpIndex = 0;  //Index for x,y-coop, for line command
	int fstPointIndex = 0;  //Temp index for first point of moveto and lineto command
	int sndPointIndex = 0;  //Temp index for second point of moveto and lineto command
	int countIndex = 0;  //Number of lines in .ps file
	string strLength = "";  //Length of pic
	string strWidth = "";  //Width of pic
	int lineType = 0;  //Type of command in .ps file. 1 for line, 2 for moveto, 3 for lineto, 4 for curveto
	double xyCurveOp[5000][8];
	int xyCurveLnIndex = 0;
	double curvePara = 0.05;

	if(argv[1] == ">" && argv[2] == "out.xpm")
	{
		//./CG_hw1 -f hw1.ps -a 0 -b 0 -c 499 -d 499 -s 1.0 -m 0 -n 0 -r 0 > out.xpm
		strInputFile = "ExtraCredit.ps";
		scalingFactor = 1.0;
		cClockwiseRotation = 0;
		xTranslation = 0;
		yTranslation = 0;
		xLowerBound = 0;
		yLowerBound = 0;
		xUpperBound = 500;
		yUpperBound = 500;
		viewXLowerBound = 0;
		viewYLowerBound = 0;
		viewXUpperBound = 200;
		viewYUpperBound = 200;
		curvePara = 0.05;
		//strOutputFile = "out.xpm";
	}
	else
	{
		for(int i = 0; i < argc; i++)
		{
			string tempStr = "";
			tempStr = argv[i];
			if(tempStr == "-f")
			{
				strInputFile = argv[i+1];
			}
			else if(tempStr == "-s")
			{
				scalingFactor = atof(argv[i+1]);
			}
			else if(tempStr == "-r")
			{
				cClockwiseRotation = atoi(argv[i+1]);
			}
			else if(tempStr == "-m")
			{
				xTranslation = atoi(argv[i+1]);
			}
			else if(tempStr == "-n")
			{
				yTranslation = atoi(argv[i+1]);
			}
			else if(tempStr == "-a")
			{
				xLowerBound = atoi(argv[i+1]);
			}
			else if(tempStr == "-b")
			{
				yLowerBound = atoi(argv[i+1]);
			}
			else if(tempStr == "-c")
			{
				xUpperBound = atoi(argv[i+1]);
			}
			else if(tempStr == "-d")
			{
				yUpperBound = atoi(argv[i+1]);
			}
			else if(tempStr == "-j")
			{
				viewXLowerBound = atoi(argv[i+1]);
			}
			else if(tempStr == "-k")
			{
				viewYLowerBound = atoi(argv[i+1]);
			}
			else if(tempStr == "-o")
			{
				viewXUpperBound = atoi(argv[i+1]);
			}
			else if(tempStr == "-p")
			{
				viewYUpperBound = atoi(argv[i+1]);
			}
			else if(tempStr == ">")
			{
				//strOutputFile = argv[i+1];
			}			
			else if(tempStr == "-L")
			{
				curvePara = atof(argv[i+1]);
			}
		}
		//cout << strInputFile << " " << scalingFactor << " "  << cClockwiseRotation << " "  << xTranslation << " "  << yTranslation << " "  << xLowerBound << " "  << yLowerBound << " "  << xUpperBound << " "  << yUpperBound << " "  << endl;
	}

	strPS = readFile(strInputFile);
	//cout << strPS << endl;

	for(int i = 0; i < strPS.length(); i++)
	{
		if(strPS[i] == '\n')
		{
			strSetPS[countIndex] = strPS.substr(stIndex, i - stIndex);
			stIndex = i + 1;
			countIndex++;
		}
	}

	for (int i = 0; i < countIndex; i++)
	{
		//cout << strSetPS[i] << endl;
		if(strSetPS[i].substr(0, 8) == "%%%BEGIN")
		{
			stIndex = i + 1;
		}
		else if(strSetPS[i].substr(0, 6) == "%%%END")
		{
			edIndex = i - 1;
		}
		else
		{
			continue;
		}
	}

	for (int i = stIndex; i <= edIndex; i++)
	{
		std::size_t posLineLow = strSetPS[i].find("line");
		std::size_t posLineUp = strSetPS[i].find("Line");
		std::size_t posLineToLow = strSetPS[i].find("lineto");
		std::size_t posLineToUp = strSetPS[i].find("Lineto");
		std::size_t posLineToUpAll = strSetPS[i].find("LineTo");
		std::size_t posMoveToLow = strSetPS[i].find("moveto");
		std::size_t posMoveToUp = strSetPS[i].find("Moveto");
		std::size_t posMoveToUpAll = strSetPS[i].find("MoveTo");
		std::size_t posCurveToLow = strSetPS[i].find("curveto");
		std::size_t posCurveToUp = strSetPS[i].find("Curveto");
		std::size_t posCurveToUpAll = strSetPS[i].find("CurveTo");
		std::size_t posStrokeLow = strSetPS[i].find("stroke");
		std::size_t posStrokeUp = strSetPS[i].find("Stroke");
		if(posLineToLow != std::string::npos || posLineToLow != std::string::npos || posLineToUpAll != std::string::npos)
		{
			lineType = 3;
		}
		else if(posMoveToLow != std::string::npos || posMoveToUp != std::string::npos || posMoveToUpAll != std::string::npos)
		{
			lineType = 2;
		}
		else if((posLineLow != std::string::npos || posLineUp != std::string::npos) && (posLineToLow == std::string::npos || posLineToLow == std::string::npos || posLineToUpAll == std::string::npos))
		{
			lineType = 1;
		}
		else if(posCurveToLow != std::string::npos || posCurveToUp != std::string::npos || posCurveToUpAll != std::string::npos)
		{
			lineType = 4;
		}
		//else if(posStrokeLow != std::string::npos || posStrokeUp != std::string::npos)
		else
		{
			fstPointIndex = 0;
			continue;
		}
		for (int j = 0; j < strSetPS[i].length(); j++)
		{			
			if(lineType == 2)
			{
				if(strSetPS[i][j] == ' ' && strSetPS[i][stIndexOp] != '-')
				{
					tempStrOp = strSetPS[i].substr(stIndexOp, j - stIndexOp);
					//xyOp[xyLnIndex][fstPointIndex] = atof(tempStrOp.c_str());
					xyCurveOp[xyCurveLnIndex][fstPointIndex] = atof(tempStrOp.c_str());
					stIndexOp = j + 1;
					fstPointIndex++;
				}
				else if(strSetPS[i][j] == ' ' && strSetPS[i][stIndexOp] == '-')
				{
					tempStrOp = strSetPS[i].substr(stIndexOp + 1, j - stIndexOp);
					//xyOp[xyLnIndex][fstPointIndex] = atof(tempStrOp.c_str());
					xyCurveOp[xyCurveLnIndex][fstPointIndex] = 0 - atof(tempStrOp.c_str());
					stIndexOp = j + 1;
					fstPointIndex++;
				}
			}
			else if(lineType == 1)
			{
				if(strSetPS[i][j] == ' ')
				{
					tempStrOp = strSetPS[i].substr(stIndexOp, j - stIndexOp);
					xyOp[xyLnIndex][xyOpIndex] = atof(tempStrOp.c_str());
					stIndexOp = j + 1;
					xyOpIndex++;
				}
			}
			else if(lineType == 4)
			{
				if(strSetPS[i][j] == ' ' && strSetPS[i][j] != '-')
				{
					tempStrOp = strSetPS[i].substr(stIndexOp, j - stIndexOp);
					xyCurveOp[xyCurveLnIndex][fstPointIndex] = atof(tempStrOp.c_str());
					stIndexOp = j + 1;
					fstPointIndex++;
				}
				if(strSetPS[i][j] == ' ' && strSetPS[i][j] == '-')
				{
					tempStrOp = strSetPS[i].substr(stIndexOp + 1, j - stIndexOp);
					xyCurveOp[xyCurveLnIndex][fstPointIndex] = atof(tempStrOp.c_str());
					stIndexOp = j + 1;
					fstPointIndex++;
				}
			}
		}
		xyOpIndex = 0;
		if(fstPointIndex != 2)
		{
			fstPointIndex = 0;
		}
		stIndexOp = 0;
		if(lineType == 1)
		{			
			xyLnIndex++;
		}
		if(lineType == 4)
		{
			xyCurveLnIndex++;
		}
	}
	
	/*for(int i = 0; i < xyLnIndex; i++)
	{
		cout << xyOp[i][0] << " " << xyOp[i][1] << " "  << xyOp[i][2] << " "  << xyOp[i][3] << " "  << endl;
	}*/
	/*for(int i = 0; i < xyCurveLnIndex; i++)
	{
		cout << i <<endl;
		cout << xyCurveOp[i][0] << " " << xyCurveOp[i][1] << " "  << xyCurveOp[i][2] << " "  << xyCurveOp[i][3] << " "  
				<< xyCurveOp[i][4] << " " << xyCurveOp[i][5] << " "  << xyCurveOp[i][6] << " "  << xyCurveOp[i][7] << " "<< endl;
	}*/

	/*xScale = static_cast<double>(viewXUpperBound - viewXLowerBound) / (xUpperBound - xLowerBound);
	yScale = static_cast<double>(viewYUpperBound - viewYLowerBound) / (yUpperBound - yLowerBound);*/

	for (int i = 0; i < xyLnIndex; i++)
	{
		double tempX, tempY = 0;			
		transform(xyOp[i][0], xyOp[i][1], 0, i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation, viewXLowerBound, viewYLowerBound, viewXUpperBound, viewYUpperBound, xLowerBound, yLowerBound, xUpperBound, yUpperBound);
		transform(xyOp[i][2], xyOp[i][3], 0, i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation, viewXLowerBound, viewYLowerBound, viewXUpperBound, viewYUpperBound, xLowerBound, yLowerBound, xUpperBound, yUpperBound);
		if(calXpoint(xyOp[i][0], xyOp[i][2]) == -1)
		{
			//xyOpSlope[i] = calSlope(xyOp[i][2], xyOp[i][3], xyOp[i][0], xyOp[i][1]);
			tempX = xyOp[i][0];
			tempY = xyOp[i][1];
			xyOp[i][0] = xyOp[i][2];
			xyOp[i][1] = xyOp[i][3];
			xyOp[i][2] = tempX;
			xyOp[i][3] = tempY;
		}
		//else
		//{
		//	xyOpSlope[i] = calSlope(xyOp[i][0], xyOp[i][1], xyOp[i][2], xyOp[i][3]);
		//}
		xyOpSlope[i] = calSlope(xyOp[i][0], xyOp[i][1], xyOp[i][2], xyOp[i][3]);
		
		//cout << xyOp[i][0] << " " << xyOp[i][1] << " "  << xyOp[i][2] << " "  << xyOp[i][3] << " " << xyOpSlope[i] << endl;
		
		if(xyOpSlope[i] >= 0 && xyOpSlope[i] <= 1)
		{
			typeOpPath[i] = 0;
			//No need for modification
			bresenham(xyOp[i][0], xyOp[i][1], xyOp[i][2], xyOp[i][3], typeOpPath[i], i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation);
		}
		else if(xyOpSlope[i] >= -1 && xyOpSlope[i] < 0)
		{
			typeOpPath[i] = 1;
			//Switch y1 and y2 to make 0<=slope<=1
			bresenham(xyOp[i][0], -xyOp[i][1], xyOp[i][2], -xyOp[i][3], typeOpPath[i], i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation);
		}
		else if(xyOpSlope[i] > 1)
		{
			typeOpPath[i] = 2;
			//(y1, x1)(y2, x2) to make 0<=slope<=1
			bresenham(xyOp[i][1], xyOp[i][0], xyOp[i][3], xyOp[i][2], typeOpPath[i], i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation);
		}
		else if(xyOpSlope[i] < -1)
		{
			typeOpPath[i] = 3;
			//(y2, -x2) and (y1, -x1) to make 0<=slope<=1
			bresenham(xyOp[i][3], -xyOp[i][2], xyOp[i][1], -xyOp[i][0], typeOpPath[i], i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation);
		}
		//else if(xyOpSlope[i] == 1000)
		//{
		//	for (int j = xyOp[i][1]; j <= xyOp[i][3]; j++)
		//	{
		//		writePixel(xyOp[i][0], j, 4, i, (j - xyOp[i][1]));
		//	}
		//	//bresenham(xyOp[i][0], xyOp[i][1], xyOp[i][2], xyOp[i][3], 4, i);
		//}
		//else if(xyOpSlope[i] == 1000)
		//{
		//	for (int j = xyOp[i][3]; j <= xyOp[i][1]; j++)
		//	{
		//		writePixel(xyOp[i][0], j, 5, i, (j - xyOp[i][3]));
		//	}
		//	//bresenham(xyOp[i][0], xyOp[i][1], xyOp[i][2], xyOp[i][3], 5, i);
		//}
	}
	for(int k = 0; k < countIndex; k++)
	{
		int xCoop = 0;
		int yCoop = 0;
		for (int l = 0; l < lenOpPath[k]; l++)
		{					
			xCoop = xOpPath[k][l];
			yCoop = yOpPath[k][l];
			realWorldCoor[xCoop][yCoop] = 1;
			if((xCoop >= viewXLowerBound && xCoop <= viewXUpperBound) && (yCoop >= viewYLowerBound && yCoop <= viewYUpperBound))
			{
				viewWorldCoor[xCoop][yCoop] = 1;
				/*if(yCoop > yMax)
				{
					yMax = yCoop;
				}
				else if(yCoop < yMin)
				{
					yMin = yCoop;
				}*/
			}
		}
	}

	for (int i = 0; i < xyCurLnOpPathIndex; i++)
	{
		lenOpPath[i] = 0;
	}

	for(int i = 0; i < xyCurveLnIndex; i++)
	{
		//Bezier curve drawing
		transform(xyCurveOp[i][0], xyCurveOp[i][1], 0, i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation, viewXLowerBound, viewYLowerBound, viewXUpperBound, viewYUpperBound, xLowerBound, yLowerBound, xUpperBound, yUpperBound);
		transform(xyCurveOp[i][2], xyCurveOp[i][3], 0, i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation, viewXLowerBound, viewYLowerBound, viewXUpperBound, viewYUpperBound, xLowerBound, yLowerBound, xUpperBound, yUpperBound);
		transform(xyCurveOp[i][4], xyCurveOp[i][5], 0, i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation, viewXLowerBound, viewYLowerBound, viewXUpperBound, viewYUpperBound, xLowerBound, yLowerBound, xUpperBound, yUpperBound);
		transform(xyCurveOp[i][6], xyCurveOp[i][7], 0, i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation, viewXLowerBound, viewYLowerBound, viewXUpperBound, viewYUpperBound, xLowerBound, yLowerBound, xUpperBound, yUpperBound);

		/*double x = 0.0;
		double y = 0.0;
		for(int j = 0; j < 1; j += 3)
		{
			double px1 = xyCurveOp[i][0];
			double py1 = xyCurveOp[i][1];
			double px2 = xyCurveOp[i][2];
			double py2 = xyCurveOp[i][3];
			double px3 = xyCurveOp[i][4];
			double py3 = xyCurveOp[i][5];
			double px4 = xyCurveOp[i][6];
			double py4 = xyCurveOp[i][7];
		}*/
		//List drawingPoints = new List();
		//for(int i = 0; i < controlPoints.Count - 3; i+=3)
		//{
		//	Vector3 p0 = controlPoints[i];
		//	Vector3 p1 = controlPoints[i + 1];
		//	Vector3 p2 = controlPoints[i + 2];
		//	Vector3 p3 = controlPoints[i + 3];    
 
		//	if(i == 0) //Only do this for the first endpoint.
		//				//When i != 0, this coincides with the end
		//				//point of the previous segment
		//	{
		//	drawingPoints.Add(CalculateBezierPoint(0, p0, p1, p2, p3));
		//	}    
 
		//	for(int j = 1; j <= SEGMENTS_PER_CURVE; j++)
		//	{
		//	float t = j / (float) SEGMENTS_PER_CURVE;
		//	drawingPoints.Add(CalculateBezierPoint(t, p0, p1, p2, p3));
		//	}
		//}

		int n = 3;
		int k = 0;
		int sigPlus = 0;
		for(double u = 0; u <= 1.00001; u += curvePara)
		{
			//cout<<": "<<u<<endl;
			int x = pow((1 - u), 3.0) * xyCurveOp[i][0] + 3 * u * pow((1 - u), 2.0) * xyCurveOp[i][2] + 3 * pow(u, 2.0) * (1 - u) * xyCurveOp[i][4] + pow(u, 3.0) * xyCurveOp[i][6];
			int y = pow((1 - u), 3.0) * xyCurveOp[i][1] + 3 * u * pow((1 - u), 2.0) * xyCurveOp[i][3] + 3 * pow(u, 2.0) * (1 - u) * xyCurveOp[i][5] + pow(u, 3.0) * xyCurveOp[i][7];
			xCurveOpPath[i][k] = x;
			yCurveOpPath[i][k] = y;
			k++;
			if((u + curvePara) > 1.00001 && sigPlus == 0)
			{
				u = 1.0 - curvePara;
				sigPlus++;
			}
		}
		lenCurveOpPath[i] = k;		
		/*x = (1-t)^3 *x0 + 3*t*(1-t)^2 *x1 + 3*t^2*(1-t) *x2 + t^3 *x3 
		y = (1-t)^3 *y0 + 3*t*(1-t)^2 *y1 + 3*t^2*(1-t) *y2 + t^3 *y3*/
	}
	/*for(int i = 0; i < xyCurveLnIndex; i++)
	{
		cout<< i << endl;
		cout << xyCurveOp[i][0] << " " << xyCurveOp[i][1] << " "  << xyCurveOp[i][2] << " "  << xyCurveOp[i][3] << " "  
				<< xyCurveOp[i][4] << " " << xyCurveOp[i][5] << " "  << xyCurveOp[i][6] << " "  << xyCurveOp[i][7] << " "<< endl;
		for (int j = 0; j < lenCurveOpPath[i]; j++)
		{
			cout << j << "::" << xCurveOpPath[i][j] << ", " << yCurveOpPath[i][j] << endl;
		}
		cout << endl;
	}*/
	double segment = 1 / curvePara;
	if(segment >= 285)
	{
		//Just output
		for(int k = 0; k < xyCurveLnIndex; k++)
		{
			int xCoop = 0;
			int yCoop = 0;
			for (int l = 0; l < lenCurveOpPath[k]; l++)
			{					
				xCoop = xCurveOpPath[k][l];
				yCoop = yCurveOpPath[k][l];
				realWorldCoor[xCoop][yCoop] = 1;
				if((xCoop >= viewXLowerBound && xCoop <= viewXUpperBound) && (yCoop >= viewYLowerBound && yCoop <= viewYUpperBound))
				{
					viewWorldCoor[xCoop][yCoop] = 1;
				}
			}
		}	
	}
	else
	{
		//Draw lines
		int k = 0;
		for(int i = 0; i < xyCurveLnIndex; i++)
		{
			for (int j = 0; j < lenCurveOpPath[i] - 1; j++)
			{
				/*xyCurLnOpPath[k][0] = xCurveOpPath[i][j];
				xyCurLnOpPath[k][1] = yCurveOpPath[i][j];*/
				if(j == 0)
				{
					xyCurLnOpPath[k][0] = xCurveOpPath[i][j];
					xyCurLnOpPath[k][1] = yCurveOpPath[i][j];
				}
				xyCurLnOpPath[k][2] = xCurveOpPath[i][j + 1];
				xyCurLnOpPath[k][3] = yCurveOpPath[i][j + 1];
				if(j != lenCurveOpPath[i] - 2)
				{
					xyCurLnOpPath[k + 1][0] = xCurveOpPath[i][j + 1];
					xyCurLnOpPath[k + 1][1] = yCurveOpPath[i][j + 1];
				}
				/*if(j == lenCurveOpPath[i] - 2 && xCurveOpPath[i][j + 1] != xyCurveOp[i][6] && yCurveOpPath[i][j + 1] != xyCurveOp[i][7])
				{
					k++;
					xyCurLnOpPath[k][0] = xCurveOpPath[i][j + 1];
					xyCurLnOpPath[k][1] = yCurveOpPath[i][j + 1];
					xyCurLnOpPath[k][2] = xyCurveOp[i][6];
					xyCurLnOpPath[k][3] = xyCurveOp[i][7];
				}*/
				/*if(xyCurLnOpPath[k][0] == 6 &&	xyCurLnOpPath[k][1] == 108 && xyCurLnOpPath[k][2] == 16 && xyCurLnOpPath[k][3] == 127)
				{
					cout<<k<<"  True"<<endl;
				}*/
				k++;
			}
		}
		/*for(int l = 0; l < k; l++)
		{
			if(xyCurLnOpPath[k][0] == 6 &&	xyCurLnOpPath[k][1] == 108 && xyCurLnOpPath[k][2] == 16 && xyCurLnOpPath[k][3] == 127)
			{
				xyCurLnOpPath[k][0] = 6;
				xyCurLnOpPath[k][1] = 108;
				xyCurLnOpPath[k][2] = 16;
				xyCurLnOpPath[k][3] = 127;
				k++;
			}
			if(xyCurLnOpPath[k][0] == 392 &&	xyCurLnOpPath[k][1] == 277 && xyCurLnOpPath[k][2] == 401 && xyCurLnOpPath[k][3] == 276)
			{
				xyCurLnOpPath[k][0] = 392;
				xyCurLnOpPath[k][1] = 277;
				xyCurLnOpPath[k][2] = 401;
				xyCurLnOpPath[k][3] = 276;
				k++;
			}
		}*/
		xyCurLnOpPathIndex = k;
		/*for(int l = 0; l < k; l++)
		{
			if(xyCurLnOpPath[l][0] == 392 &&	xyCurLnOpPath[l][1] == 277 && xyCurLnOpPath[l][2] == 401 && xyCurLnOpPath[l][3] == 276)
			{
				cout<<l<<"  Yeah!"<<endl;
			}
			cout<<l<<": ";
			cout<<xyCurLnOpPath[l][0]<<" "<<xyCurLnOpPath[l][1]<<" "<<xyCurLnOpPath[l][2]<<" "<<xyCurLnOpPath[l][3]<<endl;
		}*/
		for (int i = 0; i < xyCurLnOpPathIndex; i++)
		{
			double tempX, tempY = 0;			
			/*transform(xyOp[i][0], xyOp[i][1], 0, i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation, viewXLowerBound, viewYLowerBound, viewXUpperBound, viewYUpperBound, xLowerBound, yLowerBound, xUpperBound, yUpperBound);
			transform(xyOp[i][2], xyOp[i][3], 0, i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation, viewXLowerBound, viewYLowerBound, viewXUpperBound, viewYUpperBound, xLowerBound, yLowerBound, xUpperBound, yUpperBound);*/
			if(calXpoint(xyCurLnOpPath[i][0], xyCurLnOpPath[i][2]) == -1)
			{
				//xyOpSlope[i] = calSlope(xyOp[i][2], xyOp[i][3], xyOp[i][0], xyOp[i][1]);
				tempX = xyCurLnOpPath[i][0];
				tempY = xyCurLnOpPath[i][1];
				xyCurLnOpPath[i][0] = xyCurLnOpPath[i][2];
				xyCurLnOpPath[i][1] = xyCurLnOpPath[i][3];
				xyCurLnOpPath[i][2] = tempX;
				xyCurLnOpPath[i][3] = tempY;
			}
			//else
			//{
			//	xyOpSlope[i] = calSlope(xyOp[i][0], xyOp[i][1], xyOp[i][2], xyOp[i][3]);
			//}
			xyOpSlope[i] = calSlope(xyCurLnOpPath[i][0], xyCurLnOpPath[i][1], xyCurLnOpPath[i][2], xyCurLnOpPath[i][3]);
			//cout<< i << ": " << xyCurLnOpPath[i][0] << " " << xyCurLnOpPath[i][1] << " "  << xyCurLnOpPath[i][2] << " "  << xyCurLnOpPath[i][3] << " " << xyOpSlope[i] << endl;
		
			if(xyOpSlope[i] >= 0 && xyOpSlope[i] <= 1)
			{
				typeOpPath[i] = 0;
				//No need for modification
				bresenham(xyCurLnOpPath[i][0], xyCurLnOpPath[i][1], xyCurLnOpPath[i][2], xyCurLnOpPath[i][3], typeOpPath[i], i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation);
			}
			else if(xyOpSlope[i] >= -1 && xyOpSlope[i] < 0)
			{
				typeOpPath[i] = 1;
				//Switch y1 and y2 to make 0<=slope<=1
				bresenham(xyCurLnOpPath[i][0], -xyCurLnOpPath[i][1], xyCurLnOpPath[i][2], -xyCurLnOpPath[i][3], typeOpPath[i], i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation);
			}
			else if(xyOpSlope[i] > 1)
			{
				typeOpPath[i] = 2;
				//(y1, x1)(y2, x2) to make 0<=slope<=1
				bresenham(xyCurLnOpPath[i][1], xyCurLnOpPath[i][0], xyCurLnOpPath[i][3], xyCurLnOpPath[i][2], typeOpPath[i], i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation);
			}
			else if(xyOpSlope[i] < -1)
			{
				typeOpPath[i] = 3;
				//(y2, -x2) and (y1, -x1) to make 0<=slope<=1
				bresenham(xyCurLnOpPath[i][3], -xyCurLnOpPath[i][2], xyCurLnOpPath[i][1], -xyCurLnOpPath[i][0], typeOpPath[i], i, scalingFactor, cClockwiseRotation, xTranslation, yTranslation);
			}
			//else if(xyOpSlope[i] == 1000)
			//{
			//	for (int j = xyOp[i][1]; j <= xyOp[i][3]; j++)
			//	{
			//		writePixel(xyOp[i][0], j, 4, i, (j - xyOp[i][1]));
			//	}
			//	//bresenham(xyOp[i][0], xyOp[i][1], xyOp[i][2], xyOp[i][3], 4, i);
			//}
			//else if(xyOpSlope[i] == 1000)
			//{
			//	for (int j = xyOp[i][3]; j <= xyOp[i][1]; j++)
			//	{
			//		writePixel(xyOp[i][0], j, 5, i, (j - xyOp[i][3]));
			//	}
			//	//bresenham(xyOp[i][0], xyOp[i][1], xyOp[i][2], xyOp[i][3], 5, i);
			//}
		}
		/*for (int i = 0; i < xyCurLnOpPathIndex; i++)
		{
			cout<<i<<": ";
			cout<<xyCurLnOpPath[i][0]<<" "<<xyCurLnOpPath[i][1]<<" "<<xyCurLnOpPath[i][2]<<" "<<xyCurLnOpPath[i][3]<<endl;
		}*/
		//cout<<"******************************"<<endl;
		for (int i = 0; i < xyCurLnOpPathIndex; i++)
		{
			//cout<<i<<": "<<lenOpPath[i]<<endl;
			lenOpPathBack[i] = lenOpPath[i];
			//cout<<i<<": "<<lenOpPathBack[i]<<endl;
		}
		for(int ik = 0; ik < xyCurLnOpPathIndex; ik++)
		{
			int xCoop = 0;
			int yCoop = 0;
			int index = lenOpPath[ik];
			if(lenOpPath[ik] > lenOpPathBack[ik])
			{
				index = lenOpPath[ik];
			}
			else if(lenOpPath[ik] < lenOpPathBack[ik])
			{
				index = lenOpPathBack[ik];
			}
			/*if(ik == 115)
			{
				cout<<ik<<": "<<lenOpPath[ik]<<" "<<lenOpPathBack[ik]<<endl;
			}*/
			for (int l = 0; l < index; l++)
			{					
				xCoop = xOpPath[ik][l];
				yCoop = yOpPath[ik][l];
				realWorldCoor[xCoop][yCoop] = 1;
				if((xCoop >= viewXLowerBound && xCoop <= viewXUpperBound) && (yCoop >= viewYLowerBound && yCoop <= viewYUpperBound))
				{
					viewWorldCoor[xCoop][yCoop] = 1;
					//cout << ik << ": " << xCoop << " " << yCoop << endl;
					/*if(yCoop > yMax)
					{
						yMax = yCoop;
					}
					else if(yCoop < yMin)
					{
						yMin = yCoop;
					}*/
				}
			}
		}
	}
	/*for(int i = 0; i < xyCurveLnIndex; i++)
	{
		cout << xCurveOpPath[i][0] << " " << yCurveOpPath[i][0] << endl;
	}*/

	//scanLineFill(viewXLowerBound, viewYLowerBound, viewXUpperBound, viewYUpperBound, xyLnIndex);

	std::stringstream sl;
	sl << (xUpperBound - xLowerBound + 1);
	strLength = sl.str();
	std::stringstream sw;
	sw << (yUpperBound - yLowerBound + 1);
	strWidth = sw.str();
	
	strOutput = "/* XPM */\n";
	strOutput += "static char *quad_bw[] = {\n";
	strOutput += "/* columns rows colors chars-per-pixel */\n";
	//strOutput += "\"" + strLength + " " + strWidth + " 2 1\",\n";
	strOutput += "\" 501 501 2 1\",\n";
	strOutput += "/* pixels */\n";
	strOutput += "\"@ c #000000\",\n";
	strOutput += "\"  c #FFFFFF\",\n";

	drawPic(strOutput, countIndex, xLowerBound, xUpperBound, yLowerBound, yUpperBound, viewXLowerBound, viewYLowerBound, viewXUpperBound, viewYUpperBound);
	cout << strOutput;

	return 0;
}

string readFile(string strPath)
{
	//Read the ps file
	string strPS = "";
	char chStrPath[20];
	strcpy(chStrPath, strPath.c_str());
	char buffer[256];
	ifstream hsfile;
	hsfile.open(chStrPath);
	if(!hsfile){
        cout << "Unable to open the file";
        exit(1);  // terminate with error
	}
	while (!hsfile.eof())
	{
		//strPS = hsfile.get();
		hsfile.getline(buffer,50);
		//cout << buffer << endl;
		strPS += buffer;
		strPS += "\n";
	}
	hsfile.close();
	return strPS;
}

void bresenham(int qx, int qy, int rx, int ry, int sign, int index, double scal, int rota, int xTran, int yTran)
{
	//Bresenham's algorithm
	int length_line = 0;
	int dx, dy, D, x, y = 0;
	//cout<<qx<<","<<qy<<" "<<rx<<","<<ry<<endl;
	dx = rx - qx;
	dy = ry - qy;
	D = 2 * dy - dx;
	y = qy;
	for (x = qx; x <= rx; x++)
	{
		writePixel(x, y, sign, index, length_line);
		//cout << x << "," << y << endl;
		if(D <= 0)
		{
			D += 2 * dy;
		}
		else
		{
			D += 2 * (dy -dx) ;
			y++;
		}
		length_line++;
	}
	lenOpPath[index] = length_line;
	/*if(index == 115)
	{
		cout<<"**************"<<endl;
		cout<<lenOpPath[index]<<endl;
		cout<<"**************"<<endl;
	}*/
}

void writePixel(int x, int y, int sign, int index, int length)
{
	//Generating mapping point of the line
	switch (sign)
	{
		case 0:
			xOpPath[index][length] = x;
			yOpPath[index][length] = y;
			break;
		case 1:
			xOpPath[index][length] = x;
			yOpPath[index][length] = -y;
			break;
		case 2:	
			xOpPath[index][length] = y;
			yOpPath[index][length] = x;
			break;
		case 3:
			xOpPath[index][length] = -y;
			yOpPath[index][length] = x;
			break;
		/*case 4:
		case 5:
			xOpPath[index][length] = x;
			yOpPath[index][length] = y;
			break;*/
		default:
			xOpPath[index][length] = 0;
			yOpPath[index][length] = 0;
			cout << "Points fail" << endl;
			break;
	}
}

double calSlope(int qx, int qy, int rx, int ry)
{
	//Calculate the slope of the line
	double slope = 0.0;
	if(rx - qx != 0)
		slope = static_cast<double> (ry - qy) / (rx - qx);
	else
	{
		if(ry - qy < 0)
			slope = -1000;
		else if(ry - qy > 0)
			slope = 1000;
		else if (ry - qy == 0)
			slope = 0;
	}
	return slope;
}

int calXpoint(int qx, int rx)
{
	//Calculate the distance of points of x-coop
	if(rx - qx < 0)
	{
		return -1;
	}
	else if(rx - qx == 0)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

void drawPic(string & strOutput, int countIndex, int xLowerBound, int xUpperBound, int yLowerBound, int yUpperBound, int viewXLowerBound, int viewYLowerBound, int viewXUpperBound, int viewYUpperBound)
{
	//Building output string for xpm file
	int pExist = 0;
	int xCoop = 0;
	int yCoop = 0;
	yUpperBound = 500;
	yLowerBound = 0;
	xUpperBound = 500;
	xLowerBound = 0;
	for(int i = 0; i <= (yUpperBound - yLowerBound); i++)
	{
		//Clipping line by Sutherland-Hodgman Algorithm 
		strOutput += "\"";
		for (int j = xLowerBound; j <= xUpperBound; j++)
		{
			if(realWorldCoor[j][yUpperBound -i] == 1 && viewWorldCoor[j][yUpperBound -i] == 1 && (j >= viewXLowerBound && j <= viewXUpperBound && (yUpperBound -i) >= viewYLowerBound && (yUpperBound -i) <= viewYUpperBound))
			{
				strOutput += "@";
				//cout << ": " << j << " " << (yUpperBound -i) << endl;
			}
			else
			{
				strOutput += " ";
			}
		}
		if(i == yUpperBound)
		{
			strOutput += "\"\n};";
		}
		else
		{
			strOutput += "\",\n";
		}
	}
}

void transform(double& x, double& y, int sign, int index, double scal, int rota, int xTran, int yTran, int umin, int vmin, int umax, int vmax, int xl, int yl, int xu, int yu)
{

	//Transformation of lines
	double xScale = static_cast<double>(umax - umin) / (xu - xl);
	double yScale = static_cast<double>(vmax - vmin) / (yu - yl);
	double tempX = x * scal;
	double tempY = y * scal;;
	double cosValue = cos((rota*PI)/180.0);
	double sinValue = sin((rota*PI)/180.0);
	x = tempX * cosValue - tempY * sinValue + xTran;
	y = tempX * sinValue + tempY * cosValue + yTran;
	x = (x - xl) * xScale + umin;
	y = (y - yl) * yScale + vmin;
}

//void scanLineFill(int viewXLowerBound, int viewYLowerBound, int viewXUpperBound, int viewYUpperBound, int xyLnIndex)
//{
//	
//	for(int i = 0; i < 501; i++)
//	{
//		//Initialization
//		scanCount[i] = 0;
//		scanLineCount[i] = 0;
//		for(int j = 0; j < 501; j++)
//		{
//			scanLine[i][j] = 0;
//			scanLineSet[i][j] = -1;
//			sortedInter[i][j] = -1;
//		}
//	}
//	for(int i = yMax; i >= yMin; i--)
//	{
//		int count = 0;
//		for (int j = 0; j < xyLnIndex; j++)
//		{
//			if(xyOpSlope[j] == 0)
//			{
//				//Ignore
//			}
//			else if(xyOp[j][3] > xyOp[j][1] && i == xyOp[j][3])
//			{
//				//Ignore
//			}
//			else if(xyOp[j][3] < xyOp[j][1] && i == xyOp[j][1])
//			{
//				//Ignore
//			}
//			else if(xyOp[j][3] > xyOp[j][1] && (i < xyOp[j][3] && i >= xyOp[j][1]))
//			{
//				//Add to scan line set
//				scanLine[count][i] = j;
//				count++;
//				scanLineCount[i]++;
//			}
//			else if(xyOp[j][3] < xyOp[j][1] && (i >= xyOp[j][3] && i < xyOp[j][1]))
//			{
//				//Add to scan line set
//				scanLine[count][i] = j;
//				count++;
//				scanLineCount[i]++;
//			}
//		}
//	}
//	for(int i = yMax; i >= yMin; i--)
//	{
//		//Calculate intersections
//		int count = 0;
//		for(int j = 0; j < scanLineCount[i]; j++)
//		{
//			int lIndex = scanLine[j][i];
//			for(int k = 0; k < 2000; k++)
//			{
//				if(yOpPath[lIndex][k] == i)
//				{
//					sortedInter[count][i] = xOpPath[lIndex][k];
//					count++;
//					scanCount[i]++;
//					break;
//				}
//			}
//		}
//	}
//	for(int i = yMax; i >= yMin; i--)
//	{
//		//Sort the intersections in x
//		for(int j = 0; j < scanLineCount[i]; j++)
//		{
//			for(int k = 0; k < scanLineCount[i] - 1; k++)
//			{
//				int temp = 0;
//				if(sortedInter[k][i] > sortedInter[k + 1][i])
//				{
//					temp = sortedInter[k][i];
//					sortedInter[k][i] = sortedInter[k + 1][i];
//					sortedInter[k + 1][i] = temp;
//				}
//			}
//		}
//	}
//	/*for(int i = yMax; i >= yMin; i--)
//	{
//		cout<<i<<endl;
//		for(int j = 0; j < scanLineCount[i]; j++)
//		{
//			cout<<sortedInter[j][i]<<" ";
//		}
//		cout<<scanLineCount[i]<<endl;
//	}*/
//	for(int i = yMax; i >= yMin; i--)
//	{
//		int count = 0;
//		for (int j = 0; j < 501; j++)
//		{
//			if(realWorldCoor[j][i] == 0 && j > sortedInter[count][i] && j < sortedInter[count + 1][i])
//			{
//				realWorldCoor[j][i] = 1;
//				if((j >= viewXLowerBound && j <= viewXUpperBound) && (i >= viewYLowerBound && i <= viewYUpperBound))
//				{
//					viewWorldCoor[j][i] = 1;
//				}
//			}
//			else if(j == sortedInter[count + 1][i])
//			{
//				if(sortedInter[count][i] == sortedInter[count + 1][i])
//				{
//					j = j - 1;
//				}
//				if(scanLineCount[i] > (count + 2) || scanCount[i] > (count + 2))
//				{					
//					count = count + 2;
//				}
//				else
//				{
//					break;
//				}
//			}
//		}
//	}
//}

void CalculateBezierPoint(int &x, int &y, double t, double px1, double py1, double px2, double py2, double px3, double py3, double px4, double py4)
{
	double u = 1 - t;
	double tt = t * t;
	double uu = u * u;
	double uuu = uu * u;
	double ttt = tt * t;
	double rx = uuu * px1;
	double ry = uuu * py1;
	rx += 3 * uu * t * px2;
	ry += 3 * uu * t * py2;
	rx += 3 * u * tt * px3;
	ry += 3 * u * tt * py3;
	rx += ttt * px4;
	ry += ttt * py4;
	x = rx;
	y = ry;
 
	//float u = 1 ¨C t;
	//float tt = t*t;
	//float uu = u*u;
	//float uuu = uu * u;
	//float ttt = tt * t;
	//Vector3 p = uuu * p0; //first term
	//p += 3 * uu * t * p1; //second term
	//p += 3 * u * tt * p2; //third term
	//p += ttt * p3; //fourth term
	//return p;
}

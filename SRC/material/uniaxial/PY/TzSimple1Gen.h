#include <iostream>
#include <fstream>
#include <cmath>
#include <ID.h> 

class TzSimple1Gen
{

	// Variables used for reading input files:
	int NumNodes, NumTzEle, NumPileEle, NPile, NumLayer, NumMtLoadSp, NumLoad, 
		NumSp, NumMt, NumMat;
	double p, zground, TULT, ZULT, ru, ca, depth, stress, delta, b, Sa;
	int *NodeNum;								// Arrays for Nodes File
	double *Nodey, *Nodex;
	int *TzEleNum, *TzNode1, *TzNode2, *TzMat, *TzDir;	// Arrays for Py Elements File
	int *PileEleNum, *PileNode1, *PileNode2;			// Arrays for Pile Elements File
	int *tzType;
	double *gamma_t, *gamma_b, *z_t, *z_b, *p_t, *p_b, *c_t, *c_b, *ca_t, *ca_b, *delta_t, *delta_b,
		*zLoad_t, *zLoad_b, *load_val_t, *load_val_b, *zSp_t, *zSp_b, *sp_val_t,
		*sp_val_b, *zMt_t, *zMt_b, *mt_val_t, *mt_val_b, tribcoord[2], *Sa_b, *Sa_t, *ru_t, *ru_b,
		*tult_t, *tult_b, *zult_t, *zult_b;
	char **MatType, *PatternInfo;
	

	// Member functions for reading input files:
	void GetNodes(char *file);
	void GetTzElements(char *file);
	void GetPileElements(char *file);
	void GetSoilProperties(char *file);
	double GetTult(char *type);
	double GetZ50(char *type);
	double GetMt(double *vx, double *vy, double x, int length);
	void GetTributaryCoordsTz(int nodenum1);
	void GetTributaryCoordsPile(int nodenum1);
	int NumRows(char *file, char *begin);

	// Member functions for generating output:
	void GetTzSimple1(char *file1, char *file2, char *file3, char *file4, char *file5);
	void GetPattern(char *file6);

	// Member functions for calculating tult:
	double GetVStress(double z);
	double linterp(double x1, double x2, double y1, double y2, double x3);
	
public:

	void WriteTzSimple1(char *file1, char *file2, char *file3, char *file4, char *file5);
	void WriteTzSimple1(char *file1, char *file2, char *file3, char *file4, char *file5, char *file6);

	TzSimple1Gen();
	~TzSimple1Gen();
};

#ifndef MaterialResponse_h
#define MaterialResponse_h

#include <Response.h>
#include <Information.h>

class Material;

class ID;
class Vector;
class Matrix;
class Tensor;

class MaterialResponse : public Response
{
public:
	MaterialResponse(Material *mat, int id);
	MaterialResponse(Material *mat, int id, int val);
	MaterialResponse(Material *mat, int id, double val);
	MaterialResponse(Material *mat, int id, const ID &val);
	MaterialResponse(Material *mat, int id, const Vector &val);
	MaterialResponse(Material *mat, int id, const Matrix &val);
	MaterialResponse(Material *mat, int ID, const Tensor &val);
	~MaterialResponse();

	int getResponse(void);
	void Print(ostream &s, int flag = 0);

private:
	Material *theMaterial;
	int responseID;
	Information matInfo;
};

#endif

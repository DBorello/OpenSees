#include <MaterialResponse.h>
#include <Material.h>

MaterialResponse::MaterialResponse(Material *mat, int id):
Response(), theMaterial(mat), responseID(id), matInfo()
{

}

MaterialResponse::MaterialResponse(Material *mat, int id, int val):
Response(), theMaterial(mat), responseID(id), matInfo(val)
{

}

MaterialResponse::MaterialResponse(Material *mat, int id, double val):
Response(), theMaterial(mat), responseID(id), matInfo(val)
{

}

MaterialResponse::MaterialResponse(Material *mat, int id, const ID &val):
Response(), theMaterial(mat), responseID(id), matInfo(val)
{

}

MaterialResponse::MaterialResponse(Material *mat, int id, const Vector &val):
Response(), theMaterial(mat), responseID(id), matInfo(val)
{

}

MaterialResponse::MaterialResponse(Material *mat, int id, const Matrix &val):
Response(), theMaterial(mat), responseID(id), matInfo(val)
{

}

MaterialResponse::MaterialResponse(Material *mat, int id, const Tensor &val):
Response(), theMaterial(mat), responseID(id), matInfo(val)
{

}

MaterialResponse::~MaterialResponse()
{

}

int
MaterialResponse::getResponse(void)
{
	return theMaterial->getResponse(responseID, matInfo);
}

void
MaterialResponse::Print(ostream &s, int flag)
{
	matInfo.Print(s, flag);
}


#include <ElementResponse.h>
#include <Element.h>

ElementResponse::ElementResponse(Element *ele, int id):
Response(), theElement(ele), responseID(id), eleInfo()
{

}

ElementResponse::ElementResponse(Element *ele, int id, int val):
Response(), theElement(ele), responseID(id), eleInfo(val)
{

}

ElementResponse::ElementResponse(Element *ele, int id, double val):
Response(), theElement(ele), responseID(id), eleInfo(val)
{

}

ElementResponse::ElementResponse(Element *ele, int id, const ID &val):
Response(), theElement(ele), responseID(id), eleInfo(val)
{

}

ElementResponse::ElementResponse(Element *ele, int id, const Vector &val):
Response(), theElement(ele), responseID(id), eleInfo(val)
{

}

ElementResponse::ElementResponse(Element *ele, int id, const Matrix &val):
Response(), theElement(ele), responseID(id), eleInfo(val)
{

}

ElementResponse::ElementResponse(Element *ele, int id, const Tensor &val):
Response(), theElement(ele), responseID(id), eleInfo(val)
{

}

ElementResponse::~ElementResponse()
{

}

int
ElementResponse::getResponse(void)
{
	return theElement->getResponse(responseID, eleInfo);
}

void
ElementResponse::Print(ostream &s, int flag)
{
	eleInfo.Print(s, flag);
}

#ifndef ElementResponse_h
#define ElementResponse_h

#include <Response.h>
#include <Information.h>

class Element;

class ID;
class Vector;
class Matrix;
class Tensor;

class ElementResponse : public Response
{
public:
	ElementResponse(Element *ele, int id);
	ElementResponse(Element *ele, int id, int val);
	ElementResponse(Element *ele, int id, double val);
	ElementResponse(Element *ele, int id, const ID &val);
	ElementResponse(Element *ele, int id, const Vector &val);
	ElementResponse(Element *ele, int id, const Matrix &val);
	ElementResponse(Element *ele, int id, const Tensor &val);
	~ElementResponse();

	int getResponse(void);
	void Print(ostream &s, int flag = 0);

private:
	Element *theElement;
	int responseID;
	Information eleInfo;
};

#endif

#ifndef FiberResponse_h
#define FiberResponse_h

#include <Response.h>
#include <Information.h>

class Fiber;

class ID;
class Vector;
class Matrix;
class Tensor;

class FiberResponse : public Response
{
public:
	FiberResponse(Fiber *fib, int id);
	FiberResponse(Fiber *fib, int id, int val);
	FiberResponse(Fiber *fib, int id, double val);
	FiberResponse(Fiber *fib, int id, const ID &val);
	FiberResponse(Fiber *fib, int id, const Vector &val);
	FiberResponse(Fiber *fib, int id, const Matrix &val);
	FiberResponse(Fiber *fib, int id, const Tensor &val);
	~FiberResponse();

	int getResponse(void);
	void Print(ostream &s, int flag = 0);

private:
	Fiber *theFiber;
	int responseID;
	Information fibInfo;
};

#endif

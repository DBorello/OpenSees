#include <FiberResponse.h>
#include <Fiber.h>

FiberResponse::FiberResponse(Fiber *fib, int id):
Response(), theFiber(fib), responseID(id), fibInfo()
{

}

FiberResponse::FiberResponse(Fiber *fib, int id, int val):
Response(), theFiber(fib), responseID(id), fibInfo(val)
{

}

FiberResponse::FiberResponse(Fiber *fib, int id, double val):
Response(), theFiber(fib), responseID(id), fibInfo(val)
{

}

FiberResponse::FiberResponse(Fiber *fib, int id, const ID &val):
Response(), theFiber(fib), responseID(id), fibInfo(val)
{

}

FiberResponse::FiberResponse(Fiber *fib, int id, const Vector &val):
Response(), theFiber(fib), responseID(id), fibInfo(val)
{

}

FiberResponse::FiberResponse(Fiber *fib, int id, const Matrix &val):
Response(), theFiber(fib), responseID(id), fibInfo(val)
{

}

FiberResponse::FiberResponse(Fiber *fib, int id, const Tensor &val):
Response(), theFiber(fib), responseID(id), fibInfo(val)
{

}

FiberResponse::~FiberResponse()
{

}

int
FiberResponse::getResponse(void)
{
	return theFiber->getResponse(responseID, fibInfo);
}

void
FiberResponse::Print(ostream &s, int flag)
{
	fibInfo.Print(s, flag);
}

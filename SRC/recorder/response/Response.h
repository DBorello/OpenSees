class ostream;

class Response
{
public:
	Response();
	virtual ~Response();

	virtual int getResponse(void) = 0;
	virtual void Print(ostream &s, int flag = 0) = 0;
};
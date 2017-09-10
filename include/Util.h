#ifndef _UTIL_H_
#define _UTIL_H_

#include <cstdio>

/************************************************************************/
/*                   CTimer                                             */
/************************************************************************/
#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

template <class T>
class _CTimer
{
public:
	_CTimer();
	~_CTimer(){};

	void tic();	// time in clock
	double toc(const char* msg = NULL); // time out clock
private:
	T inT, f;
};

template <class T>
_CTimer<T>::_CTimer()
{
#ifdef WIN32
	QueryPerformanceFrequency(&f);
#endif
	tic();
}

#ifdef WIN32
typedef _CTimer<LARGE_INTEGER> CTimer; // only CTimer makes sense
#else
typedef _CTimer<struct timeval> CTimer;
#endif

template <class T>
void _CTimer<T>::tic()
{
#ifdef WIN32
	QueryPerformanceCounter(&inT);
#else
	gettimeofday(&inT,NULL);
#endif
}

template <class T>
double _CTimer<T>::toc(const char* msg)
{
#ifdef WIN32
	LARGE_INTEGER outT;
	QueryPerformanceCounter(&outT);

	double dt = (double)(outT.QuadPart-inT.QuadPart)/f.QuadPart;
#else
	struct timeval outT;
	gettimeofday(&outT,NULL);
	double dt = 1000000 * (outT.tv_sec - inT.tv_sec) + outT.tv_usec - inT.tv_usec;
	dt /= 1000000;
#endif
	if(msg){
		printf("%s %f [s]\n", msg, dt);
	}
	
	tic();
	return dt;
}

/************************************************************************/
/*                      Gray2ColorTable                                 */
/************************************************************************/
template<class T>
class _CColorTable
{
public:
	_CColorTable();
	~_CColorTable(){};

	// function to access the member variables
	inline T* operator [] (int index) { return m_colorTbl[index]; };

private:
	T m_colorTbl[256][3];
};

typedef _CColorTable<unsigned char> CColorTable;

template<class T>
_CColorTable<T>::_CColorTable()
{
	// BGR table
static unsigned char gray2ColorTb[256][3] = {
{143,  0,  0},{147,  0,  0},{151,  0,  0},{155,  0,  0},{159,  0,  0},{163,  0,  0},{167,  0,  0},{171,  0,  0},{175,  0,  0},{179,  0,  0},{183,  0,  0},{187,  0,  0},{191,  0,  0},{195,  0,  0},{199,  0,  0},{203,  0,  0},{206,  0,  0},{210,  0,  0},{214,  0,  0},{218,  0,  0},{222,  0,  0},{226,  0,  0},{230,  0,  0},{234,  0,  0},{238,  0,  0},{242,  0,  0},{246,  0,  0},{250,  0,  0},{254,  0,  0},{255,  3,  0},{255,  7,  0},{255, 11,  0},
{255, 14,  0},{255, 18,  0},{255, 22,  0},{255, 26,  0},{255, 30,  0},{255, 34,  0},{255, 38,  0},{255, 42,  0},{255, 46,  0},{255, 50,  0},{255, 54,  0},{255, 58,  0},{255, 62,  0},{255, 66,  0},{255, 70,  0},{255, 74,  0},{255, 77,  0},{255, 81,  0},{255, 85,  0},{255, 89,  0},{255, 93,  0},{255, 97,  0},{255,101,  0},{255,105,  0},{255,109,  0},{255,113,  0},{255,117,  0},{255,121,  0},{255,125,  0},{255,129,  0},{255,133,  0},{255,137,  0},
{255,140,  0},{255,144,  0},{255,148,  0},{255,152,  0},{255,156,  0},{255,160,  0},{255,164,  0},{255,168,  0},{255,172,  0},{255,176,  0},{255,180,  0},{255,184,  0},{255,188,  0},{255,192,  0},{255,196,  0},{255,200,  0},{255,203,  0},{255,207,  0},{255,211,  0},{255,215,  0},{255,219,  0},{255,223,  0},{255,227,  0},{255,231,  0},{255,235,  0},{255,239,  0},{255,243,  0},{255,247,  0},{255,251,  0},{255,255,  0},{251,255,  4},{247,255,  8},
{244,255, 11},{240,255, 15},{236,255, 19},{232,255, 23},{228,255, 27},{224,255, 31},{220,255, 35},{216,255, 39},{212,255, 43},{208,255, 47},{204,255, 51},{200,255, 55},{196,255, 59},{192,255, 63},{188,255, 67},{184,255, 71},{181,255, 74},{177,255, 78},{173,255, 82},{169,255, 86},{165,255, 90},{161,255, 94},{157,255, 98},{153,255,102},{149,255,106},{145,255,110},{141,255,114},{137,255,118},{133,255,122},{129,255,126},{125,255,130},{122,255,133},
{118,255,137},{114,255,141},{110,255,145},{106,255,149},{102,255,153},{ 98,255,157},{ 94,255,161},{ 90,255,165},{ 86,255,169},{ 82,255,173},{ 78,255,177},{ 74,255,181},{ 70,255,185},{ 66,255,189},{ 62,255,193},{ 59,255,196},{ 55,255,200},{ 51,255,204},{ 47,255,208},{ 43,255,212},{ 39,255,216},{ 35,255,220},{ 31,255,224},{ 27,255,228},{ 23,255,232},{ 19,255,236},{ 15,255,240},{ 11,255,244},{  7,255,248},{  3,255,252},{  0,254,255},{  0,251,255},
{  0,247,255},{  0,243,255},{  0,239,255},{  0,235,255},{  0,231,255},{  0,227,255},{  0,223,255},{  0,219,255},{  0,215,255},{  0,211,255},{  0,207,255},{  0,203,255},{  0,199,255},{  0,195,255},{  0,191,255},{  0,188,255},{  0,184,255},{  0,180,255},{  0,176,255},{  0,172,255},{  0,168,255},{  0,164,255},{  0,160,255},{  0,156,255},{  0,152,255},{  0,148,255},{  0,144,255},{  0,140,255},{  0,136,255},{  0,132,255},{  0,128,255},{  0,125,255},
{  0,121,255},{  0,117,255},{  0,113,255},{  0,109,255},{  0,105,255},{  0,101,255},{  0, 97,255},{  0, 93,255},{  0, 89,255},{  0, 85,255},{  0, 81,255},{  0, 77,255},{  0, 73,255},{  0, 69,255},{  0, 65,255},{  0, 62,255},{  0, 58,255},{  0, 54,255},{  0, 50,255},{  0, 46,255},{  0, 42,255},{  0, 38,255},{  0, 34,255},{  0, 30,255},{  0, 26,255},{  0, 22,255},{  0, 18,255},{  0, 14,255},{  0, 10,255},{  0,  6,255},{  0,  2,255},{  0,  0,254},
{  0,  0,250},{  0,  0,246},{  0,  0,242},{  0,  0,238},{  0,  0,234},{  0,  0,230},{  0,  0,226},{  0,  0,222},{  0,  0,218},{  0,  0,214},{  0,  0,210},{  0,  0,206},{  0,  0,202},{  0,  0,198},{  0,  0,194},{  0,  0,191},{  0,  0,187},{  0,  0,183},{  0,  0,179},{  0,  0,175},{  0,  0,171},{  0,  0,167},{  0,  0,163},{  0,  0,159},{  0,  0,155},{  0,  0,151},{  0,  0,147},{  0,  0,143},{  0,  0,139},{  0,  0,135},{  0,  0,131},{  0,  0,128}
};
	for (int i = 0; i < 256; i++){
		memcpy(m_colorTbl[i], gray2ColorTb[i], 3);
	}
}

#endif // _UTIL_H_


#ifndef IMAGE_H
#define IMAGE_H

#include <Windows.h>
#include "FreeImage.h"
#include "MATH.h"

namespace RENDERER
{
	using namespace ZMATH;

	//-------------------- Í¼ÏñÀà
	class Image
	{
	private:
		int _width;
		int _height;
		uint * _pixel;
		int _wrapType;

	public:
		Image(int w, int h, void * data);
		~Image();
		
		int width() const
		{
			return _width;
		}
		int height() const
		{
			return _height;
		}
		void setWrapType(int type)
		{
			_wrapType = type;
		}
		Rgba pixelAt(int x, int y) const
		{
			return Rgba(_pixel[y * _width + x]);
		}
		Rgba pixelUV(float u, float v)
		{
			float x = u * _width;
			float y = v * _height;
			if (_wrapType == 0)
			{
				return pixelAt((unsigned)(x) % _width, (unsigned)(y) % _height);
			}
			else
			{
				if (x >= _width)
				{
					x = _width - 1;
				}
				if (y >= _height)
				{
					y = _height - 1;
				}
				return pixelAt(x, y);
			}
		}

		static Image * loadFromFile(const char * fileName);
	};

}

#endif // !IMAGE_H
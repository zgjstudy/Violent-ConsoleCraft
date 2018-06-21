#include "Image.h"

namespace RENDERER
{
	Image::Image(int w, int h, void * data)
	{
		if (w == 0 || h == 0 || data == 0)
		{
			_width = 0;
			_height = 0;
			_pixel = 0;
		}
		else
		{
			_width = w;
			_height = h;
			_pixel = new uint[w * h];
			memcpy(_pixel, data, w * h * sizeof(uint));
		}
	}
	Image::~Image()
	{
		delete[] _pixel;
	}
	Image * RENDERER::Image::loadFromFile(const char * fileName)
	{
		//获取图片格式
		FREE_IMAGE_FORMAT fifmt = FreeImage_GetFileType(fileName, 0);
		if (fifmt == FIF_UNKNOWN)
		{
			return 0;
		}
		//加载图片
		FIBITMAP * dib = FreeImage_Load(fifmt, fileName, 0);

		FREE_IMAGE_COLOR_TYPE type = FreeImage_GetColorType(dib);

		//获取数据指针
		FIBITMAP * temp = dib;
		dib = FreeImage_ConvertTo32Bits(dib);
		FreeImage_Unload(temp);

		BYTE * pixels = (BYTE*)FreeImage_GetBits(dib);
		int width = FreeImage_GetWidth(dib);
		int height = FreeImage_GetHeight(dib);

		//图片的上下反转
		int pitch = width * 4;
		BYTE * trow = new BYTE[pitch];
		for (int j = 0; j < height / 2; j++)
		{
			memcpy(trow, pixels + j * pitch, pitch);

			memcpy(pixels + j * pitch, pixels + (height - j - 1) * pitch, pitch);

			memcpy(pixels + (height - j - 1) * pitch, trow, pitch);
		}
		delete[]trow;

		Image * image = new Image(width, height, pixels);
		FreeImage_Unload(dib);

		return image;
	}
}
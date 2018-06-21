#ifndef RASTER3D_H
#define RASTER3D_H

#include <Windows.h>
#include <iostream>
#include <tchar.h>
#include <stdio.h>
#include "MATH.h"
#include "Image.h"
#include "RasterDependence.h"

namespace RENDERER
{
	using namespace ZMATH;

	class Span
	{
	public:
		int		_xStart;
		int		_xEnd;
		Rgba	_colorStart;
		Rgba	_colorEnd;

		float2	_uvStart;
		float2	_uvEnd;

		int		_y;

		Span(int xStart, int xEnd, int y, Rgba colorStart, Rgba colorEnd, float2 uvStart, float2 uvEnd);
	};

	class Edge
	{
	public:
		int		_x1;
		int		_y1;
		float2	_uv1;
		Rgba	_color1;

		int		_x2;
		int		_y2;
		float2	_uv2;
		Rgba	_color2;

		//1号点在应2号点上方
		Edge(int x1, int y1, Rgba color1, float2 uv1, int x2, int y2, Rgba color2, float2 uv2);
	};

	//-------------------- 渲染引擎
	class Raster3D
	{
	private:
		struct Vertex
		{
			int2 p0;
			float2 uv0;
			Rgba c0;

			int2 p1;
			float2 uv1;
			Rgba c1;

			int2 p2;
			float2 uv2;
			Rgba c2;
		};

	public:
		uint * _buffer;
		int				_height;
		int				_width;
		Rgba			_color;
		Image *			_texture;
		float *			_Zbuffer;

		matrix4			_matrixModel;		//模型矩阵
		matrix4			_matrixView;		//观察矩阵
		matrix4			_matrixProj;		//投影矩阵
		matrix4			_matrixProjView;
		float2			_viewPort;			//观察位置（宽高）
		Frustum			_frust;				//截椎体

		DataElementDesc	_positionPointer;	//顶点数据
		DataElementDesc	_colorPointer;		//颜色数据
		DataElementDesc	_uvPointer;			//uv数据

		DataElementDesc	_defaultColorPointer;
		DataElementDesc	_defaultUVPointer;
		Rgba			_defaultColorArray[3];
		float2			_detaultUVArray[3];

		bool _ignoreOutside = 0;

		Raster3D(int w, int h, void * buffer);

		~Raster3D();

		void clear();



		//-------------------- 状态位操作
		void NOTignoreOutside()
		{
			_ignoreOutside = 0;
		}
		void ignoreOutside()
		{
			_ignoreOutside = 1;
		}

		//-------------------- 观察矩阵操作
		void loadViewMatrix(const matrix4 & mat);

		void resetViewMatrix();

		//生成观察矩阵
		void lookat(float3 const & eye, float3 const & center, float3 const & up);

		void setView(const matrix4 & viewMat);

		void setViewPort(int x, int y, int w, int h);

		//-------------------- 投影矩阵操作
		void loadProjMatrix(const matrix4 & mat);

		void resetProjMatrix();

		//生成投影矩阵
		void setPerspective(float fovy, float aspect, float zNear, float zFar);


		//-------------------- 模型矩阵操作
		void loadMatrix(const matrix4 & mat);

		void resetMatrix();

		//-------------------- 矩阵操作
		void setVertexPointer(int size, DATATYPE type, int stride, const void* data);

		void setColorPointer(int size, DATATYPE type, int stride, const void* data);

		void setTextureCoordPointer(int size, DATATYPE type, int stride, const void* data);

		//-------------------- 接口操作
		void drawArrays(DRAWMODE dm, int ptstart, int ptcount);

		void drawImage(int startX, int startY, const Image* image, float alpha = 1);

		void bindTexture(Image * image);

		void drawLine2D(int2 pt1, int2 pt2, Rgba color1, Rgba color2);

		void initializeZbuffer();

	private:
		//-------------------- 管线转化
		float3 piplineTransform(float3 pos);

		//-------------------- 2D三角形绘制
		void drawEdge(const Edge& e1, const Edge& e2, Image * image);

		void drawSpan(const Span& span, Image * image);

		void drawTriangle(Edge edges[3]);

	public:
		//-------------------- 像素操作
		inline Rgba getPixel(uint x, uint y)
		{
			return Rgba(_buffer[y * _width + x]);
		}

		inline void setPixel(uint x, uint y, Rgba color)
		{
			if (x >= _width || y >= _height)
			{
				return;
			}
			_buffer[y * _width + x] = color._color;
		}

		inline void setPixelEx(uint x, uint y, Rgba color)
		{
			_buffer[y * _width + x] = color._color;
		}

		//-------------------- 深度缓冲区操作
		inline void setZat(uint x, uint y, float Z)
		{
			if (x >= _width || y >= _height)
			{
				return;
			}
			_Zbuffer[y * _width + x] = Z;
		}

		inline float getZat(uint x, uint y)
		{
			if (x >= _width || y >= _height)
			{
				return 0;
			}
			return _Zbuffer[y * _width + x];
		}

		inline bool DepthTest(uint x, uint y, float testZ)
		{
			if (x >= _width || y >= _height)
			{
				return true;
			}
			float & targetZ = _Zbuffer[y * _width + x];
			if (testZ <= targetZ && testZ  > 0.0f)
			{
				targetZ = testZ;
				return TRUE;
			}
			else
			{
				return false;
			}
		}
	};

}

#endif // !RASTER3D_H
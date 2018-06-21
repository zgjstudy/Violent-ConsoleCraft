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

		//1�ŵ���Ӧ2�ŵ��Ϸ�
		Edge(int x1, int y1, Rgba color1, float2 uv1, int x2, int y2, Rgba color2, float2 uv2);
	};

	//-------------------- ��Ⱦ����
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

		matrix4			_matrixModel;		//ģ�;���
		matrix4			_matrixView;		//�۲����
		matrix4			_matrixProj;		//ͶӰ����
		matrix4			_matrixProjView;
		float2			_viewPort;			//�۲�λ�ã���ߣ�
		Frustum			_frust;				//��׵��

		DataElementDesc	_positionPointer;	//��������
		DataElementDesc	_colorPointer;		//��ɫ����
		DataElementDesc	_uvPointer;			//uv����

		DataElementDesc	_defaultColorPointer;
		DataElementDesc	_defaultUVPointer;
		Rgba			_defaultColorArray[3];
		float2			_detaultUVArray[3];

		bool _ignoreOutside = 0;

		Raster3D(int w, int h, void * buffer);

		~Raster3D();

		void clear();



		//-------------------- ״̬λ����
		void NOTignoreOutside()
		{
			_ignoreOutside = 0;
		}
		void ignoreOutside()
		{
			_ignoreOutside = 1;
		}

		//-------------------- �۲�������
		void loadViewMatrix(const matrix4 & mat);

		void resetViewMatrix();

		//���ɹ۲����
		void lookat(float3 const & eye, float3 const & center, float3 const & up);

		void setView(const matrix4 & viewMat);

		void setViewPort(int x, int y, int w, int h);

		//-------------------- ͶӰ�������
		void loadProjMatrix(const matrix4 & mat);

		void resetProjMatrix();

		//����ͶӰ����
		void setPerspective(float fovy, float aspect, float zNear, float zFar);


		//-------------------- ģ�;������
		void loadMatrix(const matrix4 & mat);

		void resetMatrix();

		//-------------------- �������
		void setVertexPointer(int size, DATATYPE type, int stride, const void* data);

		void setColorPointer(int size, DATATYPE type, int stride, const void* data);

		void setTextureCoordPointer(int size, DATATYPE type, int stride, const void* data);

		//-------------------- �ӿڲ���
		void drawArrays(DRAWMODE dm, int ptstart, int ptcount);

		void drawImage(int startX, int startY, const Image* image, float alpha = 1);

		void bindTexture(Image * image);

		void drawLine2D(int2 pt1, int2 pt2, Rgba color1, Rgba color2);

		void initializeZbuffer();

	private:
		//-------------------- ����ת��
		float3 piplineTransform(float3 pos);

		//-------------------- 2D�����λ���
		void drawEdge(const Edge& e1, const Edge& e2, Image * image);

		void drawSpan(const Span& span, Image * image);

		void drawTriangle(Edge edges[3]);

	public:
		//-------------------- ���ز���
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

		//-------------------- ��Ȼ���������
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
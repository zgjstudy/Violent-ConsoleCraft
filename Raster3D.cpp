#include "Raster3D.h"

namespace RENDERER
{

	RENDERER::Raster3D::Raster3D(int w, int h, void * buffer)
		:_height(h), _width(w)
	{
		_texture = 0;
		_buffer = (uint*)buffer;
		memset(&_positionPointer, 0, sizeof(_positionPointer));
		memset(&_colorPointer, 0, sizeof(_colorPointer));
		memset(&_uvPointer, 0, sizeof(_uvPointer));

		_defaultColorPointer._size		=	4;
		_defaultColorPointer._type		=	DT_BYTE;
		_defaultColorPointer._stride	=	sizeof(Rgba);
		_defaultColorPointer._data		=	_defaultColorArray;

		_defaultUVPointer._size		=	 2;
		_defaultUVPointer._type		=	 DT_FLOAT;
		_defaultUVPointer._stride	=	 sizeof(float2);
		_defaultUVPointer._data		=	 _detaultUVArray;
			
		_matrixModel	=	 matrix4(1);
		_matrixView		=	 matrix4(1);
		_matrixProj		=	 matrix4(1);

		_Zbuffer = new float[w * h];
		initializeZbuffer();
	}

	Raster3D::~Raster3D()
	{
		delete[] _Zbuffer;
	}

	void RENDERER::Raster3D::clear()
	{
		memset(_buffer, 0, _height * _width * sizeof(Rgba));
	}

	void Raster3D::loadViewMatrix(const matrix4 & mat)
	{
		_matrixView = mat;
	}

	void Raster3D::resetViewMatrix()
	{
		_matrixView = matrix4(1);
	}

	void Raster3D::lookat(float3 const & eye, float3 const & center, float3 const & up)
	{
		_matrixView = lookAt(eye, center, up);
	}

	void Raster3D::setView(const matrix4 & viewMat)
	{
		_matrixView = viewMat;
	}

	void Raster3D::setViewPort(int x, int y, int w, int h)
	{
		_viewPort.x = w;
		_viewPort.y = h;
	}

	void Raster3D::loadProjMatrix(const matrix4 & mat)
	{
		_matrixProj = mat;
	}

	void Raster3D::resetProjMatrix()
	{
		_matrixProj = matrix4(1);
	}

	void Raster3D::setPerspective(float fovy, float aspect, float zNear, float zFar)
	{
		_matrixProj = perspective<float>(fovy, aspect, zNear, zFar);
	}

	void Raster3D::loadMatrix(const matrix4 & mat)
	{
		_matrixModel = mat;
	}

	void Raster3D::resetMatrix()
	{
		_matrixModel = RENDERER::matrix4(1);
	}

	void Raster3D::setVertexPointer(int size, DATATYPE type, int stride, const void * data)
	{
		_positionPointer._size = size;
		_positionPointer._type = type;
		_positionPointer._stride = stride;
		_positionPointer._data = data;
	}

	void Raster3D::setColorPointer(int size, DATATYPE type, int stride, const void * data)
	{
		_colorPointer._size = size;
		_colorPointer._type = type;
		_colorPointer._stride = stride;
		_colorPointer._data = data;
	}

	void Raster3D::setTextureCoordPointer(int size, DATATYPE type, int stride, const void * data)
	{
		_uvPointer._size = size;
		_uvPointer._type = type;
		_uvPointer._stride = stride;
		_uvPointer._data = data;
	}

	void Raster3D::drawArrays(DRAWMODE dm, int start, int count)
	{
		if (_positionPointer._data == 0)
		{
			return;
		}

		DataElementDesc colorPointerdesc = _colorPointer;
		DataElementDesc uvPointerdesc = _uvPointer;

		if (colorPointerdesc._data == 0)
		{
			colorPointerdesc = _defaultColorPointer;
		}
		if (uvPointerdesc._data == 0)
		{
			uvPointerdesc = _defaultUVPointer;
		}
		char * posData = (char*)_positionPointer._data;
		char * cData = (char*)colorPointerdesc._data;
		char * uvData = (char*)uvPointerdesc._data;

		//构建截锥体
		_matrixProjView = _matrixProj * _matrixView;

		matrix4 matPVT = _matrixProjView.transpose();

		_frust.loadFrustum(matPVT);

		float * fData;

		switch (dm)
		{
		case 0: //DM_POINTS
		{
			for (int i = start; i < start + count; ++i)
			{
				fData = (float*)posData;
				float3 p01(fData[0], fData[1], fData[2]);
				posData += _positionPointer._stride;

				p01 = p01 * _matrixModel;

				//转化为屏幕坐标
				if (!_ignoreOutside || _frust.pointInFrustum(p01))
				{
					p01 = piplineTransform(p01);
					if (DepthTest((uint)p01.x, (uint)p01.y, p01.z))
					{
						int2 p0((uint)p01.x, (uint)p01.y);
						Rgba * pColor = (Rgba*)cData;
						Rgba c0(*pColor);
						cData += _colorPointer._stride;

						setPixel(p0.x, p0.y, c0);

						if (_colorPointer._data == 0)
						{
							cData = (char*)colorPointerdesc._data;
						}
					}
				}
			}
			break;
		}
		case 1: //DM_LINES
		{
			for (int i = start; i < start + count; i += 2)
			{
				fData = (float*)posData;
				float3 p01(fData[0], fData[1], fData[2]);
				posData += _positionPointer._stride;
				fData = (float*)(posData);
				float3 p11(fData[0], fData[1], fData[2]);
				posData += _positionPointer._stride;

				p01 = p01 * _matrixModel;
				p11 = p11 * _matrixModel;

				//转化为屏幕坐标

				p01 = piplineTransform(p01);
				p11 = piplineTransform(p11);
				//if (DepthTest((uint)p01.x, (uint)p01.y, p01.z) && DepthTest((uint)p11.x, (uint)p11.y, p11.z)) //需保证线间不互相穿过
				{
					int2 p0(p01.x, p01.y);
					int2 p1(p11.x, p11.y);

					Rgba * pColor = (Rgba*)cData;
					Rgba c0(*pColor);
					cData += _colorPointer._stride;
					Rgba c1(*pColor);
					cData += _colorPointer._stride;

					drawLine2D(p0, p1, c0, c1);
				}

				if (_colorPointer._data == 0)
				{
					cData = (char*)colorPointerdesc._data;
				}

			}
			break;
		}
		case 2: //DM_LINE_LOOP
		{
			for (int i = start + 1; i < start + count; ++i)
			{
				fData = (float*)posData;
				float3 p01(fData[0], fData[1], fData[2]);

				posData += _positionPointer._stride;
				fData = (float*)(posData);
				float3 p11(fData[0], fData[1], fData[2]);

				p01 = p01 * _matrixModel;
				p11 = p11 * _matrixModel;

				//转化为屏幕坐标

				p01 = piplineTransform(p01);
				p11 = piplineTransform(p11);

				//if (DepthTest((uint)p01.x, (uint)p01.y, p01.z) && DepthTest((uint)p11.x, (uint)p11.y, p11.z)) //需保证线间不互相穿过
				{
					int2 p0(p01.x, p01.y);
					int2 p1(p11.x, p11.y);

					Rgba * pColor = (Rgba*)cData;
					Rgba c0(*pColor);
					cData += _colorPointer._stride;
					Rgba c1(*pColor);
					cData += _colorPointer._stride;

					drawLine2D(p0, p1, c0, c1);
				}

				if (_colorPointer._data == 0)
				{
					cData = (char*)colorPointerdesc._data;
				}

			}

			fData = (float*)posData;
			float3 p01(fData[0], fData[1], fData[2]);

			posData = (char*)_positionPointer._data + _positionPointer._stride * start;
			fData = (float*)(posData);
			float3 p11(fData[0], fData[1], fData[2]);

			p01 = p01 * _matrixModel;
			p11 = p11 * _matrixModel;

			//转化为屏幕坐标

			p01 = piplineTransform(p01);
			p11 = piplineTransform(p11);

			//if (DepthTest((uint)p01.x, (uint)p01.y, p01.z) && DepthTest((uint)p11.x, (uint)p11.y, p11.z)) //需保证线间不互相穿过
			{
				int2 p0(p01.x, p01.y);
				int2 p1(p11.x, p11.y);

				Rgba * pColor = (Rgba*)cData;
				Rgba c0(*pColor);
				cData += _colorPointer._stride;
				Rgba c1(*pColor);
				cData += _colorPointer._stride;

				drawLine2D(p0, p1, c0, c1);
			}
			if (_colorPointer._data == 0)
			{
				cData = (char*)colorPointerdesc._data;
			}
			break;
		}
		case 3: //DM_LINE_STRIP
		{
			for (int i = start + 1; i < start + count; ++i)
			{
				fData = (float*)posData;
				float3 p01(fData[0], fData[1], fData[2]);

				posData += _positionPointer._stride;
				fData = (float*)(posData);
				float3 p11(fData[0], fData[1], fData[2]);

				p01 = p01 * _matrixModel;
				p11 = p11 * _matrixModel;

				//转化为屏幕坐标

				p01 = piplineTransform(p01);
				p11 = piplineTransform(p11);
				//if (DepthTest((uint)p01.x, (uint)p01.y, p01.z) && DepthTest((uint)p11.x, (uint)p11.y, p11.z)) //需保证线间不互相穿过
				{
					int2 p0(p01.x, p01.y);
					int2 p1(p11.x, p11.y);

					Rgba * pColor = (Rgba*)cData;
					Rgba c0(*pColor);
					cData += _colorPointer._stride;
					Rgba c1(*pColor);
					cData += _colorPointer._stride;

					drawLine2D(p0, p1, c0, c1);
				}

				if (_colorPointer._data == 0)
				{
					cData = (char*)colorPointerdesc._data;
				}

			}
			break;
		}
		case 4: //DM_TRIANGLES
		{
			for (int i = start; i < start + count; i += 3)
			{
				fData = (float*)posData;
				float3 p01(fData[0], fData[1], fData[2]);

				posData += _positionPointer._stride;
				fData = (float*)(posData);
				float3 p11(fData[0], fData[1], fData[2]);

				posData += _positionPointer._stride;
				fData = (float*)(posData);
				float3 p21(fData[0], fData[1], fData[2]);

				posData += _positionPointer._stride;

				p01 = p01 * _matrixModel;
				p11 = p11 * _matrixModel;
				p21 = p21 * _matrixModel;

				//转化为屏幕坐标
				if (!_ignoreOutside || _frust.pointInFrustum(p01)
					|| _frust.pointInFrustum(p11)
					|| _frust.pointInFrustum(p21)
					)
				{
					p01 = piplineTransform(p01);
					p11 = piplineTransform(p11);
					p21 = piplineTransform(p21);
					//if (DepthTest((uint)p01.x, (uint)p01.y, p01.z) && //需保证三角面间不互相穿过
					//	DepthTest((uint)p11.x, (uint)p11.y, p11.z) &&
					//	DepthTest((uint)p21.x, (uint)p21.y, p21.z))
					{
						int2 p0(p01.x, p01.y);
						int2 p1(p11.x, p11.y);
						int2 p2(p21.x, p21.y);

						Rgba * pColor = (Rgba*)cData;
						Rgba c0(*pColor);
						cData += _colorPointer._stride;
						Rgba c1(*pColor);
						cData += _colorPointer._stride;
						Rgba c2(*pColor);
						cData += _colorPointer._stride;

						float * pUV = (float*)uvData;
						float2 uv0(pUV[0], pUV[1]);
						uvData += _uvPointer._stride;
						pUV = (float*)uvData;
						float2 uv1(pUV[0], pUV[1]);
						uvData += _uvPointer._stride;
						pUV = (float*)uvData;
						float2 uv2(pUV[0], pUV[1]);
						uvData += _uvPointer._stride;

						Edge edges[3] =
						{
							Edge(p0.x,p0.y,c0, uv0,  p1.x,p1.y,c1, uv1),
							Edge(p1.x,p1.y,c1, uv1,  p2.x,p2.y,c2, uv2),
							Edge(p2.x,p2.y,c2, uv2,  p0.x,p0.y,c0, uv0),
						};
						drawTriangle(edges);
					}
					if (_colorPointer._data == 0)
					{
						cData = (char*)colorPointerdesc._data;
					}
					if (_uvPointer._data == 0)
					{
						uvData = (char*)uvPointerdesc._data;
					}
				}
			}
			break;
		}
			break;
		default:
			break;
		}
	}

	void Raster3D::drawImage(int startX, int startY, const Image * image, float alpha)
	{
		int left = tmax<int>(startX, 0);
		int top = tmax<int>(startY, 0);

		int right = tmin<int>(startX + image->width(), _width);
		int bottom = tmin<int>(startY + image->height(), _height);

		for (int x = left; x < right; ++x)
		{
			for (int y = top; y < bottom; ++y)
			{
				Rgba sourceColor = image->pixelAt(x - left, y - top);
				Rgba backgroundColor = getPixel(x, y);
				Rgba color = colorLerp(backgroundColor, sourceColor, sourceColor._a / 255.0f * alpha);
				setPixelEx(x, y, color);
			}
		}
	}

	void Raster3D::bindTexture(Image * image)
	{
		_texture = image;
	}

	void Raster3D::initializeZbuffer()
	{
		int count = _width * _height;
		for (int i = 0; i < count; ++i)
		{
			_Zbuffer[i] = 2.0f;
		}
	}

	void Raster3D::drawTriangle(Edge edges[3])
	{
		int iMax = 0;
		int length = edges[0]._y2 - edges[0]._y1;

		for (int i = 1; i < 3; ++i)
		{
			int len = edges[i]._y2 - edges[i]._y1;
			if (len > length)
			{
				length = len;
				iMax = i;
			}
		}
		int iShort1 = (iMax + 1) % 3;
		int iShort2 = (iMax + 2) % 3;

		drawEdge(edges[iMax], edges[iShort1], _texture);
		drawEdge(edges[iMax], edges[iShort2], _texture);

	}

	void Raster3D::drawLine2D(int2 pt1, int2 pt2, Rgba color1, Rgba color2)
	{
		float xOffset = pt1.x - pt2.x;
		float yOffset = pt1.y - pt2.y;

		if (xOffset == 0 && yOffset == 0)
		{
			setPixel(pt1.x, pt1.y, color1);
		}
		if (fabs(xOffset) > fabs(yOffset))
		{
			float xMin;
			float xMax;
			if (pt1.x < pt2.x)
			{
				xMin = pt1.x;
				xMax = pt2.x;
			}
			else
			{
				xMin = pt2.x;
				xMax = pt1.x;
			}

			float lenth = xMax - xMin;
			float slope = yOffset / xOffset;
			for (float x = xMin; x <= xMax; x += 255)
			{
				float   y = pt1.y + (x - pt1.x) * slope;
				float   scaler = (x - xMin) / lenth;
				Rgba    color = colorLerp(color1, color2, scaler);
				setPixel(x, y, color);
			}
		}
		else
		{
			float   yMin;
			float   yMax;
			if (pt1.y < pt2.y)
			{
				yMin = pt1.y;
				yMax = pt2.y;
			}
			else
			{
				yMin = pt2.y;
				yMax = pt1.y;
			}

			float   lenth = yMax - yMin;
			float   slope = xOffset / yOffset;
			for (float y = yMin; y <= yMax; y += 255)
			{
				float   x = pt1.x + (y - pt1.y) * slope;
				float   scaler = (y - yMin) / lenth;
				Rgba    color = colorLerp(color1, color2, scaler);
				setPixel(x, y, color);
			}
		}
	}


	//管线转化

	float3 Raster3D::piplineTransform(float3 pos)
	{
		float4 world(pos.x, pos.y, pos.z, 1);

		float4 screen = _matrixProjView * world;

		if (screen.w == 0.0f)
		{
			return false;
		}
		screen.x /= screen.w;
		screen.y /= screen.w;
		screen.z /= screen.w;

		// map to range 0 - 1
		screen.x = screen.x * 0.5f + 0.5f;
		screen.y = screen.y * 0.5f + 0.5f;
		screen.z = screen.z * 0.5f + 0.5f;

		// map to viewport
		screen.x = screen.x * _viewPort.x;
		screen.y = _height - screen.y * _viewPort.y; //使左下角成为0,0点

		return float3(screen.x, screen.y, screen.z);
	}

	void Raster3D::drawEdge(const Edge & e1, const Edge & e2, Image * image)
	{
		float yOffset1 = e1._y2 - e1._y1;
		if (yOffset1 == 0)
		{
			return;
		}

		float yOffset2 = e2._y2 - e2._y1;
		if (yOffset2 == 0)
		{
			return;
		}

		float xOffset = e2._x2 - e2._x1;
		float scale = 0;
		float step = 1.0f / yOffset2;


		float xOffset1 = e1._x2 - e1._x1;
		float scale1 = (float)(e2._y1 - e1._y1) / yOffset1;
		float step1 = 1.0f / yOffset1;

		for (int y = e2._y1; y < e2._y2; ++y)
		{
			int x1 = e1._x1 + (int)(scale1 * xOffset1);
			int x2 = e2._x1 + (int)(scale * xOffset);
			Rgba color2 = colorLerp(e2._color1, e2._color2, scale);
			Rgba color1 = colorLerp(e1._color1, e1._color2, scale1);

			float2 uvStart = uvLerp(e1._uv1, e1._uv2, scale1);
			float2 uvEnd = uvLerp(e2._uv1, e2._uv2, scale);

			Span span(x1, x2, y, color1, color2, uvStart, uvEnd);
			drawSpan(span, image);

			scale += step;
			scale1 += step1;
		}
	}

	void Raster3D::drawSpan(const Span & span, Image * image)
	{
		float length = span._xEnd - span._xStart;
		float scale = 0;
		float step = 1.0f / length;
		for (int x = span._xStart; x < span._xEnd; ++x)
		{
			float2 uv = uvLerp(span._uvStart, span._uvEnd, scale);
			
			//Rgba color = colorLerp(span._colorStart, span._colorEnd, scale);
		
			Rgba color = image->pixelUV(uv.x, uv.y);

			scale += step;

			setPixel(x, span._y, color);
		}
	}

	Span::Span(int xStart, int xEnd, int y, Rgba colorStart, Rgba colorEnd, float2 uvStart, float2 uvEnd)
	{
		if (xStart < xEnd)
		{
			_xStart = xStart;
			_xEnd = xEnd;
			_colorStart = colorStart;
			_colorEnd = colorEnd;
			_uvStart = uvStart;
			_uvEnd = uvEnd;
			_y = y;
		}
		else
		{
			_xStart = xEnd;
			_xEnd = xStart;
			_colorStart = colorEnd;
			_colorEnd = colorStart;
			_uvStart = uvEnd;
			_uvEnd = uvStart;
			_y = y;
		}
	}
	Edge::Edge(int x1, int y1, Rgba color1, float2 uv1, int x2, int y2, Rgba color2, float2 uv2)
	{
		if (y1 < y2)
		{
			_x1 = x1;
			_y1 = y1;
			_uv1 = uv1;
			_color1 = color1;

			_x2 = x2;
			_y2 = y2;
			_uv2 = uv2;
			_color2 = color2;
		}
		else
		{
			_x1 = x2;
			_y1 = y2;
			_uv1 = uv2;
			_color1 = color2;

			_x2 = x1;
			_y2 = y1;
			_uv2 = uv1;
			_color2 = color1;
		}
	}
	
}

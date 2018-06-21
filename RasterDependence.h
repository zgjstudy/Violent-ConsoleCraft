#pragma once

#include "MATH.h"

namespace RENDERER
{
	enum DRAWMODE
	{
		DM_POINTS = 0,	//画点
		DM_LINES = 1,	//画线段
		DM_LINE_LOOP = 2,	//所有点依次首尾相连
		DM_LINE_STRIP = 3,	//所有点依次首尾不相连
		DM_TRIANGLES = 4,	//画三角形
	};
	enum DATATYPE
	{
		DT_BYTE,			//
		DT_FLOAT,			//
		DT_DOUBLE,			//
	};
	struct DataElementDesc			//数据元素描述
	{
		int				_size;		//数据个数
		DATATYPE		_type;		//数据种类
		int				_stride;	//单位偏移量
		const void *	_data;		//数据地址
	};
	struct Vertex
	{
		float x, y, z;
		float u, v;
		Rgba color;
	};
}
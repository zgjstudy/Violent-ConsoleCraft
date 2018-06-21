#pragma once

#include "MATH.h"

namespace RENDERER
{
	enum DRAWMODE
	{
		DM_POINTS = 0,	//����
		DM_LINES = 1,	//���߶�
		DM_LINE_LOOP = 2,	//���е�������β����
		DM_LINE_STRIP = 3,	//���е�������β������
		DM_TRIANGLES = 4,	//��������
	};
	enum DATATYPE
	{
		DT_BYTE,			//
		DT_FLOAT,			//
		DT_DOUBLE,			//
	};
	struct DataElementDesc			//����Ԫ������
	{
		int				_size;		//���ݸ���
		DATATYPE		_type;		//��������
		int				_stride;	//��λƫ����
		const void *	_data;		//���ݵ�ַ
	};
	struct Vertex
	{
		float x, y, z;
		float u, v;
		Rgba color;
	};
}
#include "Block.h"

using RENDERER::float3;
using RENDERER::int3;
using RENDERER::Image;

CONSOLEWORLD::Block::Block()
{
}

CONSOLEWORLD::Block::~Block()
{
}

void CONSOLEWORLD::Block::Break()
{
	_ID = 0;
	_topImage = nullptr;
	_sideImage = nullptr;
	_bottomImage = nullptr;
}

void CONSOLEWORLD::Block::initialize(int3 coord)
{
	_coord = coord;
	_pointCoords[0] = float3(BLOCK_EDGE_LENGTH * coord.x, BLOCK_EDGE_LENGTH * coord.y, BLOCK_EDGE_LENGTH * coord.z);
	_pointCoords[1] = _pointCoords[0] + float3(0, BLOCK_EDGE_LENGTH, 0);
	_pointCoords[2] = _pointCoords[1] + float3(BLOCK_EDGE_LENGTH, 0, 0);
	_pointCoords[3] = _pointCoords[2] - float3(0, BLOCK_EDGE_LENGTH, 0);
	_pointCoords[4] = _pointCoords[2] + float3(0, 0, BLOCK_EDGE_LENGTH);
	_pointCoords[5] = _pointCoords[4] - float3(BLOCK_EDGE_LENGTH, 0, 0);
	_pointCoords[6] = _pointCoords[5] - float3(0, BLOCK_EDGE_LENGTH, 0);
	_pointCoords[7] = _pointCoords[6] + float3(BLOCK_EDGE_LENGTH, 0, 0);
	_midCoord = _pointCoords[0] + float3(BLOCK_EDGE_LENGTH/2);

	aabb.setExtents(_pointCoords[0], _pointCoords[4]);
	_CollisionBox.setExtents(_pointCoords[0] - float3(20), _pointCoords[4] + float3(20));
}

int3 CONSOLEWORLD::Block::getCoord()
{
	return _coord;
}

float3 CONSOLEWORLD::Block::getPointCoord(int i)
{
	if (i < 0 || i > 7)
	{
		return float3(-1);
	}
	return _pointCoords[i];
}

int CONSOLEWORLD::Block::getID() const
{
	return _ID;
}

void CONSOLEWORLD::Block::setID(int id)
{
	_ID = id;
}

float3 CONSOLEWORLD::Block::getMidCoord()
{
	return _midCoord;
}

void CONSOLEWORLD::Block::setTopImage(Image * image)
{
	_topImage = image;
}

Image * CONSOLEWORLD::Block::getTopImage()
{
	return _topImage;
}

void CONSOLEWORLD::Block::setSideImage(Image * image)
{
	_sideImage = image;
}

Image * CONSOLEWORLD::Block::getSideImage()
{
	return _sideImage;
}

void CONSOLEWORLD::Block::setBottomImage(Image * image)
{
	_bottomImage = image;
}

Image * CONSOLEWORLD::Block::getBottomImage()
{
	return _bottomImage;
}

void CONSOLEWORLD::Block::RenderTop(RENDERER::Raster3D & raster)
{
	RENDERER::Vertex vertexsTop[] =
	{
		{ _pointCoords[1].x, _pointCoords[1].y , _pointCoords[1].z , 0,	0,	Rgba() },
	{ _pointCoords[2].x, _pointCoords[2].y , _pointCoords[2].z , 1,	0,	Rgba() },
	{ _pointCoords[4].x, _pointCoords[4].y , _pointCoords[4].z , 1,	1,	Rgba() },

	{ _pointCoords[1].x, _pointCoords[1].y , _pointCoords[1].z , 0,	0,	Rgba() },
	{ _pointCoords[4].x, _pointCoords[4].y , _pointCoords[4].z , 1,	1,	Rgba() },
	{ _pointCoords[5].x, _pointCoords[5].y , _pointCoords[5].z , 0,	1,	Rgba() },
	};
	raster.bindTexture(_topImage);

	raster.setVertexPointer(2, RENDERER::DT_FLOAT, sizeof(RENDERER::Vertex), &vertexsTop[0].x);
	raster.setTextureCoordPointer(2, RENDERER::DT_FLOAT, sizeof(RENDERER::Vertex), &vertexsTop[0].u);
	raster.drawArrays(RENDERER::DM_TRIANGLES, 0, 6);
}

void CONSOLEWORLD::Block::RenderBottom(RENDERER::Raster3D & raster)
{
	RENDERER::Vertex vertexsBottom[] =
	{
		{ _pointCoords[0].x, _pointCoords[0].y , _pointCoords[0].z , 0,	0,	Rgba() },
	{ _pointCoords[3].x, _pointCoords[3].y , _pointCoords[3].z , 1,	0,	Rgba() },
	{ _pointCoords[3].x, _pointCoords[7].y , _pointCoords[7].z , 1,	1,	Rgba() },

	{ _pointCoords[0].x, _pointCoords[0].y , _pointCoords[0].z , 0,	0,	Rgba() },
	{ _pointCoords[7].x, _pointCoords[7].y , _pointCoords[7].z , 1,	1,	Rgba() },
	{ _pointCoords[6].x, _pointCoords[6].y , _pointCoords[6].z , 0,	1,	Rgba() },
	};
	raster.bindTexture(_bottomImage);

	raster.setVertexPointer(2, RENDERER::DT_FLOAT, sizeof(RENDERER::Vertex), &vertexsBottom[0].x);
	raster.setTextureCoordPointer(2, RENDERER::DT_FLOAT, sizeof(RENDERER::Vertex), &vertexsBottom[0].u);
	raster.drawArrays(RENDERER::DM_TRIANGLES, 0, 6);
}

void CONSOLEWORLD::Block::RenderRight(RENDERER::Raster3D & raster)
{
	RENDERER::Vertex vertex[] =
	{
		{ _pointCoords[4].x, _pointCoords[4].y , _pointCoords[4].z , 0,	0,	Rgba() },
	{ _pointCoords[2].x, _pointCoords[2].y , _pointCoords[2].z , 1,	0,	Rgba() },
	{ _pointCoords[3].x, _pointCoords[3].y , _pointCoords[3].z , 1,	1,	Rgba() },

	{ _pointCoords[4].x, _pointCoords[4].y , _pointCoords[4].z , 0,	0,	Rgba() },
	{ _pointCoords[3].x, _pointCoords[3].y , _pointCoords[3].z , 1,	1,	Rgba() },
	{ _pointCoords[7].x, _pointCoords[7].y , _pointCoords[7].z , 0,	1,	Rgba() },
	};
	raster.bindTexture(_sideImage);

	raster.setVertexPointer(2, RENDERER::DT_FLOAT, sizeof(RENDERER::Vertex), &vertex[0].x);
	raster.setTextureCoordPointer(2, RENDERER::DT_FLOAT, sizeof(RENDERER::Vertex), &vertex[0].u);
	raster.drawArrays(RENDERER::DM_TRIANGLES, 0, 6);
}

void CONSOLEWORLD::Block::RenderLeft(RENDERER::Raster3D & raster)
{
	RENDERER::Vertex vertex[] =
	{
		{ _pointCoords[1].x, _pointCoords[1].y , _pointCoords[1].z , 0,	0,	Rgba() },
	{ _pointCoords[5].x, _pointCoords[5].y , _pointCoords[5].z , 1,	0,	Rgba() },
	{ _pointCoords[6].x, _pointCoords[6].y , _pointCoords[6].z , 1,	1,	Rgba() },

	{ _pointCoords[1].x, _pointCoords[1].y , _pointCoords[1].z , 0,	0,	Rgba() },
	{ _pointCoords[6].x, _pointCoords[6].y , _pointCoords[6].z , 1,	1,	Rgba() },
	{ _pointCoords[0].x, _pointCoords[0].y , _pointCoords[0].z , 0,	1,	Rgba() },
	};
	raster.bindTexture(_sideImage);

	raster.setVertexPointer(2, RENDERER::DT_FLOAT, sizeof(RENDERER::Vertex), &vertex[0].x);
	raster.setTextureCoordPointer(2, RENDERER::DT_FLOAT, sizeof(RENDERER::Vertex), &vertex[0].u);
	raster.drawArrays(RENDERER::DM_TRIANGLES, 0, 6);
}

void CONSOLEWORLD::Block::RenderFront(RENDERER::Raster3D & raster)
{
	RENDERER::Vertex vertex[] =
	{
		{ _pointCoords[5].x, _pointCoords[5].y , _pointCoords[5].z , 0,	0,	Rgba() },
	{ _pointCoords[4].x, _pointCoords[4].y , _pointCoords[4].z , 1,	0,	Rgba() },
	{ _pointCoords[7].x, _pointCoords[7].y , _pointCoords[7].z , 1,	1,	Rgba() },

	{ _pointCoords[5].x, _pointCoords[5].y , _pointCoords[5].z , 0,	0,	Rgba() },
	{ _pointCoords[7].x, _pointCoords[7].y , _pointCoords[7].z , 1,	1,	Rgba() },
	{ _pointCoords[6].x, _pointCoords[6].y , _pointCoords[6].z , 0,	1,	Rgba() },
	};
	raster.bindTexture(_sideImage);

	raster.setVertexPointer(2, RENDERER::DT_FLOAT, sizeof(RENDERER::Vertex), &vertex[0].x);
	raster.setTextureCoordPointer(2, RENDERER::DT_FLOAT, sizeof(RENDERER::Vertex), &vertex[0].u);
	raster.drawArrays(RENDERER::DM_TRIANGLES, 0, 6);
}

void CONSOLEWORLD::Block::RenderBack(RENDERER::Raster3D & raster)
{
	RENDERER::Vertex vertex[] =
	{
		{ _pointCoords[2].x, _pointCoords[2].y , _pointCoords[2].z , 0,	0,	Rgba() },
	{ _pointCoords[1].x, _pointCoords[1].y , _pointCoords[1].z , 1,	0,	Rgba() },
	{ _pointCoords[0].x, _pointCoords[0].y , _pointCoords[0].z , 1,	1,	Rgba() },

	{ _pointCoords[2].x, _pointCoords[2].y , _pointCoords[2].z , 0,	0,	Rgba() },
	{ _pointCoords[0].x, _pointCoords[0].y , _pointCoords[0].z , 1,	1,	Rgba() },
	{ _pointCoords[3].x, _pointCoords[3].y , _pointCoords[3].z , 0,	1,	Rgba() },
	};
	raster.bindTexture(_sideImage);

	raster.setVertexPointer(2, RENDERER::DT_FLOAT, sizeof(RENDERER::Vertex), &vertex[0].x);
	raster.setTextureCoordPointer(2, RENDERER::DT_FLOAT, sizeof(RENDERER::Vertex), &vertex[0].u);
	raster.drawArrays(RENDERER::DM_TRIANGLES, 0, 6);
}

void CONSOLEWORLD::Block::RenderBlock(RENDERER::Raster3D & raster, float3 eye, int3 Bcoord)
{
	if (_ID == 0)
		return;
	if (getCoord() == Bcoord)
	{
		
	}

	RELATIVE_POSITIVE rp;

	if (eye.x >= _midCoord.x && eye.y >= _midCoord.y && eye.z >= _midCoord.z)
		rp = RP_TOP_RIGHT_FRONT;
	else if (eye.x >= _midCoord.x && eye.y >= _midCoord.y && eye.z < _midCoord.z)
		rp = RP_TOP_RIGHT_BACK;
	else if (eye.x < _midCoord.x && eye.y >= _midCoord.y && eye.z >= _midCoord.z)
		rp = RP_TOP_LEFT_FRONT;
	else if (eye.x < _midCoord.x && eye.y >= _midCoord.y && eye.z < _midCoord.z)
		rp = RP_TOP_LEFT_BACK;
	else if (eye.x < _midCoord.x && eye.y < _midCoord.y && eye.z < _midCoord.z)
		rp = RP_BOTTOM_LEFT_BACK;
	else if (eye.x < _midCoord.x && eye.y < _midCoord.y && eye.z >= _midCoord.z)
		rp = RP_BOTTOM_LEFT_FRONT;
	else if (eye.x >= _midCoord.x && eye.y < _midCoord.y && eye.z < _midCoord.z)
		rp = RP_BOTTOM_RIGHT_BACK;
	else if (eye.x >= _midCoord.x && eye.y < _midCoord.y && eye.z >= _midCoord.z)
		rp = RP_BOTTOM_RIGHT_FRONT;

	switch (rp)
	{
	case Block::RP_RIGHT:
		break;
	case Block::RP_LEFT:
		break;
	case Block::RP_TOP:
		break;
	case Block::RP_BOTTOM:
		break;
	case Block::RP_FRONT:
		break;
	case Block::RP_BACK:
		break;
	case Block::RP_TOP_RIGHT_FRONT:
		RenderRight(raster);
		RenderFront(raster);
		RenderTop(raster);
		break;
	case Block::RP_BOTTOM_RIGHT_FRONT:
		RenderRight(raster);
		RenderFront(raster);
		RenderBottom(raster);
		break;
	case Block::RP_TOP_RIGHT_BACK:
		RenderRight(raster);
		RenderBack(raster);
		RenderTop(raster);
		break;
	case Block::RP_BOTTOM_RIGHT_BACK:
		RenderRight(raster);
		RenderBack(raster);
		RenderBottom(raster);
		break;
	case Block::RP_TOP_LEFT_FRONT:
		RenderLeft(raster);
		RenderFront(raster);
		RenderTop(raster);
		break;
	case Block::RP_BOTTOM_LEFT_FRONT:
		RenderLeft(raster);
		RenderFront(raster);
		RenderBottom(raster);
		break;
	case Block::RP_TOP_LEFT_BACK:
		RenderLeft(raster);
		RenderBack(raster);
		RenderTop(raster);
		break;
	case Block::RP_BOTTOM_LEFT_BACK:
		RenderLeft(raster);
		RenderBack(raster);
		RenderBottom(raster);
		break;
	default:
		break;
	}
}

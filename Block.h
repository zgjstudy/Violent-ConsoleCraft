#pragma once

#include "Raster3D.h"
#include "Image.h"
#include "MATH.h"

namespace CONSOLEWORLD
{
	using RENDERER::Rgba;
	using RENDERER::float3;
	using RENDERER::int3;
	using RENDERER::Image;

#define BLOCK_EDGE_LENGTH 100


	//-------------------- 方块类
	class Block
	{
	private:
		int3		_coord			= int3(0,0,0);
		int			_ID				= 0;
		Image *		_topImage		= nullptr;
		Image *		_sideImage		= nullptr;
		Image *		_bottomImage	= nullptr;
		float3		_midCoord;
		float3		_pointCoords[8]	= {};

	public:
		ZMATH::AABB3D _CollisionBox;
		ZMATH::AABB3D aabb;

		enum RELATIVE_POSITIVE
		{
			RP_RIGHT					=	0,
			RP_LEFT,
			RP_TOP,
			RP_BOTTOM,
			RP_FRONT,
			RP_BACK,
			RP_TOP_RIGHT_FRONT,
			RP_BOTTOM_RIGHT_FRONT,
			RP_TOP_RIGHT_BACK,
			RP_BOTTOM_RIGHT_BACK,
			RP_TOP_LEFT_FRONT,
			RP_BOTTOM_LEFT_FRONT,
			RP_TOP_LEFT_BACK,
			RP_BOTTOM_LEFT_BACK,
		};

		/*
		   1-----2
		  /|    /|
	 	 / |   / |
		5-----4  |
		|  0--|--3
		| /   | /
		|/    |/
		6-----7
		*/

		/*
				    /\ y
				    |
					|
					|
					|
					| (0,0,0)
					|________________\ x
					/				 /
				   /
				  /
				 /
				/
			   V z
		*/

		
		Block();

		~Block();

		//-------------------- 功能函数
		void Break();

		void initialize(int3 coord);

		int3 getCoord();
		
		float3 getPointCoord(int i);

		int getID() const;
		
		void setID(int id);

		float3 getMidCoord();

		void setTopImage(Image * image);

		Image * getTopImage();

		void setSideImage(Image * image);

		Image * getSideImage();

		void setBottomImage(Image * image);

		Image * getBottomImage();
		
		void RenderBlock(RENDERER::Raster3D & raster, float3 eye, int3 Bcoord);
		
	private:
		void RenderTop(RENDERER::Raster3D & raster);

		void RenderBottom(RENDERER::Raster3D & raster);

		void RenderRight(RENDERER::Raster3D & raster);

		void RenderLeft(RENDERER::Raster3D & raster);
		
		void RenderFront(RENDERER::Raster3D & raster);
		
		void RenderBack(RENDERER::Raster3D & raster);	
	};
}
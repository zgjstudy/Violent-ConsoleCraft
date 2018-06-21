#pragma once

#include <cmath>
#include <fstream>
#include "Block.h"
#include "Raster3D.h"
#include "MaterialIndex.h"

namespace CONSOLEWORLD
{
	using namespace ZMATH;
	//-------------------- 世界信息
#define WORLD_LENGTH 20
#define WORLD_WIDTH 20
#define WORLD_HEIGHT 15

	//-------------------- 世界类
	class World
	{
	public:
		Block _blocks[WORLD_WIDTH][WORLD_HEIGHT][WORLD_LENGTH];

		World();


		Block & getBloctAt(int3 coord);

		Block & getBloctAt(int x, int y, int z);

		bool setBlockAt(int3 coord, Block & block);

		//-------------------- 移动时的碰撞检测(bug)
		bool MoveCollideCheckX(AABB3D & aabb);

		bool MoveCollideCheckY(AABB3D & aabb);

		bool MoveCollideCheckZ(AABB3D & aabb);

		//-------------------- 获取指向的方块坐标
		int3 GetPointedBlock(Ray ray);

		//-------------------- 获取待定的方块放置位置
		int3 GetSelectedPosition(Ray ray, int3 Bcoord);

		//-------------------- 使某一平面铺满某方块
		void setY(int y, Block & block);

		void setX(int x, Block & block);

		void setZ(int z, Block & block);

		//-------------------- 渲染世界
		void RenderWorld(RENDERER::Raster3D & raster, float3 eye, int3 Bcoord);

		//-------------------- 读档辅助
		void set(Block*AllBlock);

		//-------------------- 存档函数
		bool Save(char *szPath, int num);

		//-------------------- 读档函数
		bool Load(char *szPath, int num);

		bool ResetSave(char *szPath, int num);
	};
}
#pragma once

#include <cmath>
#include <fstream>
#include "Block.h"
#include "Raster3D.h"
#include "MaterialIndex.h"

namespace CONSOLEWORLD
{
	using namespace ZMATH;
	//-------------------- ������Ϣ
#define WORLD_LENGTH 20
#define WORLD_WIDTH 20
#define WORLD_HEIGHT 15

	//-------------------- ������
	class World
	{
	public:
		Block _blocks[WORLD_WIDTH][WORLD_HEIGHT][WORLD_LENGTH];

		World();


		Block & getBloctAt(int3 coord);

		Block & getBloctAt(int x, int y, int z);

		bool setBlockAt(int3 coord, Block & block);

		//-------------------- �ƶ�ʱ����ײ���(bug)
		bool MoveCollideCheckX(AABB3D & aabb);

		bool MoveCollideCheckY(AABB3D & aabb);

		bool MoveCollideCheckZ(AABB3D & aabb);

		//-------------------- ��ȡָ��ķ�������
		int3 GetPointedBlock(Ray ray);

		//-------------------- ��ȡ�����ķ������λ��
		int3 GetSelectedPosition(Ray ray, int3 Bcoord);

		//-------------------- ʹĳһƽ������ĳ����
		void setY(int y, Block & block);

		void setX(int x, Block & block);

		void setZ(int z, Block & block);

		//-------------------- ��Ⱦ����
		void RenderWorld(RENDERER::Raster3D & raster, float3 eye, int3 Bcoord);

		//-------------------- ��������
		void set(Block*AllBlock);

		//-------------------- �浵����
		bool Save(char *szPath, int num);

		//-------------------- ��������
		bool Load(char *szPath, int num);

		bool ResetSave(char *szPath, int num);
	};
}
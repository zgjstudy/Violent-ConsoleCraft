#include "ConsoleCraftWorld.h"
using namespace ZMATH;

CONSOLEWORLD::World::World()
{
	for (int x = 0; x < WORLD_WIDTH; ++x)
		for (int y = 0; y < WORLD_HEIGHT; ++y)
			for (int z = 0; z < WORLD_LENGTH; ++z)
			{
				_blocks[x][y][z].initialize(int3(x, y, z));
			}
}

CONSOLEWORLD::Block & CONSOLEWORLD::World::getBloctAt(int3 coord)
{
	if (coord.x < 0 || coord.x < 0 || coord.x < 0 || coord.x >= WORLD_WIDTH || coord.y >= WORLD_HEIGHT || coord.z >= WORLD_LENGTH)
		return _blocks[0][0][0];
	return _blocks[coord.x][coord.y][coord.z];
}

CONSOLEWORLD::Block & CONSOLEWORLD::World::getBloctAt(int x, int y, int z)
{
	return getBloctAt(int3(x, y, z));
}

bool CONSOLEWORLD::World::setBlockAt(int3 coord, Block & block)
{
	if (coord.x < 0 || coord.x < 0 || coord.x < 0 || coord.x >= WORLD_WIDTH || coord.y >= WORLD_HEIGHT || coord.z >= WORLD_LENGTH)
		return false;

	_blocks[coord.x][coord.y][coord.z].setID(block.getID());
	_blocks[coord.x][coord.y][coord.z].setTopImage(block.getTopImage());
	_blocks[coord.x][coord.y][coord.z].setSideImage(block.getSideImage());
	_blocks[coord.x][coord.y][coord.z].setBottomImage(block.getBottomImage());
	return true;
}


//-------------------- 移动时的碰撞检测(bug)

bool CONSOLEWORLD::World::MoveCollideCheckX(AABB3D & aabb)
{
	float3 pos = aabb.getCenter();

	int3 coord(floor(pos.x / BLOCK_EDGE_LENGTH), floor(pos.y / BLOCK_EDGE_LENGTH), floor(pos.z / BLOCK_EDGE_LENGTH));
	if (coord.x < -1 || coord.x > WORLD_WIDTH)
		return true;
	else if (coord.x >= 0 && coord.x < WORLD_WIDTH)
	{
		//如果aabb不相交
		if (getBloctAt(coord - int3(1, 0, 0)).getID() ? !(aabb.intersects(getBloctAt(coord - int3(1, 0, 0))._CollisionBox)) : 1)
			if (getBloctAt(coord + int3(1, 0, 0)).getID() ? !(aabb.intersects(getBloctAt(coord + int3(1, 0, 0))._CollisionBox)) : 1)
				return true;
	}
	else if (coord.x == -1)
	{
		if (getBloctAt(coord + int3(1, 0, 0)).getID() ? !(aabb.intersects(getBloctAt(coord + int3(1, 0, 0))._CollisionBox)) : 1)
			return true;
	}
	else if (coord.x == WORLD_WIDTH)
	{
		if (getBloctAt(coord - int3(1, 0, 0)).getID() ? !(aabb.intersects(getBloctAt(coord - int3(1, 0, 0))._CollisionBox)) : 1)
			return true;
	}

	return false;
}

bool CONSOLEWORLD::World::MoveCollideCheckY(AABB3D & aabb)
{
	float3 pos = aabb.getCenter();

	int3 coord(floor(pos.x / BLOCK_EDGE_LENGTH), floor(pos.y / BLOCK_EDGE_LENGTH), floor(pos.z / BLOCK_EDGE_LENGTH));

	if (coord.y < -1 || coord.y > WORLD_HEIGHT)
		return true;
	else if (coord.y >= 0 && coord.y < WORLD_HEIGHT)
	{
		//如果aabb不相交
		if (getBloctAt(coord - int3(0, 1, 0)).getID() ? !(aabb.intersects(getBloctAt(coord - int3(0, 1, 0))._CollisionBox)) : 1)
			if (getBloctAt(coord + int3(0, 1, 0)).getID() ? !(aabb.intersects(getBloctAt(coord + int3(0, 1, 0))._CollisionBox)) : 1)
				return true;
	}
	else if (coord.y == -1)
	{
		if (getBloctAt(coord + int3(0, 1, 0)).getID() ? !(aabb.intersects(getBloctAt(coord + int3(0, 1, 0))._CollisionBox)) : 1)
			return true;
	}
	else if (coord.y == WORLD_HEIGHT)
	{
		if (getBloctAt(coord - int3(0, 1, 0)).getID() ? !(aabb.intersects(getBloctAt(coord - int3(0, 1, 0))._CollisionBox)) : 1)
			return true;
	}

	return false;

}

bool CONSOLEWORLD::World::MoveCollideCheckZ(AABB3D & aabb)
{
	float3 pos = aabb.getCenter();

	int3 coord(floor(pos.x / BLOCK_EDGE_LENGTH), floor(pos.y / BLOCK_EDGE_LENGTH), floor(pos.z / BLOCK_EDGE_LENGTH));

	if (coord.z < -1 || coord.z > WORLD_LENGTH)
		return true;
	else if (coord.z >= 0 && coord.z < WORLD_LENGTH)
	{
		//如果aabb不相交
		if (getBloctAt(coord - int3(0, 0, 1)).getID() ? !(aabb.intersects(getBloctAt(coord - int3(0, 0, 1))._CollisionBox)) : 1)
			if (getBloctAt(coord + int3(0, 0, 1)).getID() ? !(aabb.intersects(getBloctAt(coord + int3(0, 0, 1))._CollisionBox)) : 1)
				return true;
	}
	else if (coord.z == -1)
	{
		if (getBloctAt(coord + int3(0, 0, 1)).getID() ? !(aabb.intersects(getBloctAt(coord + int3(0, 0, 1))._CollisionBox)) : 1)
			return true;
	}
	else if (coord.z == WORLD_LENGTH)
	{
		if (getBloctAt(coord - int3(0, 0, 1)).getID() ? !(aabb.intersects(getBloctAt(coord - int3(0, 0, 1))._CollisionBox)) : 1)
			return true;
	}

	return false;
}

//-------------------- 获取指向的方块坐标

int3 CONSOLEWORLD::World::GetPointedBlock(Ray ray)
{
	float3 origin = ray.getOrigin();
	int minD = 1e8;
	int dis;
	int3 res(-1);
	for (int x = 0; x < WORLD_WIDTH; ++x)
		for (int y = 0; y < WORLD_HEIGHT; ++y)
			for (int z = 0; z < WORLD_LENGTH; ++z)
			{
				if (int3(x, y, z) == int3(floor(origin.x / BLOCK_EDGE_LENGTH), floor(origin.y / BLOCK_EDGE_LENGTH), floor(origin.z / BLOCK_EDGE_LENGTH)))
					continue;
				if (getBloctAt(x, y, z).getID() != 0)
					if ((ray.intersects(getBloctAt(int3(x, y, z)).aabb)).first)
					{
						dis = distance(float3(x, y, z), origin);
						if (dis < minD)
						{
							res = int3(x, y, z);
						}
					}
			}
	return res;
}

//-------------------- 获取待定的方块放置位置

ZMATH::int3 CONSOLEWORLD::World::GetSelectedPosition(Ray ray, int3 Bcoord)
{
	const float eps = 1e-5;
	float3 eye = ray.getOrigin();
	float3 blockMid = getBloctAt(Bcoord).getMidCoord();
	float dis = ray.intersects(getBloctAt(Bcoord).aabb).second;
	float3 intersection = ray.getPoint(dis);
	float3 corner(
		eye.x > blockMid.x ? blockMid.x + BLOCK_EDGE_LENGTH / 2 : blockMid.x - BLOCK_EDGE_LENGTH / 2,
		eye.y > blockMid.y ? blockMid.y + BLOCK_EDGE_LENGTH / 2 : blockMid.y - BLOCK_EDGE_LENGTH / 2,
		eye.z > blockMid.z ? blockMid.z + BLOCK_EDGE_LENGTH / 2 : blockMid.z - BLOCK_EDGE_LENGTH / 2
	);

	if (fabs(intersection.x - corner.x) < eps)
	{
		if (corner.x > blockMid.x)
			return Bcoord + int3(1, 0, 0);
		else
			return Bcoord - int3(1, 0, 0);
	}
	if (fabs(intersection.y - corner.y) < eps)
	{
		if (corner.y > blockMid.y)
			return Bcoord + int3(0, 1, 0);
		else
			return Bcoord - int3(0, 1, 0);
	}
	if (fabs(intersection.z - corner.z) < eps)
	{
		if (corner.z > blockMid.z)
			return Bcoord + int3(0, 0, 1);
		else
			return Bcoord - int3(0, 0, 1);
	}
	return int3(-1);
}

//-------------------- 使某一平面铺满某方块

void CONSOLEWORLD::World::setY(int y, Block & block)
{
	for (int x = 0; x < WORLD_WIDTH; ++x)
		for (int z = 0; z < WORLD_LENGTH; ++z)
			setBlockAt(int3(x, y, z), block);
}

void CONSOLEWORLD::World::setX(int x, Block & block)
{
	for (int y = 0; y < WORLD_WIDTH; ++y)
		for (int z = 0; z < WORLD_LENGTH; ++z)
			setBlockAt(int3(x, y, z), block);
}

void CONSOLEWORLD::World::setZ(int z, Block & block)
{
	for (int y = 0; y < WORLD_WIDTH; ++y)
		for (int x = 0; x < WORLD_LENGTH; ++x)
			setBlockAt(int3(x, y, z), block);
}

//-------------------- 渲染世界

void CONSOLEWORLD::World::RenderWorld(RENDERER::Raster3D & raster, float3 eye, int3 Bcoord)
{
	int3 ec(floor(eye.x / BLOCK_EDGE_LENGTH), floor(eye.y / BLOCK_EDGE_LENGTH), floor(eye.z / BLOCK_EDGE_LENGTH));

	if (ec.x < 0)
	{
		ec.x = -1;
	}
	if (ec.y < 0)
	{
		ec.y = -1;
	}
	if (ec.z < 0)
	{
		ec.z = -1;
	}
	if (ec.x > WORLD_WIDTH - 1)
	{
		ec.x = WORLD_WIDTH;
	}
	if (ec.y > WORLD_HEIGHT - 1)
	{
		ec.y = WORLD_HEIGHT;
	}
	if (ec.z > WORLD_LENGTH - 1)
	{
		ec.z = WORLD_LENGTH;
	}

	for (int x = 0; x < ec.x; ++x)
		for (int y = 0; y < ec.y; ++y)
			for (int z = 0; z < ec.z; ++z)
			{
				getBloctAt(x, y, z).RenderBlock(raster, eye, Bcoord);
			}
	for (int x = WORLD_WIDTH - 1; x > ec.x; --x)
		for (int y = 0; y < ec.y; ++y)
			for (int z = 0; z < ec.z; ++z)
			{
				getBloctAt(x, y, z).RenderBlock(raster, eye, Bcoord);
			}
	for (int x = 0; x < ec.x; ++x)
		for (int y = WORLD_HEIGHT - 1; y > ec.y; --y)
			for (int z = 0; z < ec.z; ++z)
			{
				getBloctAt(x, y, z).RenderBlock(raster, eye, Bcoord);
			}
	for (int x = WORLD_WIDTH - 1; x > ec.x; --x)
		for (int y = WORLD_HEIGHT - 1; y > ec.y; --y)
			for (int z = 0; z < ec.z; ++z)
			{
				getBloctAt(x, y, z).RenderBlock(raster, eye, Bcoord);
			}
	for (int x = 0; x < ec.x; ++x)
		for (int y = 0; y < ec.y; ++y)
			for (int z = WORLD_LENGTH - 1; z > ec.z; --z)
			{
				getBloctAt(x, y, z).RenderBlock(raster, eye, Bcoord);
			}
	for (int x = WORLD_WIDTH - 1; x > ec.x; --x)
		for (int y = 0; y < ec.y; ++y)
			for (int z = WORLD_LENGTH - 1; z > ec.z; --z)
			{
				getBloctAt(x, y, z).RenderBlock(raster, eye, Bcoord);
			}
	for (int x = 0; x < ec.x; ++x)
		for (int y = WORLD_HEIGHT - 1; y > ec.y; --y)
			for (int z = WORLD_LENGTH - 1; z > ec.z; --z)
			{
				getBloctAt(x, y, z).RenderBlock(raster, eye, Bcoord);
			}
	for (int x = WORLD_WIDTH - 1; x > ec.x; --x)
		for (int y = WORLD_HEIGHT - 1; y > ec.y; --y)
			for (int z = WORLD_LENGTH - 1; z > ec.z; --z)
			{
				getBloctAt(x, y, z).RenderBlock(raster, eye, Bcoord);
			}

	if (ec.x >= 0 && ec.x <= WORLD_WIDTH - 1)
	{
		for (int y = 0; y < ec.y; ++y)
			for (int z = 0; z < ec.z; ++z)
				getBloctAt(ec.x, y, z).RenderBlock(raster, eye, Bcoord);
		for (int y = WORLD_HEIGHT - 1; y > ec.y; --y)
			for (int z = 0; z < ec.z; ++z)
				getBloctAt(ec.x, y, z).RenderBlock(raster, eye, Bcoord);
		for (int y = 0; y < ec.y; ++y)
			for (int z = WORLD_LENGTH - 1; z > ec.z; --z)
				getBloctAt(ec.x, y, z).RenderBlock(raster, eye, Bcoord);
		for (int y = WORLD_HEIGHT - 1; y > ec.y; --y)
			for (int z = WORLD_LENGTH - 1; z > ec.z; --z)
				getBloctAt(ec.x, y, z).RenderBlock(raster, eye, Bcoord);
	}

	if (ec.y >= 0 && ec.y <= WORLD_HEIGHT - 1)
	{
		for (int x = 0; x < ec.x; ++x)
			for (int z = 0; z < ec.z; ++z)
				getBloctAt(x, ec.y, z).RenderBlock(raster, eye, Bcoord);
		for (int x = WORLD_WIDTH - 1; x > ec.x; --x)
			for (int z = 0; z < ec.z; ++z)
				getBloctAt(x, ec.y, z).RenderBlock(raster, eye, Bcoord);
		for (int x = 0; x < ec.x; ++x)
			for (int z = WORLD_LENGTH - 1; z > ec.z; --z)
				getBloctAt(x, ec.y, z).RenderBlock(raster, eye, Bcoord);
		for (int x = WORLD_WIDTH - 1; x > ec.x; --x)
			for (int z = WORLD_LENGTH - 1; z > ec.z; --z)
				getBloctAt(x, ec.y, z).RenderBlock(raster, eye, Bcoord);
	}

	if (ec.z >= 0 && ec.z <= WORLD_LENGTH - 1)
	{
		for (int y = 0; y < ec.y; ++y)
			for (int x = 0; x < ec.x; ++x)
				getBloctAt(x, y, ec.z).RenderBlock(raster, eye, Bcoord);
		for (int y = WORLD_HEIGHT - 1; y > ec.y; --y)
			for (int x = 0; x < ec.x; ++x)
				getBloctAt(x, y, ec.z).RenderBlock(raster, eye, Bcoord);
		for (int y = 0; y < ec.y; ++y)
			for (int x = WORLD_WIDTH - 1; x > ec.x; --x)
				getBloctAt(x, y, ec.z).RenderBlock(raster, eye, Bcoord);
		for (int y = WORLD_HEIGHT - 1; y > ec.y; --y)
			for (int x = WORLD_WIDTH - 1; x > ec.x; --x)
				getBloctAt(x, y, ec.z).RenderBlock(raster, eye, Bcoord);
	}

	if (ec.y >= 0 && ec.y <= WORLD_HEIGHT - 1 && ec.z >= 0 && ec.z <= WORLD_LENGTH - 1)
	{
		for (int x = 0; x < ec.x; ++x)
			getBloctAt(x, ec.y, ec.z).RenderBlock(raster, eye, Bcoord);
		for (int x = WORLD_WIDTH - 1; x > ec.x; --x)
			getBloctAt(x, ec.y, ec.z).RenderBlock(raster, eye, Bcoord);
	}

	if (ec.x >= 0 && ec.x <= WORLD_WIDTH - 1 && ec.z >= 0 && ec.z <= WORLD_LENGTH - 1)
	{
		for (int y = 0; y < ec.y; ++y)
			getBloctAt(ec.x, y, ec.z).RenderBlock(raster, eye, Bcoord);
		for (int y = WORLD_HEIGHT - 1; y > ec.y; --y)
			getBloctAt(ec.x, y, ec.z).RenderBlock(raster, eye, Bcoord);
	}

	if (ec.y >= 0 && ec.y <= WORLD_HEIGHT - 1 && ec.x >= 0 && ec.x <= WORLD_WIDTH - 1)
	{
		for (int z = 0; z < ec.z; ++z)
			getBloctAt(ec.x, ec.y, z).RenderBlock(raster, eye, Bcoord);
		for (int z = WORLD_LENGTH - 1; z > ec.z; --z)
			getBloctAt(ec.x, ec.y, z).RenderBlock(raster, eye, Bcoord);
	}
}

//-------------------- 读档辅助

void CONSOLEWORLD::World::set(Block * AllBlock) {
	for (int x = 0; x < WORLD_WIDTH; ++x)
		for (int y = 0; y < WORLD_HEIGHT; ++y)
			for (int z = 0; z < WORLD_LENGTH; ++z) {
				int ID = _blocks[x][y][z].getID();
				if (ID != 0)setBlockAt(int3(x, y, z), AllBlock[ID]);
			}
}

//-------------------- 存档函数

bool CONSOLEWORLD::World::Save(char * szPath, int num) {
	std::ofstream fout;
	char PATH[1024];
	sprintf(PATH, "%s/save/%d.txt", szPath, num);
	fout.open(PATH);
	if (!fout.is_open())return false;
	else {
		for (int x = 0; x < WORLD_WIDTH; ++x)
			for (int y = 0; y < WORLD_HEIGHT; y++)
				for (int z = 0; z < WORLD_LENGTH; ++z) {
					fout << _blocks[x][y][z].getID();
				}
		fout.close();
		fout.clear();
		return true;
	}
}

//-------------------- 读档函数

bool CONSOLEWORLD::World::Load(char * szPath, int num) {
	std::ifstream fin;
	char PATH[1024];
	sprintf(PATH, "%s/save/%d.txt", szPath, num);
	fin.open(PATH);
	if (!fin.is_open()) { std::abort(); return false; }
	else {
		for (int x = 0; x < WORLD_WIDTH; ++x)
			for (int y = 0; y < WORLD_HEIGHT; y++)
				for (int z = 0; z < WORLD_LENGTH; ++z) {
					char ch;
					fin >> ch;
					int temp = ch - '0';
					_blocks[x][y][z].setID(temp);
				}
		fin.close();
		fin.clear();
		return true;
	}
}

bool CONSOLEWORLD::World::ResetSave(char * szPath, int num)
{
	std::ofstream fout;
	char PATH[1024];
	sprintf(PATH, "%s/save/%d.txt", szPath, num);
	fout.open(PATH);
	if (!fout.is_open())
		return false;
	else
	{
		for (int x = 0; x < WORLD_WIDTH; ++x)
			for (int y = 0; y < WORLD_HEIGHT; y++)
				for (int z = 0; z < WORLD_LENGTH; ++z)
				{
					if (y == 0)
						fout << '2';
					else
						fout << '0';
				}
		fout.close();
		fout.clear();
		return true;
	}
}

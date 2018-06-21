#pragma once

#include "Camera.h"
#include "Block.h"

namespace CONSOLEWORLD
{
	using namespace ZMATH;

#define ROLEV 100

	//-------------------- 角色类，设定为一个50x50x50的方形人类
	class Role : public RENDERER::Camera1st
	{
	private:
		AABB3D _CollisionBox;

	public:
		int3 _coord;
		Ray _ray;

		Role()
		{
		
		}

		AABB3D & getAABB()
		{
			return _CollisionBox;
		}

		void setEye(float3 pos)
		{
			_eye = pos;
			_ray.setOrigin(_eye);
			_ray.setDirection(_dir);
			_CollisionBox.setExtents(_eye.x - ROLEV, _eye.y - 4 * ROLEV, _eye.z - ROLEV, _eye.x + ROLEV, _eye.y + ROLEV, _eye.z + ROLEV);
			_coord = int3(floor(_eye.x / BLOCK_EDGE_LENGTH), floor(_eye.y / BLOCK_EDGE_LENGTH), floor(_eye.z / BLOCK_EDGE_LENGTH));
		}

		void LookAt(float3 pos, float3 lookAt, float3 worldUp)
		{
			float3 look = normalize(lookAt - pos);
			float3 right = normalize(cross(worldUp, look));
			float3 up = cross(look, right);

			setEye(pos);
			_right = right;
			_up = up;
			_dir = look;
		}

		//消除掉前后左右移动时的Y坐标变化，符合MC操作
		void Walk(float dist)
		{
			float3 pos = _eye;
			float3 look = _dir;
			pos = pos + look * float3(dist, 0, dist);

			setEye(pos);
		}

		void Strafe(float dist)
		{
			float3 pos = _eye;
			float3 right = _right;
			pos = pos + right * float3(dist, 0, dist);

			setEye(pos);
		}

		void update()
		{
			_view = ZMATH::lookAt(_eye, _dir + _eye, _up);
		}
	};
}
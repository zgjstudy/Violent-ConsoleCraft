#ifndef CAMERA_H
#define CAMERA_H

#include "MATH.h"

namespace RENDERER
{
	using namespace ZMATH;
	//-------------------- 摄像机类
	class Camera
	{
	public:
		float3		_eye;		//眼睛位置
		float3		_target;	//看向的位置
		float3		_up;		//Y轴
		float3		_right;		//X轴

		float3		_dir;		//eye和targrt间射线的方向
		matrix4		_matView;	//观察矩阵

		matrix4		_matProj;	//投影矩阵
		matrix4		_matWorld;	//世界坐标矩阵
		float2		_viewSize;	//视口大小

		Camera(const float3& target = float3(0, 0, 0), const float3& eye = float3(0, 100, 100), const float3& right = float3(1, 0, 0))
		{
			_viewSize = float2(256, 256);
			_matView = matrix4(1);
			_matProj = matrix4(1);
			_matWorld = matrix4(1);

			_target = target;
			_eye = eye;
			_dir = normalize(_target - _eye);
			_right = right;
			_up = normalize(cross(_right, _dir));
		}

		float3 getEye() const
		{
			return _eye;
		}
		void setEye(float3 val)
		{
			_eye = val;
		}

		float3 getTarget() const
		{
			return _target;
		}
		void setTarget(float3 val)
		{
			_target = val;
		}

		float3 getRight() const
		{
			return _right;
		}
		void setRight(float3 val)
		{
			_right = val;
		}

		float3 getUp() const
		{
			return _up;
		}
		void setUp(float3 val)
		{
			_up = val;
		}

		float3 getDir() const
		{
			return _dir;
		}
		void calcDir()
		{
			_dir = _target - _eye;
			_dir = normalize(_dir);
		}

		void setViewSize(const float2& viewSize)
		{
			_viewSize = viewSize;
		}
		void setViewSize(float x, float y)
		{
			_viewSize = float2(x, y);
		}
		float2 getViewSize()
		{
			return _viewSize;
		}

		const matrix4 & getProject() const
		{
			return  _matProj;
		}
		void setProject(const matrix4& proj)
		{
			_matProj = proj;
		}

		const matrix4 & getView() const
		{
			return _matView;
		}


		//正交投影
		void ortho(float left, float right, float bottom, float top, float zNear, float zFar)
		{
			_matProj = ZMATH::ortho(left, right, bottom, top, zNear, zFar);
		}

		//透视投影
		void perspective(float fovy, float aspect, float zNear, float zFar)
		{
			_matProj = ZMATH::perspective<float>(fovy, aspect, zNear, zFar);
		}

		virtual void update()
		{
			_matView = ZMATH::lookAt(_eye, _target, _up);
		}
	};
	
	//-------------------- 第三人称摄像机类
	class Camera3rd : public Camera
	{
	public:
		Camera3rd(const float3& target = float3(0, 0, 0), const float3& eye = float3(0, 100, 100), const float3& right = float3(1, 0, 0))
			:Camera(target, eye, right)
		{};

		~Camera3rd() {};

		//世界坐标转化为窗口坐标
		bool project(const float4& world, float4& screen)
		{
			screen = (_matProj * _matView * _matWorld) * world;
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
			screen.x = screen.x * _viewSize.x;
			screen.y = _viewSize.y - (screen.y * _viewSize.y);
			return true;
		}

		float2 worldToScreen(const float3& world)
		{
			float4 worlds(world.x, world.y, world.z, 1);
			float4 screens;
			project(worlds, screens);
			return float2(screens.x, screens.y);
		}

		//窗口坐标转化为世界坐标
		float3 screenToWorld(const float2& screen)
		{
			float4 screens(screen.x, screen.y, 0, 1);
			float4 world;
			unProject(screens, world);
			return float3(world.x, world.y, world.z);
		}

		float3 screenToWorld(float x, float y)
		{
			float4 screens(x, y, 0, 1);
			float4 world;
			unProject(screens, world);
			return float3(world.x, world.y, world.z);
		}

		//窗口坐标转化为世界坐标
		bool unProject(const float4& screen, float4& world)
		{
			float4 v;
			v.x = screen.x;
			v.y = screen.y;
			v.z = screen.z;
			v.w = 1.0;

			// map from viewport to 0 - 1
			v.x = (v.x) / _viewSize.x;
			v.y = (_viewSize.y - v.y) / _viewSize.y;
			//v.y = (v.y - _viewPort.Y) / _viewPort.Height;

			// map to range -1 to 1
			v.x = v.x * 2.0f - 1.0f;
			v.y = v.y * 2.0f - 1.0f;
			v.z = v.z * 2.0f - 1.0f;

			RENDERER::matrix4 inverse = (_matProj * _matView * _matWorld).inverse();

			v = v * inverse;
			if (v.w == 0.0f)
			{
				return false;
			}
			world = v / v.w;
			return true;
		}

		Ray createRayFromScreen(int x, int y)
		{
			float4 minWorld;
			float4 maxWorld;

			float4 screen(float(x), float(y), 0, 1); //最近点
			float4 screen1(float(x), float(y), 1, 1);//最远点

			unProject(screen, minWorld);
			unProject(screen1, maxWorld);
			Ray ray;
			ray.setOrigin(float3(minWorld.x, minWorld.y, minWorld.z));

			float3  dir(maxWorld.x - minWorld.x, maxWorld.y - minWorld.y, maxWorld.z - minWorld.z);
			ray.setDirection(normalize(dir));
			return  ray;
		}

		//将摄像机的观察方向绕某个方向轴旋转一定的角度
		//改变观察者的位置，目标的位置不变化
		virtual void rotateViewY(float angle)
		{
			_dir = rotateY<float>(_dir, angle);
			_up = rotateY<float>(_up, angle);
			_right = normalize(cross(_dir, _up));
			float   len = length(_eye - _target);
			_eye = _target - _dir * len;
			_matView = RENDERER::lookAt(_eye, _target, _up);
		}

		virtual void rotateViewX(float angle)
		{
			matrix4 mat(1);
			mat.rotate(angle, _right);
			_dir = _dir * mat;
			_up = _up * mat;
			_right = normalize(cross(_dir, _up));
			float   len = length(_eye - _target);
			_eye = _target - _dir * len;
			_matView = RENDERER::lookAt(_eye, _target, _up);
		}
	};

	//-------------------- 第一人称摄像机类
	class Camera1st
	{
	protected:
		float3		_right;			//X轴  
		float3		_up;			//Y轴
		float3		_dir;			//Z轴
		float3		_eye;			//观察者位置

		float		_aspect;		//投影相关参数  
		float		_fovY;
		float		_nearZ;
		float		_farZ;

		matrix4		_view;			//视角矩阵  
		matrix4		_proj;			//投影矩阵  

	public:
		Camera1st() :
			_right(1.f, 0.f, 0.f),
			_up(0.f, 1.f, 0.f),
			_dir(0.f, 0.f, 1.f),
			_eye(0.f, 0.f, 0.f),
			_aspect(800.f / 600),
			_fovY(PI*0.25),
			_nearZ(1.f),
			_farZ(1000.f)
		{
			_view = matrix4(1.0f);
			_proj = ZMATH::perspective(_fovY, _aspect, _nearZ, _farZ);
		}

		//设置摄像机位置  
		void setEye(float x, float y, float z)
		{
			_eye = float3(x, y, z);
		}

		void setEye(float3 pos)
		{
			_eye = pos;
		}

		//获得摄像位置及朝向相关参数  
		float3 getEye() const
		{
			return _eye;
		}

		float3 getRight() const
		{
			return _right;
		}

		float3 getUp() const
		{
			return _up;
		}

		float3 getLook() const
		{
			return _dir;
		}


		//获得投影相关参数  
		float getNearZ()    const
		{
			return _nearZ;
		}

		float getFarZ()     const
		{
			return _farZ;
		}

		float getFovY()     const
		{
			return _fovY;
		}

		float getFovX()     const
		{
			return atan(_aspect * tan(_fovY * 0.5f)) * 2.0f;
		}

		float getAspect()   const
		{
			return _aspect;
		}


		//获得相关矩阵  
		matrix4 View() const
		{
			return _view;
		}

		matrix4 Projection() const
		{
			return _proj;
		}

		matrix4 ViewProjection() const
		{
			return _view * _proj;
		}


		//设置投影相关参数  
		void perspective(float fovY, float ratioAspect, float nearZ, float farZ)
		{
			_fovY = fovY;
			_aspect = ratioAspect;
			_nearZ = nearZ;
			_farZ = farZ;
		}


		//通过位置+观察点来设置视角矩阵 
		void LookAt(float3 pos, float3 lookAt, float3 worldUp)
		{
			float3 look = normalize(lookAt - pos);
			float3 right = normalize(cross(worldUp, look));
			float3 up = cross(look, right);

			_eye = pos;
			_right = right;
			_up = up;
			_dir = look;
		}


		//基本操作
		//沿观察方向前进或后退。因此本质m_pos沿_dir方向进行dist的平移。
		void Walk(float dist)
		{
			float3 pos = _eye;
			float3 look = _dir;
			pos = pos + look * float3(dist);

			_eye = pos;
		}

		//Strafe函数同理，用于左右平移，即m_pos沿_right方向进行dist的平移
		void Strafe(float dist)
		{
			float3 pos = _eye;
			float3 right = _right;
			pos = pos + right * float3(dist);

			_eye = pos;
		}

		//用于上下点头来调整视角，本质上即把_up向量沿_right向量进行angle的旋转
		void Pitch(float angle)
		{
			matrix4 rotation;
			rotation.rotate(angle, _right);

			_up = _up * rotation;
			_dir = _dir * rotation;
		}

		//用于摇头以旋转视角，本质上即所有三个坐标轴(_right, _up, _dir)沿Y轴方向进行angle的旋转
		void RotateY(float angle)
		{
			matrix4 rotation;
			rotation.rotateY(angle);

			_right =  _right * rotation;
			_up = _up * rotation;
			_dir =  _dir * rotation;
		}

		//更新相关矩阵
		virtual void update()
		{
			_view = ZMATH::lookAt(_eye, _dir + _eye, _up);
		}
	};

}

#endif // !CAMERA_H
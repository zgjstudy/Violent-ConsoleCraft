#pragma once
#include <Windows.h>
#include <tchar.h>
#include "Raster3D.h"
#include "TIME.h"
#include "Role.h"
#include "ConsoleCraftWorld.h"


inline int KeyDown(int vKey)
{
	return GetAsyncKeyState(vKey) & 0x8000;
}

using ZMATH::float2;
using ZMATH::float3;
using ZMATH::int2;
using ZMATH::int3;
using ZMATH::Rgba;
using ZMATH::matrix2;
using ZMATH::matrix3;
using ZMATH::matrix4;
using RENDERER::Vertex;

//-------------------- 生成移动速度信息
float speed = 8.0f;
float3 Xspeed = float3(speed, 0, 0);
float3 Yspeed = float3(0, speed, 0);
float3 Zspeed = float3(0, 0, speed);

//-------------------- 生成摄像机
int cameraMidX = WORLD_WIDTH * BLOCK_EDGE_LENGTH / 2;
int cameraMidZ = WORLD_LENGTH * BLOCK_EDGE_LENGTH / 2;
int cameraMidY = WORLD_HEIGHT * BLOCK_EDGE_LENGTH / 2;
CONSOLEWORLD::Role role;
CONSOLEWORLD::World BlockWorld;
CONSOLEWORLD::Block block0;

//-------------------- 生成状态
int width = 800;
int height = 600;
POINT lastPos{ width / 2, height / 2 };
int3 Bcoord;
int3 Scoord;
int selected = 2;
bool F3Mode = 0;
bool leave = 0;
bool loaded = 0;
int openedFileNum = 1;
bool autoSave = 0;

//-------------------- 生成屏幕信息
int frameX_Width = GetSystemMetrics(SM_CXFIXEDFRAME);
int frameY_Height = GetSystemMetrics(SM_CYFIXEDFRAME);
int frameY_Caption = GetSystemMetrics(SM_CYCAPTION);
int scrWidth = GetSystemMetrics(SM_CXSCREEN);
int scrHeight = GetSystemMetrics(SM_CYSCREEN);

//-------------------- 计算射线与XOZ平面的交点
float3 calcIntersectPoint(ZMATH::Ray & ray)
{
	float3 pos = ray.getOrigin();
	float tm = abs((pos.y) / ray.getDirection().y);
	float3 target = ray.getPoint(tm);
	return float3(target.x, 0, target.z);
}

//-------------------- 获取控制台句柄
HWND GetConsoleHwnd(void)
{
#define MY_BUFSIZE 1024 // 控制台标题的缓冲区大小
	HWND hwndFound; // 返回给调用者的窗口句柄
	TCHAR NewWindowTitle[MY_BUFSIZE]; // Contains fabricated 
	TCHAR OldWindowTitle[MY_BUFSIZE]; // Contains original 

	GetConsoleTitle(OldWindowTitle, MY_BUFSIZE);

	wsprintf(NewWindowTitle, L"%d/%d",
		GetTickCount(),
		GetCurrentProcessId());

	SetConsoleTitle(NewWindowTitle);

	Sleep(40);

	hwndFound = FindWindow(NULL, NewWindowTitle);

	SetConsoleTitle(OldWindowTitle);

	return(hwndFound);
}

//-------------------- 获取当前文件路径
void getResourcePath(HINSTANCE hInstance, char pPath[1024])
{
	char szPathName[1024];
	char szDriver[64];
	char szPath[1024];
	GetModuleFileNameA(hInstance, szPathName, sizeof(szPathName));
	_splitpath(szPathName, szDriver, szPath, 0, 0);
	sprintf(pPath, "%s%s", szDriver, szPath);
}

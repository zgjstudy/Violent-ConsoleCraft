#include "Running.h"

int main(HINSTANCE hInstance)
{
	//-------------------- 生成缓冲窗口句柄
	HWND hWnd = CreateWindowExA(
		NULL,
		"Violent ConsoleCraft",
		"Violent ConsoleCraft",
		WS_POPUPWINDOW,
		0,
		0,
		800,
		600,
		0,
		0,
		hInstance,
		0);

	//-------------------- 生成画板信息
	int width = 800;
	int height = 600;

	//-------------------- 生成双缓冲
	void * buffer = nullptr;
	HDC hDC = GetDC(hWnd);
	HDC hMem = ::CreateCompatibleDC(hDC);

	//-------------------- 生成画布信息
	BITMAPINFO bmpInfo;
	bmpInfo.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	bmpInfo.bmiHeader.biWidth = width;
	bmpInfo.bmiHeader.biHeight = -height;
	bmpInfo.bmiHeader.biPlanes = 1;
	bmpInfo.bmiHeader.biBitCount = 32;
	bmpInfo.bmiHeader.biCompression = BI_RGB;
	bmpInfo.bmiHeader.biSizeImage = 0;
	bmpInfo.bmiHeader.biXPelsPerMeter = 0;
	bmpInfo.bmiHeader.biYPelsPerMeter = 0;
	bmpInfo.bmiHeader.biClrUsed = 0;
	bmpInfo.bmiHeader.biClrImportant = 0;

	HBITMAP hBmp = CreateDIBSection(hDC, &bmpInfo, DIB_RGB_COLORS, (void**)&buffer, 0, 0);
	SelectObject(hMem, hBmp);

	//-------------------- 配置控制台窗口
	SetConsoleTitle(L"Violent ConsoleCraft");
	HWND consoleHwnd = GetConsoleHwnd();
	void * consoleBuffer = nullptr;
	HDC consoleDC = GetDC(consoleHwnd);
	HDC consoleMem = ::CreateCompatibleDC(consoleDC);
	HANDLE consoleOut = ::GetStdHandle(STD_OUTPUT_HANDLE);
	HANDLE consoleIn = ::GetStdHandle(STD_INPUT_HANDLE);
	CONSOLE_CURSOR_INFO CursorInfo;
	DWORD state = 0, res;
	INPUT_RECORD rec;

	//-------------------- 加载纹理
	char szPath[1024];
	getResourcePath(0, szPath);
	char szImage[1024];

	//planks_birch
	sprintf(szImage, "%s/blocks/planks_birch.png", szPath);
	RENDERER::Image * planks_birch = RENDERER::Image::loadFromFile(szImage);

	//melon_top
	sprintf(szImage, "%s/blocks/melon_top.png", szPath);
	RENDERER::Image * melon_top = RENDERER::Image::loadFromFile(szImage);

	//melon_side
	sprintf(szImage, "%s/blocks/melon_side.png", szPath);
	RENDERER::Image * melon_side = RENDERER::Image::loadFromFile(szImage);

	//cobblestone
	sprintf(szImage, "%s/blocks/cobblestone.png", szPath);
	RENDERER::Image * cobblestone = RENDERER::Image::loadFromFile(szImage);

	//brick
	sprintf(szImage, "%s/blocks/brick.png", szPath);
	RENDERER::Image * brick = RENDERER::Image::loadFromFile(szImage);

	//grass_side
	sprintf(szImage, "%s/blocks/grass_side.png", szPath);
	RENDERER::Image * grass_side = RENDERER::Image::loadFromFile(szImage);

	//grass_top
	sprintf(szImage, "%s/blocks/grass_top.png", szPath);
	RENDERER::Image * grass_top = RENDERER::Image::loadFromFile(szImage);

	//dirt
	sprintf(szImage, "%s/blocks/dirt.png", szPath);
	RENDERER::Image * dirt = RENDERER::Image::loadFromFile(szImage);

	//hay_block_top
	sprintf(szImage, "%s/blocks/hay_block_top.png", szPath);
	RENDERER::Image * hay_block_top = RENDERER::Image::loadFromFile(szImage);

	//hay_block_side
	sprintf(szImage, "%s/blocks/hay_block_side.png", szPath);
	RENDERER::Image * hay_block_side = RENDERER::Image::loadFromFile(szImage);

	//stonebrick
	sprintf(szImage, "%s/blocks/stonebrick.png", szPath);
	RENDERER::Image * stonebrick = RENDERER::Image::loadFromFile(szImage);

	//leaves_oak
	sprintf(szImage, "%s/blocks/leaves_oak.png", szPath);
	RENDERER::Image * leaves_oak = RENDERER::Image::loadFromFile(szImage);

	//log_oak
	sprintf(szImage, "%s/blocks/log_oak.png", szPath);
	RENDERER::Image * log_oak = RENDERER::Image::loadFromFile(szImage);

	//log_oak_top
	sprintf(szImage, "%s/blocks/log_oak_top.png", szPath);
	RENDERER::Image * log_oak_top = RENDERER::Image::loadFromFile(szImage);

	//UI
	sprintf(szImage, "%s/image/bottom.png", szPath);
	RENDERER::Image * Pbottom = RENDERER::Image::loadFromFile(szImage);

	sprintf(szImage, "%s/image/widgets.png", szPath);
	RENDERER::Image * widgets = RENDERER::Image::loadFromFile(szImage);

	sprintf(szImage, "%s/image/cursor.png", szPath);
	RENDERER::Image * cursor = RENDERER::Image::loadFromFile(szImage);

	sprintf(szImage, "%s/image/First.jpg", szPath);
	RENDERER::Image * SMBG = RENDERER::Image::loadFromFile(szImage);

	sprintf(szImage, "%s/image/Second.jpg", szPath);
	RENDERER::Image * SSBG = RENDERER::Image::loadFromFile(szImage);

	sprintf(szImage, "%s/image/Press Enter To Play.png", szPath);
	RENDERER::Image * PETP = RENDERER::Image::loadFromFile(szImage);
	

	//Debug
	sprintf(szImage, "%s/image/555.jpg", szPath);
	RENDERER::Image * image = RENDERER::Image::loadFromFile(szImage);

	sprintf(szImage, "%s/blocks/debug.png", szPath);
	RENDERER::Image * debug = RENDERER::Image::loadFromFile(szImage);


	//-------------------- 构造方块
	CONSOLEWORLD::Block AllBlock[10];
	int i = 1;
	CONSOLEWORLD::Block B_planks_birch;
	B_planks_birch.setID(1);
	B_planks_birch.setTopImage(planks_birch);
	B_planks_birch.setSideImage(planks_birch);
	B_planks_birch.setBottomImage(planks_birch);
	AllBlock[i] = B_planks_birch;  ++i;

	CONSOLEWORLD::Block B_grass;
	B_grass.setID(2);
	B_grass.setTopImage(grass_top);
	B_grass.setSideImage(grass_side);
	B_grass.setBottomImage(dirt);
	AllBlock[i] = B_grass;  ++i;

	CONSOLEWORLD::Block B_cobblestone;
	B_cobblestone.setID(3);
	B_cobblestone.setTopImage(cobblestone);
	B_cobblestone.setSideImage(cobblestone);
	B_cobblestone.setBottomImage(cobblestone);
	AllBlock[i] = B_cobblestone;  ++i;

	CONSOLEWORLD::Block B_melon;
	B_melon.setID(4);
	B_melon.setTopImage(melon_top);
	B_melon.setSideImage(melon_side);
	B_melon.setBottomImage(melon_top);
	AllBlock[i] = B_melon;  ++i;


	CONSOLEWORLD::Block B_brick;
	B_brick.setID(5);
	B_brick.setTopImage(brick);
	B_brick.setSideImage(brick);
	B_brick.setBottomImage(brick);
	AllBlock[i] = B_brick;  ++i;


	CONSOLEWORLD::Block B_stonebrick;
	B_stonebrick.setID(6);
	B_stonebrick.setTopImage(stonebrick);
	B_stonebrick.setSideImage(stonebrick);
	B_stonebrick.setBottomImage(stonebrick);
	AllBlock[i] = B_stonebrick;  ++i;

	CONSOLEWORLD::Block B_hay_block_side;
	B_hay_block_side.setID(7);
	B_hay_block_side.setTopImage(hay_block_top);
	B_hay_block_side.setSideImage(hay_block_side);
	B_hay_block_side.setBottomImage(hay_block_top);
	AllBlock[i] = B_hay_block_side;  ++i;

	CONSOLEWORLD::Block B_leaves_oak;
	B_leaves_oak.setID(8);
	B_leaves_oak.setTopImage(leaves_oak);
	B_leaves_oak.setSideImage(leaves_oak);
	B_leaves_oak.setBottomImage(leaves_oak);
	AllBlock[i] = B_leaves_oak;  ++i;

	CONSOLEWORLD::Block B_log_oak;
	B_log_oak.setID(9);
	B_log_oak.setTopImage(log_oak_top);
	B_log_oak.setSideImage(log_oak);
	B_log_oak.setBottomImage(log_oak_top);
	AllBlock[i] = B_log_oak;  ++i;

	block0 = B_grass;

	//-------------------- 生成渲染器
	RENDERER::Raster3D raster(width, height, buffer);

	raster.setViewPort(0, 0, width, height);
	raster.setPerspective(60, (float)(width) / (float)(height), 0.1, 10000);
	raster.ignoreOutside();

	role.LookAt(float3(cameraMidX, cameraMidY, 1500), float3(cameraMidX, cameraMidY, cameraMidZ), float3(0, 1, 0));
	role.perspective(45, (float)(width) / (float)(height), 1.f, 1000.f);
	role.update();

	enum GAME_STATE
	{
		GS_STARTMENU = 0,
		GS_SELECTSAVE,
		GS_CRAFT,
		GS_LEAVE
	} gamestate = GS_STARTMENU;


	//-------------------- 一个迷你状态机
	while (true)
	{
		::MoveWindow(
			consoleHwnd,
			(scrWidth - width) / 2,
			(scrHeight - height) / 2,
			width * 0.845,
			height * 0.87,
			true);

		switch (gamestate)
		{
		case GS_STARTMENU:
		{
			raster.drawImage(0, 0, SMBG);

			static float alphaBlend = 0;
			raster.drawImage(0, 0, SMBG);
			raster.drawImage(width / 2 + 42, height / 2 + 50, PETP, sin(alphaBlend += 0.01));
			BitBlt(consoleDC, 0, 0, width, height, hMem, 0, 0, SRCCOPY);
			if (KeyDown(VK_RETURN))
			{
				gamestate = GS_SELECTSAVE;
			}
			if (KeyDown(VK_ESCAPE) && KeyDown('1'))	// 按ESC + 1键退出程序
			{
				gamestate = GS_LEAVE;
				leave = 1;
				break;
			}
			break;
		}
		case GS_SELECTSAVE:
		{
			if (loaded == 1)
			{

				BlockWorld.Load(szPath, openedFileNum);
				BlockWorld.set(AllBlock);

				//负载测试
				//for (int i = 0; i < WORLD_HEIGHT; ++i)
				//{
				//	BlockWorld.setY(i, B_brick);
				//}
				

				gamestate = GS_CRAFT;
			}

			if (KeyDown('R'))
			{
				if (KeyDown('1'))
				{
					BlockWorld.ResetSave(szPath, 1);
					MessageBoxA(0, "Success!", "提示", 0);
				}
				if (KeyDown('2'))
				{
					BlockWorld.ResetSave(szPath, 2);
					MessageBoxA(0, "Success!", "提示", 0);
				}
				if (KeyDown('3'))
				{
					BlockWorld.ResetSave(szPath, 3);
					MessageBoxA(0, "Success!", "提示", 0);
				}
			}
			if (KeyDown(VK_LSHIFT))
			{
				if (KeyDown('1'))
				{
					openedFileNum = 1;
					loaded = 1;
				}
				if (KeyDown('2'))
				{
					openedFileNum = 2;
					loaded = 1;
				}
				if (KeyDown('3'))
				{
					openedFileNum = 3;
					loaded = 1;
				}
			}
			raster.drawImage(0, 0, SSBG);
			BitBlt(consoleDC, 0, 0, width, height, hMem, 0, 0, SRCCOPY);
			break;
		}
		case GS_CRAFT:
		{
			loaded = 0;
			autoSave = 0;

			//-------------------- 消息循环
			while (true)
			{
				::MoveWindow(
					consoleHwnd,
					(scrWidth - width) / 2,
					(scrHeight - height) / 2,
					width * 0.845,
					height * 0.87,
					true);

				//-------------------- Console消息获取
				ReadConsoleInput(consoleIn, &rec, 1, &res);

				if (KeyDown('A'))
				{
					role.Strafe(speed);
					if (!BlockWorld.MoveCollideCheckX(role.getAABB()) || !BlockWorld.MoveCollideCheckZ(role.getAABB()))
						role.Strafe(-speed);
				}
				else if (KeyDown('D'))
				{
					role.Strafe(-speed);
					if (!BlockWorld.MoveCollideCheckX(role.getAABB()) || !BlockWorld.MoveCollideCheckZ(role.getAABB()))
						role.Strafe(speed);
				}
				if (KeyDown('W'))
				{
					role.Walk(speed);
					if (!BlockWorld.MoveCollideCheckX(role.getAABB()) || !BlockWorld.MoveCollideCheckZ(role.getAABB()))
						role.Walk(-speed);
				}
				else if (KeyDown('S'))
				{
					role.Walk(-speed);
					if (!BlockWorld.MoveCollideCheckX(role.getAABB()) || !BlockWorld.MoveCollideCheckZ(role.getAABB()))
						role.Walk(speed);
				}
				if (KeyDown(VK_SPACE))
				{
					role.setEye(role.getEye() + Yspeed);
					if (!BlockWorld.MoveCollideCheckY(role.getAABB()))
						role.setEye(role.getEye() - Yspeed);
				}
				else if (KeyDown(VK_LSHIFT))
				{
					role.setEye(role.getEye() - Yspeed);
					if (!BlockWorld.MoveCollideCheckY(role.getAABB()))
						role.setEye(role.getEye() + Yspeed);
				}
				if (KeyDown(VK_F3))
				{
					F3Mode = !F3Mode;
				}
				if (KeyDown('G'))
				{
					BlockWorld.Save(szPath, openedFileNum);
				}
				if (KeyDown('1'))
				{
					selected = 1;
					block0 = B_planks_birch;
				}
				if (KeyDown('2'))
				{
					selected = 2;
					block0 = B_grass;
				}
				if (KeyDown('3'))
				{
					selected = 3;
					block0 = B_cobblestone;
				}
				if (KeyDown('4'))
				{
					selected = 4;
					block0 = B_melon;
				}
				if (KeyDown('5'))
				{
					selected = 5;
					block0 = B_brick;
				}
				if (KeyDown('6'))
				{
					selected = 6;
					block0 = B_stonebrick;
				}
				if (KeyDown('7'))
				{
					selected = 7;
					block0 = B_hay_block_side;
				}
				if (KeyDown('8'))
				{
					selected = 8;
					block0 = B_leaves_oak;
				}
				if (KeyDown('9'))
				{
					selected = 9;
					block0 = B_log_oak;
				}
				if (KeyDown('H'))
				{
					autoSave = 1;
				}
				if (KeyDown('J'))
				{
					autoSave = 0;
				}
				if (KeyDown(VK_ESCAPE))	// 按ESC键退出程序
				{
					gamestate = GS_STARTMENU;
					break;
				}
				if (KeyDown(VK_LCONTROL))
				{
					speed = 20.0f;
					Xspeed = float3(speed, 0, 0);
					Yspeed = float3(0, speed, 0);
					Zspeed = float3(0, 0, speed);
				}
				else
				{
					speed = 8.0f;
					Xspeed = float3(speed, 0, 0);
					Yspeed = float3(0, speed, 0);
					Zspeed = float3(0, 0, speed);
				}
				if (rec.EventType == MOUSE_EVENT)
				{
					if (rec.Event.MouseEvent.dwEventFlags == MOUSE_MOVED)
					{
						POINT currentPos;
						::GetCursorPos(&currentPos);

						float dx = currentPos.x - lastPos.x;
						float dy = currentPos.y - lastPos.y;

						role.Pitch(dy);
						role.RotateY(-dx);

						lastPos = currentPos;

						role.update();
					}
					if (rec.Event.MouseEvent.dwButtonState == FROM_LEFT_1ST_BUTTON_PRESSED)
					{
						if (Bcoord != int3(-1, -1, -1))
						{
							BlockWorld.getBloctAt(Bcoord).Break();
						}
					}
					if (rec.Event.MouseEvent.dwButtonState == RIGHTMOST_BUTTON_PRESSED)
					{
						if (Scoord != int3(-1, -1, -1))
						{
							if(Scoord != role._coord)
								BlockWorld.setBlockAt(Scoord, block0);
						}
					}
				}

				GetConsoleCursorInfo(consoleOut, &CursorInfo);	//获取控制台光标信息
				CursorInfo.bVisible = false;					//隐藏控制台光标
				SetConsoleCursorInfo(consoleOut, &CursorInfo);	//设置控制台光标状态

				//-------------------- 刷新摄像机
				role.update();
				Bcoord = BlockWorld.GetPointedBlock(role._ray);
				Scoord = BlockWorld.GetSelectedPosition(role._ray, Bcoord);

				//-------------------- 刷新帧率计数器
				RENDERER::Timestamp tms;
				tms.update();
				raster.clear();

				raster.setView(role.View());

				//-------------------- 渲染世界
				BlockWorld.RenderWorld(raster, role.getEye(), Bcoord);
				raster.drawImage(width / 2 - 4, height / 2 - 4, cursor);
				raster.drawImage(100, 525, Pbottom);
				raster.drawImage(100 + 66 * (selected - 1), 525, widgets);

				//-------------------- 计算并输出FPS
				double ms = tms.getElapsedTimeInMilliSec();
				TCHAR szBuf[128];
				swprintf_s(szBuf, L"FPS: %.2f", 1000 / ms);
				if (F3Mode)
					TextOut(hMem, 0, 0, szBuf, wcslen(szBuf));
				swprintf_s(szBuf, L"(%d,%d,%d)", role._coord.x, role._coord.y, role._coord.z);
				if (F3Mode)
					TextOut(hMem, 0, 20, szBuf, wcslen(szBuf));
				swprintf_s(szBuf, L"(%d,%d,%d)", Bcoord.x, Bcoord.y, Bcoord.z);
				if (F3Mode)
					TextOut(hMem, 0, 40, szBuf, wcslen(szBuf));

				//-------------------- 拷贝图像到DC
				BitBlt(consoleDC, 0, 0, width, height, hMem, 0, 0, SRCCOPY);

				//-------------------- 刷新深度缓冲区
				raster.initializeZbuffer();
			}

			if(autoSave)
				BlockWorld.Save(szPath, openedFileNum);

			break;
		}
		default:
			break;
		}
		
		if (leave)
		{
			//-------------------- 释放内存
			delete planks_birch;
			delete melon_top;
			delete melon_side;
			delete cobblestone;
			delete brick;
			delete grass_side;
			delete grass_top;
			delete dirt;
			delete hay_block_top;
			delete hay_block_side;
			delete stonebrick;
			delete leaves_oak;
			delete log_oak;
			delete log_oak_top;
			delete Pbottom;
			delete widgets;
			delete cursor;

			break;
		}
	}
	return 0;
}
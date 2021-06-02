#include<iostream>
#include <map>

#include"Tool.h"
#include"ToolMesh.h"


using namespace std;
using namespace MeshLib;



void main(int argc, char** argv)
{
	CTMesh mesh;
	CTMesh mesh2;
	if(argc<=1)
	{
		exit(EXIT_SUCCESS);
	}
	mesh.read_m(argv[1]);
	//CTool<CTMesh,CToolVertex> tool(&mesh);
	//Gauss-Bonet Theory
	//tool.GaussBonetTheory();
	//HarmonyMapping<CTMesh,CToolVertex> mappingTool(&mesh);
	//mappingTool.StartMapping();
	
	//复变函数积分
	//运行命令:
	//E:
	//cd \Projects\PolyCube\PolyCube\zzzbin\basic_example
	//g3dogl.exe Out.m -texturemap check.bmp
	//ComplexIntegral(&mesh);
	ComplexIntegral<CTMesh,CToolVertex> intergralMesh;
	intergralMesh.Integral(&mesh);
	mesh.write_m(argv[2]);
	cout << "finished !!! press any key to continue!!!" << endl;
	getchar();
}
#include<iostream>
#include <map>

#include"Tool.h"
#include"ToolMesh.h"
#include <opencv2/opencv.hpp>

#define PICTURE_NAME "blocks.png"


using namespace std;
using namespace MeshLib;
using namespace cv;

string GetExecPath(string raw_argv) //��ȡ��ִ���ļ�Ŀ¼
{
	string res;
	for (int i = raw_argv.length()-1; i > 0; --i)
	{
		if(raw_argv[i]!='\\')
		{
			res = raw_argv.substr(0,i);
		}
		else
		{
			return res;
		}
	}
	return res;
}


void main(int argc, char** argv)
{
	CTMesh mesh;
	CTMesh mesh2;
	if(argc<=1)
	{
		exit(EXIT_SUCCESS);
	}
	mesh.read_m(argv[1]);
	cout<<GetExecPath(argv[0]);
	Mat test = imread(GetExecPath(argv[0])+PICTURE_NAME); //����ͼ��Mat
    imshow("test", test);
    waitKey(0);
	//CTool<CTMesh,CToolVertex> tool(&mesh);
	//Gauss-Bonet Theory
	//tool.GaussBonetTheory();
	//HarmonyMapping<CTMesh,CToolVertex> mappingTool(&mesh);
	//mappingTool.StartMapping();
	
	//���亯������
	//��������:
	//E:
	//cd \Projects\PolyCube\PolyCube\zzzbin\basic_example
	//g3dogl.exe Out.m -texturemap check.bmp
	//ComplexIntegral(&mesh);
	//ComplexIntegral<CTMesh,CToolVertex> intergralMesh;
	//intergralMesh.Integral(&mesh);
	mesh.write_m(argv[2]);
	cout << "finished !!! press any key to continue!!!" << endl;
	getchar();
}
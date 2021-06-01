#include<iostream>
#include <map>

#include"Tool.h"
#include"ToolMesh.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace MeshLib;
using namespace Eigen;

map<int,int>* pVMap;
#define VMAP (*pVMap)

Eigen::MatrixXcd ComplexFunc(Eigen::MatrixXcd z)
{
	Eigen::MatrixXcd res = Eigen::MatrixXcd::Ones(1,1);
	Eigen::MatrixXcd z0 = Eigen::MatrixXcd::Ones(1,1);
	Eigen::MatrixXcd z1 = Eigen::MatrixXcd::Ones(1,1);
	z0.real() << 1;
	z0.imag() << 2;
	z1.real() << 2;
	z1.imag() << 4;
	//res = (z-z1).conjugate();
	//res(0) /= pow( abs(z(0)),2);

	//res = (z-z0)*res(0);
	res = (z-z0)*(z-z1);
	return res*0.01;


	
	//return z*0.01;
}

Eigen::MatrixXcd CompNum(CPoint cpoint)
{
	Eigen::MatrixXcd point = Eigen::MatrixXcd::Zero(1,1);
	point.real() << cpoint[0];
	point.imag() << cpoint[1];
	return point;
}

CPoint2 getPoint2(Eigen::MatrixXcd c)
{
	return CPoint2(c.real()(0,0),c.imag()(0,0));
}

CPoint2 getDeltaIntegral(CToolVertex* v1,CToolVertex* v2)
{
	return getPoint2((CompNum(v1->point())-CompNum(v2->point()))
		*ComplexFunc((CompNum(v1->point()/2)+CompNum(v2->point()/2))));
}


void IntegralRound(CToolVertex* vertex)
{
	bool flag = false;
	for (VertexVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> vertex_iterator(vertex);
		!vertex_iterator.end();
		vertex_iterator++)
	{
		if(!vertex_iterator.value()->isInteg)
		{
			vertex_iterator.value()->uv() = getDeltaIntegral(vertex_iterator.value(),vertex)+vertex->uv();
			vertex_iterator.value()->isInteg = true;
			flag = true;
		}
	}
	if(!flag)
	{
		return;
	}
	for (VertexVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> vertex_iterator(vertex);
		!vertex_iterator.end();
		vertex_iterator++)
	{
		IntegralRound(vertex_iterator.value());
	}
}

void ComplexIntegral(CTMesh* mesh)
{
	CPoint2 zero(0,0);
	CToolVertex* zeroPoint = NULL;
	for (MeshVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> vertex_iterator(mesh);
		!vertex_iterator.end();
		vertex_iterator++)
	{
		double distance = 1.0;
		if(distance>vertex_iterator.value()->point().norm())
		{
			distance = vertex_iterator.value()->point().norm();
			zeroPoint = vertex_iterator.value();
		}
	}
	zeroPoint->uv() = zero;
	zeroPoint->isInteg=true;
	IntegralRound(zeroPoint);
}

CHalfEdge* FindNextBoundaryHE(CHalfEdge* he)
{
	he = he->he_next();
	while (!he->edge()->boundary())
	{
		he = he->he_sym()->he_next();
	}
	return he;
}

void SetUV_CircleBorder(CTMesh* mesh)
{
	int boundarySum = 0;
	CHalfEdge* initHE = NULL;
	for(MeshHalfEdgeIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> he_iterator(mesh);
		!he_iterator.end();
		++he_iterator)
	{
		if(he_iterator.value()->vertex()->boundary()&&he_iterator.value()->source()->boundary())
		{
			initHE = he_iterator.value();
			break;
		}
	}
	CHalfEdge* half_edge = initHE;
	CHalfEdge* half_iter = FindNextBoundaryHE(half_edge);
	for (;half_iter!=half_edge;half_iter = FindNextBoundaryHE(half_iter))
	{
		boundarySum++;
	}
	half_iter = FindNextBoundaryHE(half_iter);
	boundarySum++;
	auto boundaryNorm = new double[boundarySum];
	for (int i = 0; i<boundarySum;i++)
	{
		boundaryNorm[i] = half_iter->source()->point() * half_iter->vertex()->point();
		if(i!=0)
		{
			boundaryNorm[i] += boundaryNorm[i - 1];
		}
		half_iter = FindNextBoundaryHE(half_iter);
	}
	for (int i = 0; i < boundarySum; i++)
	{
		half_iter->vertex()->uv() = CPoint2(
			cos(2 * PI * boundaryNorm[i] / boundaryNorm[boundarySum - 1]),
			sin(2 * PI * boundaryNorm[i] / boundaryNorm[boundarySum - 1])
		);
		half_iter = FindNextBoundaryHE(half_iter);
	}
}

double GetVertexWeight(CToolVertex* v1, CToolVertex* v2,CTMesh* mesh)
{
	double res = 0.0;
	CHalfEdge* half_edge = nullptr;
	
	for (VertexOutHalfedgeIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> 
		half_edge_iter(mesh,v1);
		!half_edge_iter.end();
		++half_edge_iter)
	{
		
		if(half_edge_iter.value()->vertex()->id()==v2->id())
		{
			half_edge = half_edge_iter.value();
		}
	}
	CVertex* v3 = half_edge->he_next()->vertex();
	res+= (v1->point() - v3->point()) * (v2->point() - v3->point())
		/ ((v1->point() - v3->point()) ^ (v2->point() - v3->point())).norm();
	if (v1->boundary() && v2->boundary())
	{
		return res;
	}
	half_edge = half_edge->he_sym();
	CVertex* v4 = half_edge->he_next()->vertex();
	res+= (v1->point() - v4->point()) * (v2->point() - v4->point())
		/ ((v1->point() - v4->point()) ^ (v2->point() - v4->point())).norm();
	return res;
}

void Set_InsideUV(CToolMesh<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge>* mesh)
{
	const int numVertices = mesh->numVertices()+1;
	Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_x = Eigen::MatrixXd::Zero(numVertices, numVertices);
	Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_y = Eigen::MatrixXd::Zero(numVertices, numVertices);
	Matrix<double, Eigen::Dynamic, 1> b_x = Eigen::MatrixXd::Zero(numVertices, 1);
	Matrix<double, Eigen::Dynamic, 1> b_y = Eigen::MatrixXd::Zero(numVertices, 1);
	Matrix<double, Eigen::Dynamic, 1> res_x = Eigen::MatrixXd::Zero(numVertices, 1);
	Matrix<double, Eigen::Dynamic, 1> res_y = Eigen::MatrixXd::Zero(numVertices, 1);
	for(MeshVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> vertex_iterator(mesh);
		!vertex_iterator.end();
		++vertex_iterator)
	{
		if (vertex_iterator.value()->boundary())
		{
			b_x(VMAP[vertex_iterator.value()->id()]) = vertex_iterator.value()->uv()[0];
			b_y(VMAP[vertex_iterator.value()->id()]) = vertex_iterator.value()->uv()[1];
			
			matrix_x(VMAP[vertex_iterator.value()->id()], VMAP[vertex_iterator.value()->id()]) = 1;
			matrix_y(VMAP[vertex_iterator.value()->id()], VMAP[vertex_iterator.value()->id()]) = 1;
			continue;
		}
		CToolVertex* vertex = vertex_iterator.value();
		double totalWeight = 0.0;
		for (VertexVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> vertex_iter(vertex);
			!vertex_iter.end();
			vertex_iter++)
		{
			totalWeight += GetVertexWeight(vertex_iter.value(), vertex,mesh);
		}
		matrix_x(VMAP[vertex_iterator.value()->id()], VMAP[vertex_iterator.value()->id()]) = -totalWeight;
		matrix_y(VMAP[vertex_iterator.value()->id()], VMAP[vertex_iterator.value()->id()]) = -totalWeight;
		for (VertexVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> vertex_iter(vertex);
			!vertex_iter.end();
			vertex_iter++)
		{
			matrix_x(VMAP[vertex->id()], VMAP[vertex_iter.value()->id()]) = GetVertexWeight(vertex_iter.value(),vertex, mesh);
			matrix_y(VMAP[vertex->id()], VMAP[vertex_iter.value()->id()]) = GetVertexWeight(vertex_iter.value(),vertex, mesh);
		}
	}
	res_x = matrix_x.colPivHouseholderQr().solve(b_x);
	res_y = matrix_y.colPivHouseholderQr().solve(b_y);
	for (MeshVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> vertex_iterator(mesh);
		!vertex_iterator.end();
		++vertex_iterator)
	{
		vertex_iterator.value()->uv()[0] = res_x[VMAP[vertex_iterator.value()->id()]];
		vertex_iterator.value()->uv()[1] = res_y[VMAP[vertex_iterator.value()->id()]];
		
	}
}

double SumOfAngles(CToolVertex* pV)
{
	CEdge* edge1 = NULL;
	CEdge* edge2 = NULL;
	CVertex* v1 = NULL;
	CVertex* v2 = NULL;
	double angleSum = 0.0;
	for (VertexFaceIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> faceIter(pV); !faceIter.end(); faceIter++)
	{
		CHalfEdge* he = faceIter.value()->halfedge();
		for (int i = 0; i < 3; i++)
		{
			if (he->vertex() != pV)
			{
				if (v1 == NULL)
				{
					v1 = he->vertex();
				}
				else
				{
					v2 = he->vertex();
				}
			}
			he = he->he_next();
		}
		angleSum+= acos(
			((v1->point() - pV->point()) * (v2->point() - pV->point()))
			/ ((v1->point() - pV->point()).norm() * (v2->point() - pV->point()).norm())
		);

		edge1 = NULL;
		edge2 = NULL;
		v1 = NULL;
		v2 = NULL;
	}
	return angleSum;
}

void SetTransfer(CTMesh* mesh)
{
	int i = 1;
	pVMap = new std::map<int,int>();
	for(MeshVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> vertex_iterator(mesh);
		!vertex_iterator.end();
		vertex_iterator++)
	{
		VMAP[vertex_iterator.value()->id()]=i++;
	}
}


void main(int argc, char** argv)
{
	cout << "Test:hello word!" << endl;

	CTMesh mesh;
	CTMesh mesh2;
	if(argc<=1)
	{
		exit(EXIT_SUCCESS);
	}
	mesh.read_m(argv[1]);
	CTool<CTMesh> tool(&mesh);
	//VMap = new std::map<int,int>();
	//test nextHE function
	/*CHalfEdge* initHE = NULL;
	for(MeshHalfEdgeIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> he_iterator(&mesh);
		!he_iterator.end();
		++he_iterator)
	{
		if(he_iterator.value()->vertex()->boundary()&&he_iterator.value()->source()->boundary())
		{
			initHE = he_iterator.value();
			break;
		}
	}
	cout<<initHE->vertex()->id()<<endl;
	CHalfEdge* he = FindNextBoundaryHE(initHE);
	while (he->vertex()->id()!=initHE->vertex()->id())
	{
		cout<<he->vertex()->id()<<endl;
		he = FindNextBoundaryHE(he);
	}*/

	//Gauss-Bonet Theory
	/*double res = 0.0;
	for (MeshVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> v(&mesh); !v.end(); v++)
	{
		if (!v.value()->boundary())
		{
			res += PI;
		}
		res += PI - SumOfAngles(v.value());
	}
	cout << "Gauss-Bonet Theory:" << endl;
	if (abs(res - 2 * PI * (mesh.numFaces() - mesh.numEdges() + mesh.numVertices())) > 0.001)
	{
		cout << "Not Match" << endl;
		cout << "Left Value " << res << endl;
		cout << "Right Value " << 2 * PI * (mesh.numFaces() - mesh.numEdges() + mesh.numVertices()) << endl;
	}
	else
	{
		cout << "Match!" << endl;
		cout << "Left Value " << res << endl;
		cout << "Right Value " << 2 * PI * (mesh.numFaces() - mesh.numEdges() + mesh.numVertices()) << endl;
	}*/

	//调和映照
	SetTransfer(&mesh);
	SetUV_CircleBorder(&mesh);
	Set_InsideUV(&mesh);
	mesh.write_m(argv[2]);

	//复变函数积分
	//运行命令:
	//E:
	//cd \Projects\PolyCube\PolyCube\zzzbin\basic_example
	//g3dogl.exe Out.m -texturemap check.bmp
	ComplexIntegral(&mesh);
	mesh.write_m(argv[2]);
	cout << "finished !!! press any key to continue!!!" << endl;
	
	getchar();
}
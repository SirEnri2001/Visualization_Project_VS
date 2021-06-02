
#pragma once

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <map>


#ifndef VMAP
#define VMAP (*pVMap)
#endif

#ifndef PI
#define PI 3.1415926535
#endif

namespace MeshLib
{
	using namespace std;
	using namespace Eigen;

	template<typename M,typename V>
	class HarmonyMapping
	{
		map<int,int>* pVMap = NULL; //矩阵下标和VertexID构成的映射
		M* p_mesh;//网格
		void SetTransfer(M* mesh);//初始化矩阵下标和VertexID的映射
		CHalfEdge* FindNextBoundaryHE(CHalfEdge* he);//找到边界半边所指向的下一个半边
	public:
		void Set_InsideUV(M* mesh);//边界映射设定好的前提下，设置内部的映射
		void SetUV_CircleBorder(M* mesh);//设定边界映射为圆形的简单映射
		double GetVertexWeight(V* v1, V* v2,M* mesh);//获得两点间的权重
		void StartMapping();//调和映照启动函数
		HarmonyMapping(M* pMesh);
		~HarmonyMapping()
		{
			delete pVMap;
		}
		
	};

	template <typename M,typename V>
	void HarmonyMapping<M,V>::StartMapping(){
		SetTransfer(p_mesh);
		SetUV_CircleBorder(p_mesh);
		Set_InsideUV(p_mesh);
	}

	template <typename M,typename V>
	HarmonyMapping<M,V>::HarmonyMapping(M* pMesh)
	{
		p_mesh = pMesh;
	}

	
	template<typename M,typename V>
	void HarmonyMapping<M,V>::SetTransfer(M* mesh)
	{
		int i = 1;
		pVMap = new std::map<int,int>();
		for(MeshVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> vertex_iterator(mesh);
			!vertex_iterator.end();
			vertex_iterator++)
		{
			VMAP[vertex_iterator.value()->id()]=i++; //根据网格内点遍历顺序构建map
		}
	}


	template<typename M,typename V>
	void HarmonyMapping<M,V>::Set_InsideUV(M* mesh)
	{
		const int numVertices = mesh->numVertices()+1; //矩阵总维数
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_x = Eigen::MatrixXd::Zero(numVertices, numVertices); //映射u坐标方程系数矩阵
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_y = Eigen::MatrixXd::Zero(numVertices, numVertices); //映射v坐标方程系数矩阵
		Eigen::Matrix<double, Eigen::Dynamic, 1> b_x = Eigen::MatrixXd::Zero(numVertices, 1); //映射u坐标方程常数向量
		Eigen::Matrix<double, Eigen::Dynamic, 1> b_y = Eigen::MatrixXd::Zero(numVertices, 1); //映射v坐标方程常数向量
		Eigen::Matrix<double, Eigen::Dynamic, 1> res_x = Eigen::MatrixXd::Zero(numVertices, 1); //映射u坐标结果
		Eigen::Matrix<double, Eigen::Dynamic, 1> res_y = Eigen::MatrixXd::Zero(numVertices, 1); //映射v坐标结果

		//对每个点构建一条方程
		for(MeshVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> vertex_iterator(mesh);
			!vertex_iterator.end();
			++vertex_iterator)
		{
			//如果v是边界点，则构建方程 v = v.uv
			if (vertex_iterator.value()->boundary())
			{
				b_x(VMAP[vertex_iterator.value()->id()]) = vertex_iterator.value()->uv()[0];
				b_y(VMAP[vertex_iterator.value()->id()]) = vertex_iterator.value()->uv()[1];
				
				matrix_x(VMAP[vertex_iterator.value()->id()], VMAP[vertex_iterator.value()->id()]) = 1;
				matrix_y(VMAP[vertex_iterator.value()->id()], VMAP[vertex_iterator.value()->id()]) = 1;
				continue;
			}

			//如果是内部点i,且点i周围点为j，则构建方程 Σ(Vj.weight*Vj)-Vi*Σ(weight*Vj) = 0
			CToolVertex* vertex = vertex_iterator.value();
			double totalWeight = 0.0;

			//求出Σ(weight*Vj)
			for (VertexVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> vertex_iter(vertex);
				!vertex_iter.end();
				vertex_iter++)
			{
				totalWeight += GetVertexWeight(vertex_iter.value(), vertex,mesh);
			}

			//Vi系数-Σ(weight*Vj)
			matrix_x(VMAP[vertex_iterator.value()->id()], VMAP[vertex_iterator.value()->id()]) = -totalWeight;
			matrix_y(VMAP[vertex_iterator.value()->id()], VMAP[vertex_iterator.value()->id()]) = -totalWeight;
			//Vj系数Vj.weight
			for (VertexVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> vertex_iter(vertex);
				!vertex_iter.end();
				vertex_iter++)
			{
				matrix_x(VMAP[vertex->id()], VMAP[vertex_iter.value()->id()]) = GetVertexWeight(vertex_iter.value(),vertex, mesh);
				matrix_y(VMAP[vertex->id()], VMAP[vertex_iter.value()->id()]) = GetVertexWeight(vertex_iter.value(),vertex, mesh);
			}
		}

		//使用QR分解求方程组的解
		res_x = matrix_x.colPivHouseholderQr().solve(b_x);
		res_y = matrix_y.colPivHouseholderQr().solve(b_y);

		//赋值UV坐标
		for (MeshVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> vertex_iterator(mesh);
			!vertex_iterator.end();
			++vertex_iterator)
		{
			vertex_iterator.value()->uv()[0] = res_x[VMAP[vertex_iterator.value()->id()]];
			vertex_iterator.value()->uv()[1] = res_y[VMAP[vertex_iterator.value()->id()]];
		}
	}

	template<typename M,typename V>
	double HarmonyMapping<M,V>::GetVertexWeight(V* v1, V* v2,M* mesh)
	{
		double res = 0.0;
		CHalfEdge* half_edge = nullptr;

		//求出两点之间的半边
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

		//求出三角形的第三个点
		CVertex* v3 = half_edge->he_next()->vertex();
		res+= (v1->point() - v3->point()) * (v2->point() - v3->point())
			/ ((v1->point() - v3->point()) ^ (v2->point() - v3->point())).norm();

		//如果是边界点，则权重只有一边的三角形
		if (v1->boundary() && v2->boundary())
		{
			return res;
		}

		//指向另外的三角形
		half_edge = half_edge->he_sym();
		CVertex* v4 = half_edge->he_next()->vertex();
		res+= (v1->point() - v4->point()) * (v2->point() - v4->point())
			/ ((v1->point() - v4->point()) ^ (v2->point() - v4->point())).norm();
		return res;
	}

	template<typename M,typename V>
	CHalfEdge* HarmonyMapping<M,V>::FindNextBoundaryHE(CHalfEdge* he)
	{
		he = he->he_next();
		while (!he->edge()->boundary())
		{
			he = he->he_sym()->he_next();
		}
		return he;
	}

	template<typename M,typename V>
	void HarmonyMapping<M,V>::SetUV_CircleBorder(M* mesh)
	{
		int boundarySum = 0;//计算边界边数量
		CHalfEdge* initHE = NULL;

		//找到第一个位于边界的半边
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

		//计算半边总数
		CHalfEdge* half_iter = FindNextBoundaryHE(half_edge);
		for (;half_iter!=half_edge;half_iter = FindNextBoundaryHE(half_iter))
		{
			boundarySum++;
		}
		half_iter = FindNextBoundaryHE(half_iter);
		boundarySum++;
		
		//边界半边长度数组
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

		//计算边界点映射uv坐标
		for (int i = 0; i < boundarySum; i++)
		{
			half_iter->vertex()->uv() = CPoint2(
				cos(2 * PI * boundaryNorm[i] / boundaryNorm[boundarySum - 1]),
				sin(2 * PI * boundaryNorm[i] / boundaryNorm[boundarySum - 1])
			);
			half_iter = FindNextBoundaryHE(half_iter);
		}
	}

}

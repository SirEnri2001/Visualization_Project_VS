
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
		map<int,int>* pVMap = NULL;
		M* p_mesh;
		void SetTransfer(M* mesh);
	public:
		void Set_InsideUV(M* mesh);
		void SetUV_CircleBorder(M* mesh);
		double GetVertexWeight(V* v1, V* v2,M* mesh);
		void StartMapping();
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
			VMAP[vertex_iterator.value()->id()]=i++;
		}
	}


	template<typename M,typename V>
	void HarmonyMapping<M,V>::Set_InsideUV(M* mesh)
	{
		const int numVertices = mesh->numVertices()+1;
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_x = Eigen::MatrixXd::Zero(numVertices, numVertices);
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_y = Eigen::MatrixXd::Zero(numVertices, numVertices);
		Eigen::Matrix<double, Eigen::Dynamic, 1> b_x = Eigen::MatrixXd::Zero(numVertices, 1);
		Eigen::Matrix<double, Eigen::Dynamic, 1> b_y = Eigen::MatrixXd::Zero(numVertices, 1);
		Eigen::Matrix<double, Eigen::Dynamic, 1> res_x = Eigen::MatrixXd::Zero(numVertices, 1);
		Eigen::Matrix<double, Eigen::Dynamic, 1> res_y = Eigen::MatrixXd::Zero(numVertices, 1);
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

	template<typename M,typename V>
	double HarmonyMapping<M,V>::GetVertexWeight(V* v1, V* v2,M* mesh)
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

	template<typename M,typename V>
	void HarmonyMapping<M,V>::SetUV_CircleBorder(M* mesh)
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

}

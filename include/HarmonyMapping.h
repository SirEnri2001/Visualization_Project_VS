
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
		map<int,int>* pVMap = NULL; //�����±��VertexID���ɵ�ӳ��
		M* p_mesh;//����
		void SetTransfer(M* mesh);//��ʼ�������±��VertexID��ӳ��
		CHalfEdge* FindNextBoundaryHE(CHalfEdge* he);//�ҵ��߽�����ָ�����һ�����
	public:
		void Set_InsideUV(M* mesh);//�߽�ӳ���趨�õ�ǰ���£������ڲ���ӳ��
		void SetUV_CircleBorder(M* mesh);//�趨�߽�ӳ��ΪԲ�εļ�ӳ��
		double GetVertexWeight(V* v1, V* v2,M* mesh);//���������Ȩ��
		void StartMapping();//����ӳ����������
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
			VMAP[vertex_iterator.value()->id()]=i++; //���������ڵ����˳�򹹽�map
		}
	}


	template<typename M,typename V>
	void HarmonyMapping<M,V>::Set_InsideUV(M* mesh)
	{
		const int numVertices = mesh->numVertices()+1; //������ά��
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_x = Eigen::MatrixXd::Zero(numVertices, numVertices); //ӳ��u���귽��ϵ������
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_y = Eigen::MatrixXd::Zero(numVertices, numVertices); //ӳ��v���귽��ϵ������
		Eigen::Matrix<double, Eigen::Dynamic, 1> b_x = Eigen::MatrixXd::Zero(numVertices, 1); //ӳ��u���귽�̳�������
		Eigen::Matrix<double, Eigen::Dynamic, 1> b_y = Eigen::MatrixXd::Zero(numVertices, 1); //ӳ��v���귽�̳�������
		Eigen::Matrix<double, Eigen::Dynamic, 1> res_x = Eigen::MatrixXd::Zero(numVertices, 1); //ӳ��u������
		Eigen::Matrix<double, Eigen::Dynamic, 1> res_y = Eigen::MatrixXd::Zero(numVertices, 1); //ӳ��v������

		//��ÿ���㹹��һ������
		for(MeshVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> vertex_iterator(mesh);
			!vertex_iterator.end();
			++vertex_iterator)
		{
			//���v�Ǳ߽�㣬�򹹽����� v = v.uv
			if (vertex_iterator.value()->boundary())
			{
				b_x(VMAP[vertex_iterator.value()->id()]) = vertex_iterator.value()->uv()[0];
				b_y(VMAP[vertex_iterator.value()->id()]) = vertex_iterator.value()->uv()[1];
				
				matrix_x(VMAP[vertex_iterator.value()->id()], VMAP[vertex_iterator.value()->id()]) = 1;
				matrix_y(VMAP[vertex_iterator.value()->id()], VMAP[vertex_iterator.value()->id()]) = 1;
				continue;
			}

			//������ڲ���i,�ҵ�i��Χ��Ϊj���򹹽����� ��(Vj.weight*Vj)-Vi*��(weight*Vj) = 0
			CToolVertex* vertex = vertex_iterator.value();
			double totalWeight = 0.0;

			//�����(weight*Vj)
			for (VertexVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> vertex_iter(vertex);
				!vertex_iter.end();
				vertex_iter++)
			{
				totalWeight += GetVertexWeight(vertex_iter.value(), vertex,mesh);
			}

			//Viϵ��-��(weight*Vj)
			matrix_x(VMAP[vertex_iterator.value()->id()], VMAP[vertex_iterator.value()->id()]) = -totalWeight;
			matrix_y(VMAP[vertex_iterator.value()->id()], VMAP[vertex_iterator.value()->id()]) = -totalWeight;
			//Vjϵ��Vj.weight
			for (VertexVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> vertex_iter(vertex);
				!vertex_iter.end();
				vertex_iter++)
			{
				matrix_x(VMAP[vertex->id()], VMAP[vertex_iter.value()->id()]) = GetVertexWeight(vertex_iter.value(),vertex, mesh);
				matrix_y(VMAP[vertex->id()], VMAP[vertex_iter.value()->id()]) = GetVertexWeight(vertex_iter.value(),vertex, mesh);
			}
		}

		//ʹ��QR�ֽ��󷽳���Ľ�
		res_x = matrix_x.colPivHouseholderQr().solve(b_x);
		res_y = matrix_y.colPivHouseholderQr().solve(b_y);

		//��ֵUV����
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

		//�������֮��İ��
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

		//��������εĵ�������
		CVertex* v3 = half_edge->he_next()->vertex();
		res+= (v1->point() - v3->point()) * (v2->point() - v3->point())
			/ ((v1->point() - v3->point()) ^ (v2->point() - v3->point())).norm();

		//����Ǳ߽�㣬��Ȩ��ֻ��һ�ߵ�������
		if (v1->boundary() && v2->boundary())
		{
			return res;
		}

		//ָ�������������
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
		int boundarySum = 0;//����߽������
		CHalfEdge* initHE = NULL;

		//�ҵ���һ��λ�ڱ߽�İ��
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

		//����������
		CHalfEdge* half_iter = FindNextBoundaryHE(half_edge);
		for (;half_iter!=half_edge;half_iter = FindNextBoundaryHE(half_iter))
		{
			boundarySum++;
		}
		half_iter = FindNextBoundaryHE(half_iter);
		boundarySum++;
		
		//�߽��߳�������
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

		//����߽��ӳ��uv����
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

#ifndef _TOOL_H_
#define _TOOL_H_

#include<vector>
#include<math.h>

#include "ToolMesh.h"
#include "HarmonyMapping.h"
#include "ComplexIntegral.h"
#include "SphereMapping.h"

#ifndef PI
#define PI 3.1415926535
#endif




namespace MeshLib
{
	using namespace std;

	

	template<typename M,typename V>
	class CTool
	{
	public:
		CTool(M* pMesh);
		~CTool(){};

		void test();
		void _change_color();
		void GaussBonetTheory();
		double GetVertexWeight(V* v1,V * v2,M* mesh);
		double SumOfAngles(V* pV);
		
	
	protected:
		typename M* m_pMesh;
	};

	template<typename M,typename V>
	double CTool<M,V>::SumOfAngles(V* pV)
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

	
	template<typename M,typename V>
	double GetVertexWeight(V* v1, V* v2,M* mesh)
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
	CTool<M,V>::CTool(M* pMesh)
	{
		m_pMesh = pMesh;
	}

	template<typename M,typename V>
	void CTool<M,V>::test()
	{
		cout << "mesh vertex num: " << m_pMesh->numVertices() << endl;

	}

	template<typename M,typename V>
	void CTool<M,V>::_change_color()
	{
		for (M::MeshVertexIterator mv(m_pMesh); !mv.end(); mv++)
		{
			M::CVertex* pVertex = mv.value();
			pVertex->rgb()[0] = 1;
			pVertex->rgb()[1] = 1;
			pVertex->rgb()[2] = 0;
		}
	}

	template <typename M,typename V>
	void CTool<M,V>::GaussBonetTheory()
	{
		M* mesh = m_pMesh;
		double res = 0.0;
		for (MeshVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> v(mesh); !v.end(); v++)
		{
			if (!v.value()->boundary())
			{
				res += PI;
			}
			res += PI - SumOfAngles(v.value());
		}
		cout << "Gauss-Bonet Theory:" << endl;
		if (abs(res - 2 * PI * (mesh->numFaces() - mesh->numEdges() + mesh->numVertices())) > 0.001)
		{
			cout << "Not Match" << endl;
			cout << "Left Value " << res << endl;
			cout << "Right Value " << 2 * PI * (mesh->numFaces() - mesh->numEdges() + mesh->numVertices()) << endl;
		}
		else
		{
			cout << "Match!" << endl;
			cout << "Left Value " << res << endl;
			cout << "Right Value " << 2 * PI * (mesh->numFaces() - mesh->numEdges() + mesh->numVertices()) << endl;
		}
	}

}
#endif
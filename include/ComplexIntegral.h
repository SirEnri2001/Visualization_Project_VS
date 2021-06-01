#pragma once
#include "MeshLib/Geometry/Point2.H"
#include "MeshLib/Geometry/Point.h"
#include <complex>

namespace MeshLib
{
	using namespace Eigen;

	template<typename M,typename V>
	class ComplexIntegral
	{
	public:
		complex<double>  ComplexFunc(complex<double> z);
		complex<double> CompNum(CPoint cpoint);
		CPoint2 getPoint2(complex<double> c);
		CPoint2 getDeltaIntegral(V* v1, V* v2);
		void IntegralRound(V* vertex);
		void Integral(M* mesh);
	};
	


	template<typename M,typename V>
	complex<double> ComplexIntegral<M,V>::ComplexFunc(complex<double>  z)
	{
		complex<double> res;
		complex<double> z0(1.0,2.0);
		complex<double> z1(2.0,4.0);
		//res = (z-z1).conjugate();
		//res(0) /= pow( abs(z(0)),2);
		//res = (z-z0)*res(0);
		res = (z - z0) * (z - z1);
		return res * 0.01;
		//return z*0.01;
	}

	template<typename M,typename V>
	complex<double> ComplexIntegral<M,V>::CompNum(CPoint cpoint)
	{
		complex<double> c(cpoint[0],cpoint[1]);
		return c;
	}

	template<typename M,typename V>
	CPoint2 ComplexIntegral<M,V>::getPoint2(complex<double> c)
	{
		return CPoint2(c.real(),c.imag());
	}

	template<typename M,typename V>
	CPoint2 ComplexIntegral<M,V>::getDeltaIntegral(V* v1, V* v2)
	{
		return getPoint2((CompNum(v1->point()) - CompNum(v2->point()))
			* ComplexFunc((CompNum(v1->point() / 2) + CompNum(v2->point() / 2))));
	}

	template<typename M,typename V>
	void ComplexIntegral<M,V>::IntegralRound(V* vertex)
	{
		bool flag = false;
		for (VertexVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> vertex_iterator(vertex);
			!vertex_iterator.end();
			vertex_iterator++)
		{
			if (!vertex_iterator.value()->isInteg)
			{
				vertex_iterator.value()->uv() = getDeltaIntegral(vertex_iterator.value(), vertex) + vertex->uv();
				vertex_iterator.value()->isInteg = true;
				flag = true;
			}
		}
		if (!flag)
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

	template<typename M,typename V>
	void ComplexIntegral<M,V>::Integral(M* mesh)
	{
		CPoint2 zero(0, 0);
		CToolVertex* zeroPoint = NULL;
		for (MeshVertexIterator<CToolVertex, CToolEdge, CToolFace, CToolHalfEdge> vertex_iterator(mesh);
			!vertex_iterator.end();
			vertex_iterator++)
		{
			double distance = 1.0;
			if (distance > vertex_iterator.value()->point().norm())
			{
				distance = vertex_iterator.value()->point().norm();
				zeroPoint = vertex_iterator.value();
			}
		}
		zeroPoint->uv() = zero;
		zeroPoint->isInteg = true;
		IntegralRound(zeroPoint);
	}

}

#pragma once
#include <iostream>
#include <vector>
#include "HarmonyMapping.h"

namespace MeshLib
{
	template<typename M,typename V>
	class SphereMapping
	{
		HarmonyMapping<M,V> harmony_mapping; //平面调和映射实例
		M* p_mesh; //亏格为零、边界为零的网格
		std::vector<V> mesh_border;

		//TODO:构建高斯映射
		void GaussMapping();
		
		//TODO:边界的选取(模型的剪开)
		void CutBorder();
		
		//TODO:模型两部分的调和映射构建，以及球极投影
		void BetterMappingProceed();
		
		//TODO:球面映射调和优化
		void OptimizeMapping();
		
		void InitialMapping();
	public:
		void StartMapping();
		
		SphereMapping(M* pMesh);
		~SphereMapping()
		{
			delete pVMap;
		}
		
	};

	template <typename M, typename V>
	void SphereMapping<M, V>::StartMapping()
	{
		InitialMapping();
		OptimizeMapping();
	}


	template <typename M, typename V>
	SphereMapping<M, V>::SphereMapping(M* pMesh)
	{
		p_mesh = pMesh;
		harmony_mapping = HarmonyMapping<M,V>(p_mesh);
	}

	template <typename M, typename V>
	void SphereMapping<M, V>::InitialMapping()
	{
		GaussMapping();

		/*
		 * CutBorder();
		 * BetterMappingProceed();
		 */

		
	}

	template <typename M, typename V>
	void SphereMapping<M, V>::GaussMapping()
	{
		//TODO:构建高斯映射
		
	}


	template <typename M, typename V>
	void SphereMapping<M, V>::CutBorder()
	{
		//TODO:边界的选取(模型的剪开)
		mesh_border = new std::vector<V>();
	}

	template <typename M, typename V>
	void SphereMapping<M, V>::BetterMappingProceed()
	{
		//TODO:模型两部分的调和映射构建，以及球极投影
		
	}

	template <typename M, typename V>
	void SphereMapping<M, V>::OptimizeMapping()
	{
		//TODO:球面映射调和优化
		
	}

	
}

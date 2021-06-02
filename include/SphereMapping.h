#pragma once
#include <iostream>
#include <vector>
#include "HarmonyMapping.h"

namespace MeshLib
{
	template<typename M,typename V>
	class SphereMapping
	{
		HarmonyMapping<M,V> harmony_mapping; //ƽ�����ӳ��ʵ��
		M* p_mesh; //����Ϊ�㡢�߽�Ϊ�������
		std::vector<V> mesh_border;

		//TODO:������˹ӳ��
		void GaussMapping();
		
		//TODO:�߽��ѡȡ(ģ�͵ļ���)
		void CutBorder();
		
		//TODO:ģ�������ֵĵ���ӳ�乹�����Լ���ͶӰ
		void BetterMappingProceed();
		
		//TODO:����ӳ������Ż�
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
		//TODO:������˹ӳ��
		
	}


	template <typename M, typename V>
	void SphereMapping<M, V>::CutBorder()
	{
		//TODO:�߽��ѡȡ(ģ�͵ļ���)
		mesh_border = new std::vector<V>();
	}

	template <typename M, typename V>
	void SphereMapping<M, V>::BetterMappingProceed()
	{
		//TODO:ģ�������ֵĵ���ӳ�乹�����Լ���ͶӰ
		
	}

	template <typename M, typename V>
	void SphereMapping<M, V>::OptimizeMapping()
	{
		//TODO:����ӳ������Ż�
		
	}

	
}

/*
Copyright (c) 2026, Michael Kazhdan and Hongyi Liu
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef TEXTURE_ADVECTION_INCLUDED
#define TEXTURE_ADVECTION_INCLUDED

#include "Misha/Exceptions.h"
#include "Misha/FEM.h"
#include "Misha/MultiThreading.h"
#include "Misha/Texels.h"
#include "Include/SimplicialMeshProcessing.h"
#include "TexturedMesh.h"
#include "TextureAndConfidence.h"

namespace MishaK
{
	namespace SimplicialMesh
	{
		template< typename Getter , unsigned int K , typename T >
		concept HasTextureAdvectionGetter = requires( const Getter getter , size_t idx , Point< double , K > p ) { { getter(idx,p) } -> std::same_as< T >; };
		template< typename Setter , typename T >
		concept HasTextureAdvectionSetter = requires( const Setter setter , unsigned int i , size_t idx , T t ) { { setter(i,idx,t) } -> std::same_as< void >; };

		template< unsigned int K >
		struct TextureAdvection
		{
			struct Parameters
			{
				double minStepSize , time;
				explicit Parameters( double minStepSize , double time=1. ) : minStepSize(minStepSize) , time(time) {}
			};

			bool verbose;

			TextureAdvection( EmbeddedPhongTexturedMesh< K , K > mesh , const unsigned int res[K] , bool verbose=false );
			~TextureAdvection( void );

			template< typename T >
			RegularGrid< K , T > operator()( const std::vector< Point< double , K+1 > > &flow , const RegularGrid< K , T > &input , Parameters params ) const;

			template< typename T >
			TextureAndConfidence < K, T > operator()(const std::vector< Point< double, K + 1 > >& flow, const TextureAndConfidence < K, T >& input, double time = 1., double minStepSize = 0.) const;

			template< typename ... Ts >
			void process( const std::vector< Point< double , K+1 > > &flow , const RegularGrid< K , Ts >& ... in , RegularGrid< K , Ts >& ... out , Parameters params ) const;
		protected:
			EmbeddedPhongTexturedMesh< K , K > _mesh;
			FEM::RiemannianMesh< double , unsigned int > _rMesh;
			Pointer( FEM::CoordinateXForm< double > ) _xForms;
			RegularGrid< K , typename Texels< true >::TexelInfo< K > > _texelInfo;

			struct _TangentVectorField : public FEM::TangentVectorField< double >
			{
				_TangentVectorField( const EmbeddedPhongMesh< K > &mesh , const std::vector< Point< double , K+1 > > & flowField ) : _vf( mesh , flowField ) {}
				FEM::TangentVector< double > operator()( const FEM::SamplePoint< double > &s ) const { return _vf[s.tIdx](s.p); }
			protected:
				PhongRodriguesIntrinsicVectorField< K > _vf;
			};

			template< typename T , HasTextureAdvectionGetter< K , T > Getter , HasTextureAdvectionSetter< T > Setter >
			void _process( const std::vector< Point< double , K+1 > > &flow , Getter && getF , Setter && setF , Parameters params ) const;

			template< bool First , typename T , typename ... Ts >
			static void _SetRes( const RegularGrid< K , T > &in , const RegularGrid< K , Ts >& ... ins , RegularGrid< K , T > &out , RegularGrid< K , Ts >& ...outs , unsigned int res[K] );

			template< unsigned int I , typename ... _T , typename T , typename ...Ts >
			static void _ToTuple( Point< double , K > p , std::tuple< _T... > &values , const RegularGrid< K , T > &t , const RegularGrid< K , Ts >& ... ts );

			template< unsigned int I , typename ... _T , typename T , typename ...Ts >
			static void _FromTuple( size_t i , const std::tuple< _T... > &values , RegularGrid< K , T > &t , RegularGrid< K , Ts >& ... ts );

			template< typename T , typename ... Ts >
			static void _Copy( const RegularGrid< K , T > &in , const RegularGrid< K , Ts >& ... ins , RegularGrid< K , T > &out , RegularGrid< K , Ts >& ...outs );
		};

#include "TextureAdvection.inl"
	}
}
#endif // TEXTURE_ADVECTION_INCLUDED
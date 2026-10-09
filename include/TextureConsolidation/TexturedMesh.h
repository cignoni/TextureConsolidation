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

// Texture-atlas support layered on top of PhongRodriguesVF's mesh classes:
//   EmbeddedTexturedMesh      - an embedded mesh with per-corner texture coordinates
//   EmbeddedPhongTexturedMesh - a Phong mesh with per-corner texture coordinates
//   TextureScalarField        - a texture-image signal, queryable (with differential) per simplex
#ifndef TEXTURED_MESH_INCLUDED
#define TEXTURED_MESH_INCLUDED

#include <vector>
#include "Misha/Geometry.h"
#include "Misha/RegularGrid.h"
#include "Include/SimplicialMeshProcessing.h"

namespace MishaK
{
	namespace SimplicialMesh
	{
		template< unsigned int K , unsigned int Dim , unsigned int TDim >
		struct EmbeddedTexturedMesh : public virtual EmbeddedMesh< K , Dim >
		{
			EmbeddedTexturedMesh( std::vector< Point< double , Dim > > &vertices , std::vector< Point< double , TDim > > &textureCoordinates , std::vector< SimplexIndex< K > > &simplices )
				: EmbeddedMesh< K , Dim >( vertices , simplices ) , _textureCoordinates( textureCoordinates )
			{}

			Point< double , TDim > &textureCoordinate( size_t t ){ return _textureCoordinates[t]; }
			std::vector< Point< double , TDim > > &textureCoordinates( void ){ return _textureCoordinates; }
			const Point< double , TDim > &textureCoordinate( size_t t ) const { return _textureCoordinates[t]; }
			const std::vector< Point< double , TDim > > &textureCoordinates( void ) const { return _textureCoordinates; }

		protected:
			using EmbeddedMesh< K , Dim >::_vertices;
			using EmbeddedMesh< K , Dim >::_simplices;
			std::vector< Point< double , TDim > > &_textureCoordinates;
		};

		template< unsigned int K , unsigned int TDim >
		struct EmbeddedPhongTexturedMesh : public EmbeddedPhongMesh< K >
		{
			static const unsigned int Dim = K+1;
			EmbeddedPhongTexturedMesh( std::vector< Point< double , Dim > > &vertices , std::vector< Point< double , Dim > > &normals , std::vector< Point< double , TDim > > &textureCoordinates , std::vector< SimplexIndex< K > > &simplices )
				: EmbeddedPhongMesh< K >( vertices , normals , simplices ) , _textureCoordinates( textureCoordinates )
			{}

			operator       EmbeddedTexturedMesh< K , K+1 , TDim >()       { return EmbeddedTexturedMesh< K , K+1 , TDim >( _vertices , _textureCoordinates , _simplices ); }
			operator const EmbeddedTexturedMesh< K , K+1 , TDim >() const { return EmbeddedTexturedMesh< K , K+1 , TDim >( _vertices , _textureCoordinates , _simplices ); }

			Point< double , TDim > &textureCoordinate( size_t t ){ return _textureCoordinates[t]; }
			std::vector< Point< double , TDim > > &textureCoordinates( void ){ return _textureCoordinates; }
			const Point< double , TDim > &textureCoordinate( size_t t ) const { return _textureCoordinates[t]; }
			const std::vector< Point< double , TDim > > &textureCoordinates( void ) const { return _textureCoordinates; }

		protected:
			using EmbeddedPhongMesh< K >::_vertices;
			using EmbeddedPhongMesh< K >::_simplices;
			std::vector< Point< double , TDim > > &_textureCoordinates;
		};

		// The function evaluating a texture-based signal
		template< unsigned int K , typename T , unsigned int TDim >
		struct TextureScalarField
		{
			TextureScalarField( const std::vector< Point< double , TDim > > &textureCoordinates , const RegularGrid< TDim , T > & signal ) : _textureCoordinates( textureCoordinates ) , _signal(signal) {}

			// Given a simplex, returns a function returning the evaluation at a point in the simplex
			auto operator[]( size_t sIdx ) const;

			struct Differential
			{
				Differential( const TextureScalarField & sf ) : _tsf(sf) {}
				Differential( TextureScalarField && sf ) : _tsf(sf) {}

				// Given a simplex, returns a function returning the differential at a point in the simplex
				auto operator[]( size_t sIdx ) const;
			protected:
				TextureScalarField _tsf;
			};

		protected:
			const std::vector< Point< double , TDim > > & _textureCoordinates;
			const RegularGrid< TDim , T > &_signal;
		};

#include "TexturedMesh.inl"
	}
}
#endif // TEXTURED_MESH_INCLUDED

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

#ifndef TEXTURE_COMPARISON_INCLUDED
#define TEXTURE_COMPARISON_INCLUDED

#include <vector>
#include "Misha/RegularGrid.h"
#include "Misha/Geometry.h"
#include "Include/SimplicialMeshProcessing.h"
#include "TexturedMesh.h"
#include "TextureAndConfidence.h"

namespace MishaK
{
	namespace SimplicialMesh
	{
		// Integrals of texture/confidence images over a textured mesh, used to order the registration and to report its error
		template< unsigned int Quadrature , unsigned int K , unsigned int Dim >
		struct TextureComparer
		{
			TextureComparer( EmbeddedTexturedMesh< K , Dim , K > mesh ) : _mesh(mesh) {}

			// The coverage of a confidence, < w > = \int w
			double coverage( const RegularGrid< K , double > &confidence ) const;

			// The overlap of two confidences, < w1 * w2 > = \int w1 * w2
			double overlap( const RegularGrid< K , double > &confidence1 , const RegularGrid< K , double > &confidence2 ) const;

			// The RMS color difference over all pairs of textures, weighted by the product of their confidences
			template< unsigned int Channels >
			double error( const std::vector< TextureAndConfidence< K , Point< double , Channels > > > &textureAndConfidences ) const;

		protected:
			EmbeddedTexturedMesh< K , Dim , K > _mesh;
		};

#include "TextureComparison.inl"
	}
}
#endif // TEXTURE_COMPARISON_INCLUDED

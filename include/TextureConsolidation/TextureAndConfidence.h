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

#ifndef TEXTURE_AND_CONFIDENCE_INCLUDED
#define TEXTURE_AND_CONFIDENCE_INCLUDED

#include <functional>
#include <set>
#include <vector>
#include <Misha/RegularGrid.h>
#include <Misha/Geometry.h>
#include <Misha/MultiThreading.h>
#include <Misha/Miscellany.h>

namespace MishaK
{
	namespace SimplicialMesh
	{
		template< unsigned int K , typename T >
		struct TextureAndConfidence
		{
			RegularGrid< K , T > texture;
			RegularGrid< K , double > confidence;

			TextureAndConfidence( void ){}
			TextureAndConfidence( const unsigned int res[K] ) : texture( res ) , confidence( res ){}
		};

		template<  unsigned int K, unsigned int Channels >
		TextureAndConfidence< K, Point< double, Channels > > StitchByWeightedConfidence(
			const TextureAndConfidence< K, Point< double, Channels > >& TC1,
			const TextureAndConfidence< K, Point< double, Channels > >& TC2
		)
		{
			TextureAndConfidence< K, Point< double, Channels > > outTC = TC1;

			// stitch the source and target texture by taking confidence weighted average in all color channels
			RegularGrid< K, Point< double, Channels > >& stitched = outTC.texture;
			for (unsigned int i = 0; i < stitched.resolution(); i++)
			{
				for (int c = 0; c < Channels; c++)
					stitched[i][c] = (TC1.texture[i][c] * TC1.confidence[i] + TC2.texture[i][c] * TC2.confidence[i]) / (TC1.confidence[i] + TC2.confidence[i] + 1e-10);
			}

			// stitch the confidence texture 
			RegularGrid< K, double >& confidence = outTC.confidence;
			for (unsigned int i = 0; i < stitched.resolution(); i++)
				confidence[i] = std::max(TC1.confidence[i], TC2.confidence[i]);

			return outTC;
		}

		template<  unsigned int K, unsigned int Channels >
		TextureAndConfidence< K, Point< double, Channels > > StitchByWeightedConfidence(const std::vector<TextureAndConfidence< K, Point< double, Channels > > >& TCs, const std::set<unsigned int>& usedTextures)
		{
			TextureAndConfidence< K, Point< double, Channels > > outTC;
			if (TCs.empty() || usedTextures.empty())
				return outTC;

			outTC.texture.resize(TCs[0].texture.res());
			outTC.confidence.resize(TCs[0].confidence.res());

			Point < double, Channels > zero_point;
			for (unsigned int c = 0; c < Channels; c++)
				zero_point[c] = 0.0;

			// For each texel, compute weighted average of used textures
			for (unsigned int i = 0; i < outTC.texture.resolution(); i++)
			{
				double sum_conf = 0;
				Point<double, Channels> weighted_sum; // is this zero by default
				for (auto idx : usedTextures)
				{
					double conf = TCs[idx].confidence[i];
					sum_conf += conf;
					for (int c = 0; c < Channels; c++)
						weighted_sum[c] += TCs[idx].texture[i][c] * conf;
				}
				if (sum_conf > 1e-10)
				{
					outTC.texture[i] = weighted_sum / sum_conf;
					outTC.confidence[i] = std::min< double >( sum_conf , 1. ); // the sum of the confidences of the used textures, capped at one
				}
				else
				{
					outTC.texture[i] = zero_point;
					outTC.confidence[i] = 0.0;
				}
			}
			return outTC;
		}

		template<  unsigned int K, unsigned int Channels >
		TextureAndConfidence< K, Point< double, Channels > > StitchByWeightedConfidence(const std::vector<TextureAndConfidence< K, Point< double, Channels > > >& TCs)
		{
			std::set <unsigned int> usedTextures;
			for (unsigned int i = 0; i < TCs.size(); i++)
				usedTextures.insert(i);
			return StitchByWeightedConfidence(TCs, usedTextures);
		}

	}
}
#endif // TEXTURE_AND_CONFIDENCE_INCLUDED
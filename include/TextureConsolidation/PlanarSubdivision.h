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

#ifndef PLANAR_SUBDIVISION_INCLUDED
#define PLANAR_SUBDIVISION_INCLUDED

#include <vector>
#include <map>
#include <utility>
#include <array>
#include "Misha/Geometry.h"

namespace MishaK
{
	// One-to-four planar subdivision of a triangle mesh, inserting a vertex at the mid-point of every edge.
	// Unlike Loop-type schemes it does not need vertex one-rings, so meshes with boundaries (and any manifold
	// or non-manifold edge connectivity) are supported.
	// Triangle t is split into triangles 4t ... 4t+3 as:
	//		( e0 , e1 , e2 ) , ( v0 , e2 , e1 ) , ( e2 , v1 , e0 ) , ( e1 , e0 , v2 )
	// with v_k the corners of t and e_k the mid-point of the edge opposite v_k.
	struct PlanarSubdivider
	{
		PlanarSubdivider( const std::vector< SimplexIndex< 2 > > &triangles , size_t vNum ) : _triangles( triangles ) , _vNum( vNum )
		{
			std::map< std::pair< size_t , size_t > , size_t > edgeIndex;
			_edgeVertex.resize( triangles.size() );
			for( size_t t=0 ; t<triangles.size() ; t++ ) for( unsigned int k=0 ; k<3 ; k++ )
			{
				size_t v1 = triangles[t][(k+1)%3] , v2 = triangles[t][(k+2)%3];
				std::pair< size_t , size_t > key( std::min( v1 , v2 ) , std::max( v1 , v2 ) );
				auto iter = edgeIndex.find( key );
				if( iter==edgeIndex.end() )
				{
					size_t idx = _vNum + _edges.size();
					edgeIndex[key] = idx;
					_edges.push_back( key );
					_edgeVertex[t][k] = idx;
				}
				else _edgeVertex[t][k] = iter->second;
			}
		}

		// The subdivided triangles
		std::vector< SimplexIndex< 2 > > operator()( void ) const
		{
			std::vector< SimplexIndex< 2 > > triangles( _triangles.size()*4 );
			for( size_t t=0 ; t<_triangles.size() ; t++ )
			{
				size_t v[] = { (size_t)_triangles[t][0] , (size_t)_triangles[t][1] , (size_t)_triangles[t][2] } , e[] = { _edgeVertex[t][0] , _edgeVertex[t][1] , _edgeVertex[t][2] };
				_Set( triangles[4*t+0] , e[0] , e[1] , e[2] );
				_Set( triangles[4*t+1] , v[0] , e[2] , e[1] );
				_Set( triangles[4*t+2] , e[2] , v[1] , e[0] );
				_Set( triangles[4*t+3] , e[1] , e[0] , v[2] );
			}
			return triangles;
		}

		// Per-vertex data, with new vertices set to the average of the edge end-points
		template< typename T >
		std::vector< T > vertexData( const std::vector< T > &data ) const
		{
			std::vector< T > _data( _vNum + _edges.size() );
			for( size_t i=0 ; i<_vNum ; i++ ) _data[i] = data[i];
			for( size_t i=0 ; i<_edges.size() ; i++ ) _data[_vNum+i] = ( data[ _edges[i].first ] + data[ _edges[i].second ] ) / 2.;
			return _data;
		}

		// Per-corner data (e.g. texture coordinates), with new corners set to the average along the triangle's edge
		template< typename T >
		std::vector< T > cornerData( const std::vector< T > &data ) const
		{
			std::vector< T > _data( _triangles.size()*12 );
			for( size_t t=0 ; t<_triangles.size() ; t++ )
			{
				T v[3] , e[3];
				for( unsigned int k=0 ; k<3 ; k++ ) v[k] = data[3*t+k] , e[k] = ( data[3*t+(k+1)%3] + data[3*t+(k+2)%3] ) / 2.;
				T sub[] = { e[0] , e[1] , e[2] , v[0] , e[2] , e[1] , e[2] , v[1] , e[0] , e[1] , e[0] , v[2] };
				for( unsigned int i=0 ; i<12 ; i++ ) _data[12*t+i] = sub[i];
			}
			return _data;
		}

	protected:
		const std::vector< SimplexIndex< 2 > > &_triangles;
		size_t _vNum;
		std::vector< std::pair< size_t , size_t > > _edges;
		std::vector< std::array< size_t , 3 > > _edgeVertex;

		static void _Set( SimplexIndex< 2 > &t , size_t a , size_t b , size_t c ){ t[0] = (unsigned int)a , t[1] = (unsigned int)b , t[2] = (unsigned int)c; }
	};
}
#endif // PLANAR_SUBDIVISION_INCLUDED

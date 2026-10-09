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

template< unsigned int Quadrature , unsigned int K , unsigned int Dim >
double TextureComparer< Quadrature , K , Dim >::coverage( const RegularGrid< K , double > &confidence ) const
{
	TextureScalarField< K , double , K > field( _mesh.textureCoordinates() , confidence );
	return _mesh.template integral< Quadrature , double >( SimplexProcessing::ArrayWrapper( [&]( size_t sIdx ){ return field[sIdx]; } ) );
}

template< unsigned int Quadrature , unsigned int K , unsigned int Dim >
double TextureComparer< Quadrature , K , Dim >::overlap( const RegularGrid< K , double > &confidence1 , const RegularGrid< K , double > &confidence2 ) const
{
	TextureScalarField< K , double , K > field1( _mesh.textureCoordinates() , confidence1 ) , field2( _mesh.textureCoordinates() , confidence2 );
	auto Weight = [&]( size_t sIdx )
		{
			return [w1=field1[sIdx],w2=field2[sIdx]]( SimplexProcessing::Position< K > p ){ return w1(p) * w2(p); };
		};
	return _mesh.template integral< Quadrature , double >( SimplexProcessing::ArrayWrapper( Weight ) );
}

template< unsigned int Quadrature , unsigned int K , unsigned int Dim >
template< unsigned int Channels >
double TextureComparer< Quadrature , K , Dim >::error( const std::vector< TextureAndConfidence< K , Point< double , Channels > > > &textureAndConfidences ) const
{
	double l2 = 0 , measure = 0;
	for( unsigned int i=0 ; i<textureAndConfidences.size() ; i++ ) for( unsigned int j=0 ; j<i ; j++ )
	{
		TextureScalarField< K , Point< double , Channels > , K > t1( _mesh.textureCoordinates() , textureAndConfidences[i].texture ) , t2( _mesh.textureCoordinates() , textureAndConfidences[j].texture );
		TextureScalarField< K , double , K > w1( _mesh.textureCoordinates() , textureAndConfidences[i].confidence ) , w2( _mesh.textureCoordinates() , textureAndConfidences[j].confidence );
		auto SquaredDifference = [&]( size_t sIdx )
			{
				return [_t1=t1[sIdx],_t2=t2[sIdx],_w1=w1[sIdx],_w2=w2[sIdx]]( SimplexProcessing::Position< K > p ){ return Point< double , Channels >::SquareNorm( _t1(p)-_t2(p) ) * _w1(p) * _w2(p); };
			};
		l2 += _mesh.template integral< Quadrature , double >( SimplexProcessing::ArrayWrapper( SquaredDifference ) );
		measure += overlap( textureAndConfidences[i].confidence , textureAndConfidences[j].confidence );
	}
	return sqrt( l2 / measure );
}

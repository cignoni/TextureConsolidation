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

template< unsigned int Quadrature , typename Solver , unsigned int K , unsigned int Dim >
TextureNormalizer< Quadrature , Solver , K , Dim >::TextureNormalizer
(
	EmbeddedTexturedMesh< K , Dim , K > mesh ,
	unsigned int width ,
	unsigned int height ,
	bool verbose
)
	: _textureToVertex( mesh , verbose ) , _mesh(mesh) , verbose( verbose )
{
	resize( width , height );
}

template< unsigned int Quadrature , typename Solver , unsigned int K , unsigned int Dim >
void TextureNormalizer< Quadrature , Solver , K , Dim >::resize
(
	unsigned int width ,
	unsigned int height
)
{
	Miscellany::PerformanceMeter pMeter( '.' );
	unsigned int res[] = { width , height };
	auto TextureSimplexFunctor = [&]( size_t sIdx )
		{
			Simplex< double , K , K > simplex;
			for( unsigned int k=0 ; k<=K ; k++ ) simplex[k] = _mesh.textureCoordinate( sIdx*(K+1) + k );
			return simplex;
		};

	_activeTexels = Texels< true >::template GetSupportedTexelInfo< Dim , false >
		(
			_mesh.simplexNum() ,
			[&]( size_t v ){ return _mesh.vertex(v); } ,
			[&]( size_t s ){ return _mesh.simplex(s); } ,
			TextureSimplexFunctor ,
			res ,
			0 ,
			verbose
		);
	if( verbose ) std::cout << pMeter( "Got active texels" ) << std::endl;
}


template< unsigned int Quadrature , typename Solver , unsigned int K , unsigned int Dim >
void TextureNormalizer< Quadrature , Solver , K , Dim >::update
(
	double diffusionTime ,
	const RegularGrid< K , double > *confidence
)
{
	_textureToVertex.update( diffusionTime , confidence );
}

template< unsigned int Quadrature , typename Solver , unsigned int K , unsigned int Dim >
void TextureNormalizer< Quadrature , Solver , K , Dim >::update
(
	double diffusionTime
)
{
	return update( diffusionTime , nullptr );
}

template< unsigned int Quadrature , typename Solver , unsigned int K , unsigned int Dim >
template< unsigned int Channels >
void TextureNormalizer< Quadrature , Solver , K , Dim >::operator()
(
	RegularGrid< K , Point< double , Channels > > &texture ,
	const RegularGrid< K , double > *confidence
)
const
{
	using ColorField = TextureScalarField< K , Point< double , Channels > , K >;
	using SmoothBrightnessField = ScalarField< K , Point< double , 1 > >;
	using SmoothColorField = ScalarField< K , Point< double , Channels > >;
	using ConfidenceField = TextureScalarField< K , double , K >;

	auto BC = [&]( Point< double , Dim > bc )
		{
			Point< double , K > _bc;
			for( unsigned int k=0 ; k<K ; k++ ) _bc[k] = bc[k+1];
			return _bc;
		};

	struct BilinearForm
	{
		BilinearForm( double s=1. ) : scale(s) {}
		double operator()( double v1 , double v2 ) const { return scale * v1 * v2; }
		double operator()( double v ) const { return scale * v; }

		double scale;
	};

	if( confidence && texture.size()!=confidence->size() ) MK_ERROR_OUT( "Confidence size does not match texture size: " , confidence->size() , " != " , texture.size() );

	if( texture.res(0)!=_activeTexels.res(0) || texture.res(1)!=_activeTexels.res(1) )
		MK_ERROR_OUT( "Texture and active texel resolutions don't match: " , texture.res(0) , " x " , texture.res(1) , " != " , _activeTexels.res(0) , " x " , _activeTexels.res(1) );

	Miscellany::PerformanceMeter pMeter( '.' );

	Point< double , Channels > avg;
	for( unsigned int c=0 ; c<Channels ; c++ ) avg[c] = 0.5;

	// Adjust the translational offset
	{
		std::vector< Point< double , Channels > > smoothedSignal = _textureToVertex( texture , confidence );

		SmoothColorField colorField( _mesh , smoothedSignal );
		ThreadPool::ParallelFor
		(
			0 , _activeTexels.size() ,
			[&]( size_t i )
			{
				if( _activeTexels[i].sIdx!=-1 && ( !confidence || (*confidence)[i] ) )
				{
					texture[i] += avg - colorField[ _activeTexels[i].sIdx ]( BC( _activeTexels[i].bc ) );
				}
			}
		);
	}
	if( verbose ) std::cout << pMeter( "Shifted colors" ) << std::endl;

	// Adjust the multiplicative offset
	{
		RegularGrid< 2 , Point< double , 1 > > brightness( texture.res() );

		ThreadPool::ParallelFor
		(
			0 , _activeTexels.size() ,
			[&]( size_t i )
			{
				if( _activeTexels[i].sIdx!=-1 && ( !confidence || (*confidence)[i] ) ) brightness[i][0] = sqrt( Point< double , Channels >::SquareNorm( texture[i] - avg ) / Channels );
			}
		);
		if( verbose ) std::cout << pMeter( "Got brightness" ) << std::endl;

		std::vector< Point< double , 1 > > smoothedBrightness = _textureToVertex( brightness , confidence );

		SmoothBrightnessField scalarField( _mesh , smoothedBrightness );
		ThreadPool::ParallelFor
		(
			0 , _activeTexels.size() ,
			[&]( size_t i )
			{
				if( _activeTexels[i].sIdx!=-1 && ( !confidence || (*confidence)[i] ) )
					texture[i] = avg + ( texture[i] - avg ) * 0.125 / scalarField[ _activeTexels[i].sIdx ]( BC( _activeTexels[i].bc ) )[0];
			}
		);
	}
	if( verbose ) std::cout << pMeter( "Shifted brightness" ) << std::endl;
}

template< unsigned int Quadrature , typename Solver , unsigned int K , unsigned int Dim >
template< unsigned int Channels >
void TextureNormalizer< Quadrature , Solver , K , Dim >::operator()
(
	RegularGrid< K , Point< double , Channels > > &texture
)
const
{
	return operator()( texture , nullptr );
}


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
TextureToVertex< Quadrature , Solver , K , Dim >::TextureToVertex
(
	EmbeddedTexturedMesh< K , Dim , K > mesh ,
	bool verbose
)
	: _mesh( mesh ) , verbose( verbose )
{
	Miscellany::PerformanceMeter pMeter( '.' );
	_M = _mesh.template massMatrix< Quadrature >();
	_S = _mesh.template stiffnessMatrix< Quadrature >();
	if( verbose ) std::cout << pMeter( "Got system" ) << std::endl;

	_solver.analyzePattern( _M + _S );
	if( verbose ) std::cout << pMeter( "Analyzed" ) << std::endl;
}

template< unsigned int Quadrature , typename Solver , unsigned int K , unsigned int Dim >
void TextureToVertex< Quadrature , Solver , K , Dim >::update
(
	double diffusionTime ,
	const RegularGrid< K , double > *confidence
)
{
	Miscellany::PerformanceMeter pMeter( '.' );
	if( confidence )
	{
		TextureScalarField< K , double , K > confidenceField( _mesh.textureCoordinates() , *confidence );
		_M = _mesh.template massMatrix< Quadrature >( confidenceField );
	}
	else _M = _mesh.template massMatrix< Quadrature >();
	if( verbose ) std::cout << pMeter( "Got mass" ) << std::endl;

	_solver.factorize( _M + _S*diffusionTime );
	if( _solver.info()!=Eigen::Success ) MK_ERROR_OUT( "Failed to factorize system" );
	if( verbose ) std::cout << pMeter( "Factorized" ) << std::endl;
}

template< unsigned int Quadrature , typename Solver , unsigned int K , unsigned int Dim >
void TextureToVertex< Quadrature , Solver , K , Dim >::update
(
	double diffusionTime
)
{
	update( diffusionTime , nullptr );
}

template< unsigned int Quadrature , typename Solver , unsigned int K , unsigned int Dim >
template< unsigned int Channels >
std::vector< Point< double , Channels > > TextureToVertex< Quadrature , Solver , K , Dim >::operator()
(
	const RegularGrid< K , Point< double , Channels > > &texture ,
	const RegularGrid< K , double > *confidence
)
const
{
	std::vector< Point< double , Channels > > signal( _mesh.vertexNum() );

	TextureScalarField< K , Point< double , Channels > , K > textureField( _mesh.textureCoordinates() , texture );

	Miscellany::PerformanceMeter pMeter( '.' );
	if( confidence )
	{
		TextureScalarField< K , double , K > confidenceField( _mesh.textureCoordinates() , *confidence );
		for( unsigned int c=0 ; c<Channels ; c++ )
		{
			auto ScalarField = [&]( size_t sIdx ){ return [c,tex=textureField[sIdx],conf=confidenceField[sIdx]]( SimplexProcessing::Position< K > p ){ return tex(p)[c] * conf(p); }; };
			Eigen::VectorXd x = _solver.solve( _mesh.template massVector< Quadrature >( SimplexProcessing::ArrayWrapper( ScalarField ) ) );
			for( unsigned int i=0 ; i<_mesh.vertexNum() ; i++ ) signal[i][c] = x[i];
		}
	}
	else
	{
		for( unsigned int c=0 ; c<Channels ; c++ )
		{
			auto ScalarField = [&]( size_t sIdx ){ return [c,tex=textureField[sIdx]]( SimplexProcessing::Position< K > p ){ return tex(p)[c]; }; };
			Eigen::VectorXd x = _solver.solve( _mesh.template massVector< Quadrature >( SimplexProcessing::ArrayWrapper( ScalarField ) ) );
			for( unsigned int i=0 ; i<_mesh.vertexNum() ; i++ ) signal[i][c] = x[i];
		}
	}
	if( verbose ) std::cout << pMeter( "Solved" ) << std::endl;

	return signal;
}

template< unsigned int Quadrature , typename Solver , unsigned int K , unsigned int Dim >
template< unsigned int Channels >
std::vector< Point< double , Channels > > TextureToVertex< Quadrature , Solver , K , Dim >::operator()
(
	const RegularGrid< K , Point< double , Channels > > &texture
)
const
{
	return operator()( texture , nullptr );
}


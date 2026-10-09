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

////////////////////////
// TextureScalarField //
////////////////////////
template< unsigned int K , typename T , unsigned int TDim >
auto TextureScalarField< K , T , TDim >::operator[]( size_t sIdx ) const
{
	Simplex< double , K , TDim > s;
	for( unsigned int k=0 ; k<=K ; k++ )
	{
		s[k] = _textureCoordinates[ sIdx*(K+1) + k ];
		for( unsigned int d=0 ; d<TDim ; d++ ) s[k][d] *= _signal.res(d)-1;
	}

	return [s,&signal=_signal]( SimplexProcessing::Position< K > p )
		{
			Point< double , TDim > tCoordinate = s[0];
			for( unsigned int k=0 ; k<K ; k++ ) tCoordinate += ( s[k+1]-s[0] ) * p[k];

			// Adjust for the fact that RegularGrid places samples at corners not centers.
			return signal( tCoordinate );
		};
}

//////////////////////////////////////
// TextureScalarField::Differential //
//////////////////////////////////////
template< unsigned int K , typename T , unsigned int TDim >
auto TextureScalarField< K , T , TDim >::Differential::operator[]( size_t sIdx ) const
{
	Simplex< double , K , TDim > s;
	for( unsigned int k=0 ; k<=K ; k++ )
	{
		s[k] = _tsf._textureCoordinates[ sIdx*(K+1) + k ];
		for( unsigned int d=0 ; d<TDim ; d++ ) s[k][d] *= _tsf._signal.res(d)-1;
	}

	Matrix< double , TDim , K > t2s;
	{
		Matrix< double , K , TDim > s2t;
		for( unsigned int k=0 ; k<K ; k++ ) for( unsigned int d=0 ; d<TDim ; d++ ) s2t( k , d ) = s[k+1][d] - s[0][d];
		t2s = s2t.transpose();
	}

	return [s,t2s,&signal=_tsf._signal]( SimplexProcessing::Position< K > p )
		{
			Point< double , TDim > tCoordinate = s[0];
			for( unsigned int k=0 ; k<K ; k++ ) tCoordinate += ( s[k+1]-s[0] ) * p[k];

			Point< T , TDim , double > diff;
			for( unsigned d=0 ; d<TDim ; d++ ) diff[d] = signal.partial( d , tCoordinate );
			return t2s * diff;
		};
}

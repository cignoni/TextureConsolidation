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

template< unsigned int Dim >
SquareDistanceToTextureBoundary< Dim >::SquareDistanceToTextureBoundary
(
	const unsigned int res[Dim] ,
	bool verbose
)
	: verbose(verbose)
{
	for( unsigned int d=0 ; d<Dim ; d++ ) _res[d] = res[d];
}

template< unsigned int Dim >
template< SimplexProcessing::HasFunction< bool , typename RegularGrid< Dim >::Index > ActiveTexelFunctor >
RegularGrid< Dim , unsigned int > SquareDistanceToTextureBoundary< Dim >::operator()
(
	ActiveTexelFunctor && F
)
const
{
	using Index = typename RegularGrid< Dim >::Index;
	using Range = typename RegularGrid< Dim >::Range; 


	Miscellany::PerformanceMeter pMeter;
	RegularGrid< Dim , unsigned int > d2( _res );
	RegularGrid< Dim , unsigned char > raster( _res ) , boundary( _res );

	Range range;
	for( unsigned int d=0 ; d<Dim ; d++ ) range.second[d] = _res[d];

	// Identify all active texels
	{
		range.processParallel( [&]( unsigned int , Index I ){ raster( I ) = F( I ) ? 1 : 0; } );
	}
	if( verbose ) std::cout << pMeter( "Active texels" ) << std::endl;

	// Identify the boundary texels
	{
		auto Kernel = [&]( unsigned int , Index I )
			{
				if( raster( I ) ) boundary( I ) = 0;
				else
				{
					bool hasActiveNeighbors = false;
					Range::Intersect( range , Range(I).dilate(1) ).process( [&]( Index I ){ hasActiveNeighbors |= raster(I)!=0; } );
					boundary( I ) = hasActiveNeighbors ? 1 : 0;
				}
			};
		range.processParallel( Kernel );
	}
	if( verbose ) std::cout << pMeter( "Boundary texels" ) << std::endl;

	d2 = SquaredEDT< double , Dim >::Saito( boundary , verbose );
	if( verbose ) std::cout << pMeter( "Saito" ) << std::endl;

	ThreadPool::ParallelFor( 0 , d2.size() , [&]( size_t i ){ if( !raster[i] ) d2[i] = 0; } );

	return d2;
}

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

#ifndef TEXTURE_IO_INCLUDED
#define TEXTURE_IO_INCLUDED

#include <string>
#include <filesystem>
#include "Misha/Exceptions.h"
#include "Misha/Image.h"
#include "Misha/RegularGrid.h"
#include "Misha/Geometry.h"

namespace MishaK
{
	// Substitutes the index into a printf-style file-name pattern containing a single "%d" (e.g. "textures/%d.png")
	inline std::string IndexedFileName( const std::string &pattern , unsigned int index )
	{
		size_t pos = pattern.find( "%d" );
		if( pos==std::string::npos || pattern.find( '%' , pos+2 )!=std::string::npos || pattern.find( '%' )!=pos ) MK_ERROR_OUT( "File-name pattern must contain exactly one %d: " , pattern );
		return pattern.substr( 0 , pos ) + std::to_string( index ) + pattern.substr( pos+2 );
	}

	// Returns the number of consecutively indexed files, starting at index zero, that exist
	inline unsigned int IndexedFileCount( const std::string &pattern )
	{
		unsigned int count = 0;
		while( std::filesystem::exists( IndexedFileName( pattern , count ) ) ) count++;
		return count;
	}

	inline RegularGrid< 2 , Point< double , 3 > > ReadTexture( std::string fileName )
	{
		unsigned int res[2];
		unsigned char *pixels = ImageReader< 8 >::ReadColor( fileName , res[0] , res[1] );

		RegularGrid< 2 , Point< double , 3 > > img;
		img.resize( res );

		for( unsigned int j=0 , idx=0 ; j<img.res(1) ; j++ ) for( unsigned int i=0 ; i<img.res(0) ; i++ ) for( unsigned int k=0 ; k<3 ; k++ , idx++ )
			img(i,res[1]-1-j)[k] = pixels[idx]/255.;

		delete[] pixels;

		return img;
	};

	inline RegularGrid< 2 , double > ReadConfidenceTexture( std::string fileName )
	{
		RegularGrid< 2 , double > confidence;
		RegularGrid< 2 , Point< double , 3 > > _confidence = ReadTexture( fileName );
		confidence.resize( _confidence.res() );
		for( size_t i=0 ; i<_confidence.resolution() ; i++ ) confidence[i] = Point< double , 3 >::Dot( _confidence[i] , Point< double , 3 >(1./3,1./3,1./3) );
		return confidence;
	}

	inline void WriteTexture( std::string fileName , const RegularGrid< 2 , Point< double , 3 > > &img )
	{
		// Create the output folder if it does not exist yet
		std::filesystem::path folder = std::filesystem::path( fileName ).parent_path();
		if( !folder.empty() ) std::filesystem::create_directories( folder );

		unsigned char * pixels = new unsigned char[ img.resolution() * 3 ];
		for( unsigned int j=0 , idx=0 ; j<img.res(1) ; j++ ) for( unsigned int i=0 ; i<img.res(0) ; i++ ) for( unsigned int k=0 ; k<3 ; k++ , idx++ )
			pixels[idx] = static_cast< unsigned char >( std::min< double >( 255. , std::max< double >( 0. , img(i,img.res(1)-1-j)[k]*255. ) ) );

		ImageWriter< 8 >::Write( fileName , pixels , img.res(0) , img.res(1) , 3 );

		delete[] pixels;
	}

	inline void WriteConfidenceTexture( std::string fileName , const RegularGrid< 2 , double > &confidence )
	{
		RegularGrid< 2 , Point< double , 3 > > _confidence;
		_confidence.resize( confidence.res() );
		for( size_t i=0 ; i<confidence.resolution() ; i++ ) _confidence[i] = Point< double , 3 >(1.,1.,1.) * confidence[i];
		return WriteTexture( fileName , _confidence );
	}
}
#endif // TEXTURE_IO_INCLUDED

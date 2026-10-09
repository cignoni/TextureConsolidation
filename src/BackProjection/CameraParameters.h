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

#ifndef CAMERA_PARAMETERS_INCLUDED
#define CAMERA_PARAMETERS_INCLUDED

#include <string>
#include <vector>
#include "Misha/Exceptions.h"
#include "Misha/Geometry.h"
#include "Misha/CmdLineParser.h"

namespace MishaK
{
	// Pinhole camera read from a cameras.json file of the form
	//   { "frames": [ { "w": W, "h": H, "fx": FX, "fy": FY, "cx": CX, "cy": CY, "w2c": [[...],[...],[...],[...]], "file_path": "image.png" }, ... ] }
	// The parser is token-based: keys must appear in this order, and brackets/braces must be separated by white-space.
	// An optional trailing "blender_camera_location" entry per frame is accepted and ignored.
	struct CameraParameters
	{
		unsigned int res[2];
		XForm< double , 4 > w2c;
		XForm< double , 3 > k;
		std::string path;

		// Projects a world-space point to pixel coordinates
		Point< double , 2 > operator()( Point< double , 3 > p ) const
		{
			Point< double , 3 > _p = w2c(p);
			return k( Point< double , 2 >( _p[0]/_p[2] , _p[1]/_p[2] ) );
		}

		Point< double , 3 > position( void ) const { return w2c.inverse()( Point< double , 3 >() ); }

		static std::vector< CameraParameters > Read( std::string fileName )
		{
			std::vector< CameraParameters > cameraParameters;
			std::vector< std::string > words = ReadWords( fileName );
			if( !words.size() || words[0]!=std::string( "{" ) || words.back()!=std::string( "}" ) ) MK_ERROR_OUT( "Poorly formed file: " , fileName );
			words.erase( words.begin() , words.begin()+1 );
			words.pop_back();
			_ReadVector( std::string( "frames" ) , words );

			while( true )
			{
				bool hasMore;
				cameraParameters.push_back( _Read( words , hasMore ) );
				if( !hasMore ) break;
			}
			return cameraParameters;
		}

	protected:
		static std::string _StripKey( std::string word )
		{
			if( word.size()<3 ) MK_ERROR_OUT( "Key too short: " , word );
			if( word.back()!=':' ) MK_ERROR_OUT( "Key poorly formed: " , word );
			word.pop_back();
			if( word[0]!='\"' || word.back()!='\"' ) MK_ERROR_OUT( "Key poorly formed: " , word );
			return word.substr( 1 , word.size()-2 );
		};

		static bool _ReadVector( std::string key , std::vector< std::string > & words )
		{
			if( words.size()<3 ) MK_ERROR_OUT( "Insufficient words: " , words.size() );
			std::string _key = _StripKey( words[0] );
			if( _key!=key ) MK_ERROR_OUT( "Key mismatch: "  , key , " != " , _key );
			if( words[1]!=std::string( "[" ) ) MK_ERROR_OUT( "Failed to parse vector start: " , words[1] );
			if( words.back()!=std::string( "]" ) && words.back()!=std::string( "]," ) ) MK_ERROR_OUT( "Failed to parse vector end: " , words.back() );
			bool end = words.back()==std::string( "]" );
			words.erase( words.begin() , words.begin()+2 );
			words.pop_back();
			return !end;
		}

		template< typename T >
		static T _ReadScalar( std::string key , std::vector< std::string > &words , bool & hasMore )
		{
			T value;

			if( _StripKey( words[0] )!=key ) MK_ERROR_OUT( "Could not find key: " , key , " <- " , words[0] );
			words.erase( words.begin() , words.begin()+1 );

			if( words[0].back()==',' )
			{
				hasMore = true;
				words[0].pop_back();
			}
			else hasMore = false;

			try
			{
				if      constexpr( std::is_same_v< T , unsigned int > ) value = std::stoi( words[0] );
				else if constexpr( std::is_same_v< T , double > ) value = std::stod( words[0] );
				else MK_ERROR_OUT( "Unsupported type: " , typeid(T).name() );
			}
			catch( ... ){ MK_ERROR_OUT( "Failed to parse: " , key ); }
			words.erase( words.begin() , words.begin()+1 );

			return value;
		}

		template< unsigned int Dim >
		static SquareMatrix< double , Dim > _ReadMatrix( std::string key , std::vector< std::string > &words , bool & hasMore )
		{
			if( _StripKey( words[0] )!=key ) MK_ERROR_OUT( "Could not find key: " , key , " <- " , words[0] );
			words.erase( words.begin() , words.begin()+1 );

			SquareMatrix< double , Dim > value;

			if( words[0]!=std::string( "[" ) ) MK_ERROR_OUT( "Failed to parse matrix start: " , words[0] );
			words.erase( words.begin() , words.begin()+1 );

			for( unsigned int r=0 ; r<Dim ; r++ )
			{
				if( words[0]!=std::string( "[" ) ) MK_ERROR_OUT( "Failed to parse rows start: " , words[0] );
				words.erase( words.begin() , words.begin()+1 );

				for( unsigned int c=0 ; c<Dim ; c++ )
				{
					if( c+1<Dim )
					{
						if( words[0].back()!=',' ) MK_ERROR_OUT( "Expected more entries" );
						else words[0].pop_back();
					}
					try{ value(c,r) = std::stod( words[0] ); }
					catch( ... ){ MK_ERROR_OUT( "Failed to parse matrix entry: " , words[0] ); }
					words.erase( words.begin() , words.begin()+1 );
				}

				if( r+1<Dim ){ if( words[0]!=std::string( "]," ) ) MK_ERROR_OUT( "Failed to parse rows end: " , words[0] ); }
				else         { if( words[0]!=std::string( "]" ) ) MK_ERROR_OUT( "Failed to parse rows end: " , words[0] ); }
				words.erase( words.begin() , words.begin()+1 );
			}

			if( words[0]!=std::string( "]" ) && words[0]!=std::string( "]," ) ) MK_ERROR_OUT( "Failed to parse matrix end: " , words[0] );
			hasMore = words[0]==std::string( "]," );
			words.erase( words.begin() , words.begin()+1 );

			return value;
		}

		template< unsigned int Dim >
		static Point< double , Dim > _ReadPoint( std::string key , std::vector< std::string > &words , bool & hasMore )
		{
			if( _StripKey( words[0] )!=key ) MK_ERROR_OUT( "Could not find key: " , key , " <- " , words[0] );
			words.erase( words.begin() , words.begin()+1 );

			Point< double , Dim > value;

			if( words[0]!=std::string( "[" ) ) MK_ERROR_OUT( "Failed to parse point start: " , words[0] );
			words.erase( words.begin() , words.begin()+1 );

			for( unsigned int c=0 ; c<Dim ; c++ )
			{
				if( c+1<Dim )
				{
					if( words[0].back()!=',' ) MK_ERROR_OUT( "Expected more entries" );
					else words[0].pop_back();
				}
				try{ value[c] = std::stod( words[0] ); }
				catch( ... ){ MK_ERROR_OUT( "Failed to parse point entry: " , words[0] ); }
				words.erase( words.begin() , words.begin()+1 );
			}

			if( words[0]!=std::string( "]" ) && words[0]!=std::string( "]," ) ) MK_ERROR_OUT( "Failed to parse point end: " , words[0] );
			hasMore = words[0]==std::string( "]," );
			words.erase( words.begin() , words.begin()+1 );

			return value;
		}

		static std::string _ReadString( std::string key , std::vector< std::string > &words , bool & hasMore )
		{
			if( _StripKey( words[0] )!=key ) MK_ERROR_OUT( "Could not find key: " , key , " <- " , words[0] );
			words.erase( words.begin() , words.begin()+1 );

			if( words[0].back()==',' )
			{
				hasMore = true;
				words[0].pop_back();
			}
			else hasMore = false;

			if( words[0][0]!='\"' ||  words[0].back()!='\"' ) MK_ERROR_OUT( "String poorly formed: " , words[0] );
			std::string value = words[0].substr( 1 , words[0].size()-2 );
			words.erase( words.begin() , words.begin()+1 );

			return value;
		}

		static CameraParameters _Read( std::vector< std::string > & words , bool & hasMore )
		{
			CameraParameters param;
			bool _hasMore;

			if( words[0]!=std::string( "{" ) ) MK_ERROR_OUT( "Failed to parse parameters start: " , words[0] );
			words.erase( words.begin() , words.begin()+1 );

			param.res[0] = _ReadScalar< unsigned int >( "w" , words , _hasMore );
			if( !_hasMore ) MK_ERROR_OUT( "Expected more content" );

			param.res[1] = _ReadScalar< unsigned int >( "h" , words , _hasMore );
			if( !_hasMore ) MK_ERROR_OUT( "Expected more content" );

			param.k(0,0) = _ReadScalar< double >( "fx" , words , _hasMore );
			if( !_hasMore ) MK_ERROR_OUT( "Expected more content" );

			param.k(1,1) = _ReadScalar< double >( "fy" , words , _hasMore );
			if( !_hasMore ) MK_ERROR_OUT( "Expected more content" );

			param.k(2,0) = _ReadScalar< double >( "cx" , words , _hasMore );
			if( !_hasMore ) MK_ERROR_OUT( "Expected more content" );

			param.k(2,1) = _ReadScalar< double >( "cy" , words , _hasMore );
			if( !_hasMore ) MK_ERROR_OUT( "Expected more content" );

			param.k(2,2) = 1.;

			param.w2c = _ReadMatrix< 4 >( "w2c" , words , _hasMore );
			if( !_hasMore ) MK_ERROR_OUT( "Expected more content" );

			param.path = _ReadString( "file_path" , words , _hasMore );
			if( _hasMore )
			{
				_ReadPoint< 3 >( "blender_camera_location" , words , _hasMore );
				if( _hasMore ) MK_ERROR_OUT( "Did not expect more content" );
			}

			if( words[0]!=std::string( "}" ) && words[0]!=std::string( "}," ) ) MK_ERROR_OUT( "Failed to parse parameters end: " , words[0] );
			hasMore = words[0]==std::string( "}," );
			words.erase( words.begin() , words.begin()+1 );

			return param;
		}
	};
}
#endif // CAMERA_PARAMETERS_INCLUDED

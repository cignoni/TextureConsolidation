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
#include <fstream>
#include <sstream>
#include <cctype>
#include <cstdlib>
#include "Misha/Exceptions.h"
#include "Misha/Geometry.h"

namespace MishaK
{
	// Pinhole camera read from a cameras.json file of the form
	//   { "frames": [ { "w": W, "h": H, "fx": FX, "fy": FY, "cx": CX, "cy": CY, "w2c": [[...],[...],[...],[...]], "file_path": "image.png" }, ... ] }
	// Keys may appear in any order, and other keys are ignored.
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
			std::ifstream file( fileName );
			if( !file ) MK_ERROR_OUT( "Could not open cameras file: " , fileName );
			std::stringstream stream;
			stream << file.rdbuf();
			std::string text = stream.str();

			size_t pos = 0;
			_JSON root = _JSON::Parse( text , pos );
			_JSON::SkipWhiteSpace( text , pos );
			if( pos!=text.size() ) MK_ERROR_OUT( "Unexpected content after the JSON object in: " , fileName );

			const _JSON &frames = root.get( "frames" );
			if( frames.type!=_JSON::ARRAY ) MK_ERROR_OUT( "Expected \"frames\" to be an array" );

			std::vector< CameraParameters > cameraParameters( frames.values.size() );
			for( size_t i=0 ; i<frames.values.size() ; i++ )
			{
				const _JSON &frame = frames.values[i];
				CameraParameters &param = cameraParameters[i];

				param.res[0] = static_cast< unsigned int >( frame.get( "w" ).toNumber( "w" ) );
				param.res[1] = static_cast< unsigned int >( frame.get( "h" ).toNumber( "h" ) );
				param.k(0,0) = frame.get( "fx" ).toNumber( "fx" );
				param.k(1,1) = frame.get( "fy" ).toNumber( "fy" );
				param.k(2,0) = frame.get( "cx" ).toNumber( "cx" );
				param.k(2,1) = frame.get( "cy" ).toNumber( "cy" );
				param.k(2,2) = 1.;

				const _JSON &w2c = frame.get( "w2c" );
				if( w2c.type!=_JSON::ARRAY || w2c.values.size()!=4 ) MK_ERROR_OUT( "Expected \"w2c\" to be a 4x4 matrix" );
				for( unsigned int r=0 ; r<4 ; r++ )
				{
					if( w2c.values[r].type!=_JSON::ARRAY || w2c.values[r].values.size()!=4 ) MK_ERROR_OUT( "Expected \"w2c\" to be a 4x4 matrix" );
					for( unsigned int c=0 ; c<4 ; c++ ) param.w2c(c,r) = w2c.values[r].values[c].toNumber( "w2c" );
				}

				const _JSON &path = frame.get( "file_path" );
				if( path.type!=_JSON::STRING ) MK_ERROR_OUT( "Expected \"file_path\" to be a string" );
				param.path = path.string;
			}
			return cameraParameters;
		}

	protected:
		// A minimal JSON reader, sufficient for cameras.json
		struct _JSON
		{
			enum Type { NUL , BOOLEAN , NUMBER , STRING , ARRAY , OBJECT };
			Type type = NUL;
			double number = 0;
			std::string string;
			std::vector< std::string > keys;	// For objects
			std::vector< _JSON > values;		// For arrays and objects

			const _JSON &get( const std::string &key ) const
			{
				if( type!=OBJECT ) MK_ERROR_OUT( "Expected an object containing: " , key );
				for( size_t i=0 ; i<keys.size() ; i++ ) if( keys[i]==key ) return values[i];
				MK_ERROR_OUT( "Missing key: " , key );
				return *this;
			}

			double toNumber( const std::string &key ) const
			{
				if( type!=NUMBER ) MK_ERROR_OUT( "Expected a number for: " , key );
				return number;
			}

			static void SkipWhiteSpace( const std::string &text , size_t &pos ){ while( pos<text.size() && std::isspace( static_cast< unsigned char >( text[pos] ) ) ) pos++; }

			static _JSON Parse( const std::string &text , size_t &pos )
			{
				_JSON value;
				SkipWhiteSpace( text , pos );
				if( pos>=text.size() ) MK_ERROR_OUT( "Unexpected end of JSON" );
				char c = text[pos];
				if( c=='{' || c=='[' )
				{
					value.type = c=='{' ? OBJECT : ARRAY;
					char close = c=='{' ? '}' : ']';
					pos++;
					SkipWhiteSpace( text , pos );
					if( pos<text.size() && text[pos]==close ){ pos++ ; return value; }
					while( true )
					{
						if( value.type==OBJECT )
						{
							SkipWhiteSpace( text , pos );
							value.keys.push_back( _ParseString( text , pos ) );
							SkipWhiteSpace( text , pos );
							if( pos>=text.size() || text[pos]!=':' ) MK_ERROR_OUT( "Expected ':' after key: " , value.keys.back() );
							pos++;
						}
						value.values.push_back( Parse( text , pos ) );
						SkipWhiteSpace( text , pos );
						if( pos>=text.size() ) MK_ERROR_OUT( "Unexpected end of JSON" );
						if( text[pos]==',' ){ pos++ ; continue; }
						if( text[pos]==close ){ pos++ ; break; }
						MK_ERROR_OUT( "Expected ',' or '" , close , "' at position " , pos );
					}
				}
				else if( c=='"' ) value.type = STRING , value.string = _ParseString( text , pos );
				else if( text.compare( pos , 4 , "true" )==0 ) value.type = BOOLEAN , value.number = 1 , pos += 4;
				else if( text.compare( pos , 5 , "false" )==0 ) value.type = BOOLEAN , value.number = 0 , pos += 5;
				else if( text.compare( pos , 4 , "null" )==0 ) value.type = NUL , pos += 4;
				else
				{
					const char *start = text.c_str() + pos;
					char *end;
					value.type = NUMBER;
					value.number = std::strtod( start , &end );
					if( end==start ) MK_ERROR_OUT( "Could not parse JSON value at position " , pos );
					pos += end - start;
				}
				return value;
			}

		protected:
			static std::string _ParseString( const std::string &text , size_t &pos )
			{
				if( pos>=text.size() || text[pos]!='"' ) MK_ERROR_OUT( "Expected a string at position " , pos );
				pos++;
				std::string str;
				while( pos<text.size() && text[pos]!='"' )
				{
					if( text[pos]=='\\' && pos+1<text.size() )
					{
						pos++;
						switch( text[pos] )
						{
						case 'n': str += '\n' ; break;
						case 't': str += '\t' ; break;
						case 'r': str += '\r' ; break;
						case 'b': str += '\b' ; break;
						case 'f': str += '\f' ; break;
						case 'u': MK_ERROR_OUT( "Unicode escapes are not supported" );
						default: str += text[pos];
						}
					}
					else str += text[pos];
					pos++;
				}
				if( pos>=text.size() ) MK_ERROR_OUT( "Unterminated string" );
				pos++;
				return str;
			}
		};
	};
}
#endif // CAMERA_PARAMETERS_INCLUDED

/*
    This file is part of the Tetrahedral Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The Tetrahedral Trees library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The Tetrahedral Trees library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Tetrahedral Trees library.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "string_management.h"

namespace string_management
{

string get_file_name(string path)
{
    // Strip the path, isolate file name
    string filename = strip_path(path);

    // remove the extension if any
    size_t point_index = filename.find_last_of('.');
    if (point_index == string::npos)
        return filename;
    else
        return filename.substr(0, point_index);
}

string strip_path(string path)
{
    string filename;
    size_t name_index_ms = path.find_last_of('\\');
    size_t name_index_unix = path.find_last_of('/');
    if (name_index_ms != string::npos)
        filename = path.substr(name_index_ms + 1);
    else if (name_index_unix != string::npos)
        filename = path.substr(name_index_unix + 1);
    else
        filename = path;
    return filename;
}

void tokenize(const string& str, vector<string>& tokens, const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

}

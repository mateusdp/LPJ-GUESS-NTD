////////////////////////////////////////////////////////////////////////////////
/// \file guessstring.h
/// \brief Utility functions for working with strings (std::string and char*)
///
/// \author Joe Siltberg
/// $Date: 2022-11-22 12:55:59 +0100 (Tue, 22 Nov 2022) $
///
////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_GUESSSTRING_H
#define LPJ_GUESS_GUESSSTRING_H

#include <string>
#include <string.h>

/// Removes leading and trailing whitespace from a string
std::string trim(const std::string& str);

/// Converts a string to upper case
std::string to_upper(const std::string& str);

/// Converts a string to lower case
std::string to_lower(const std::string& str);

/// Creates a string with printf style formatting
std::string format_string(const char* format, ...);

/// Help function that splits string into "words"
int split_string(char* str);

/// Help function that finds substring in string
bool issubstring(const char* string, const char* substring);

#endif // LPJ_GUESS_GUESSSTRING_H

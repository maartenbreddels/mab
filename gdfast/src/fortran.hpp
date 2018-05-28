#pragma once
#include <stdlib.h>
#include <inttypes.h>

template<typename T>
void fileread(FILE* file, T* value)
{
	fread(value, sizeof T, 1, file);
}

template<typename T, typename... Tail>
void read_fortran(FILE* file, T* var, Tail... tail)
{
	uint32_t struct_size1, struct_size2 ;
	fileread(file, &struct_size1);
	_read_fortran(file, var, tail...);
	fileread(file, &struct_size2);
	if(struct_size1 != struct_size2) {
		fprintf(stderr, "error reading fortran block: size before and after say: %d %d\n", struct_size1, struct_size2)
		exit(1);
	}
}

template<typename T, typename... Tail>
void _read_fortran(FILE* file, T* var, Tail... tail)
{
	fileread(file, var);
	_read_fortran(tail...);
}

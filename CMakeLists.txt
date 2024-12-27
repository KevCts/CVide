cmake_minimum_required(VERSION 3.30)

project(libcvide)

set(SRCS
    coomat.c
    )
    
set(HEADERS
    coomat.h
    )

add_library(libcvide ${SRCS} ${HEADERS})

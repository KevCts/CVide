cmake_minimum_required(VERSION 3.30)

project(libcvide)

set(SRCS
    coomat.c
    csrmat.c
    solve.c
    )
    
set(HEADERS
    coomat.h
    csrmat.h
    solve.h
    )

add_library(libcvide ${SRCS} ${HEADERS})

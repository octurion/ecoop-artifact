cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(paper_charts)

find_package(LATEX COMPONENTS PDFLATEX)

if (LATEX_PDFLATEX_FOUND)
    add_custom_target(charts_pdf
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}
        COMMAND ${PDFLATEX_COMPILER} -interaction=nonstopmode ${PROJECT_BINARY_DIR}/charts.tex
        WORKING_DIR ${PROJECT_BINARY_DIR})
else()
    message(WARNING "No installation of pdflatex found, it will not be possible to generate the charts")
endif()


set(BINDINGS_SRC
    ${CMAKE_SOURCE_DIR}/src/bindings.cpp
    ${CMAKE_SOURCE_DIR}/src/pybind_OpenLPT/pyMatrix.cpp
)

# Create pybind11 module
pybind11_add_module(pyOpenLPT ${BINDINGS_SRC})

# Link the pybind11 module with the imported libraries
target_link_libraries(pyOpenLPT PRIVATE Matrix)
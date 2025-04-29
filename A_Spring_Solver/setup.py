from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np
import os

ext = Extension(
    name="spring_solver",
    sources=["spring_solver.pyx"],
    include_dirs=[np.get_include()],
    language="c",
    extra_compile_args=["/O2","/MD"],
    export_symbols=["solve_spring_mass_c"],  # this is key for the .lib/.dll
    extra_link_args=["/DEF:spring_solver.def"] if os.name == 'nt' else [],
)

setup(
    name="spring_solver",
    ext_modules=cythonize([ext], compiler_directives={'language_level': "3"}),
)
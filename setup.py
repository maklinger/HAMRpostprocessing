from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from Cython.Build import cythonize
import numpy as np

ext_modules = cythonize(
    Extension(
        name="pp_c2",  # can be used as `import pp_c`
        sources=["c_ext/pp_c2.pyx", "c_ext/functions2.c"],
        include_dirs=[np.get_include()],
        extra_compile_args=["-fopenmp"],
        extra_link_args=["-fopenmp"]
    )
)

setup(
    name="HAMRpostprocessing",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    ext_modules=ext_modules,
    include_dirs=[np.get_include()],
)

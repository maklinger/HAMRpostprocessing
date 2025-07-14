import os
from setuptools import setup
from setuptools.extension import Extension
from Cython.Distutils import build_ext
import numpy as np
import sys

def build_pp_c():
    print("Building pp_c...")
    setup(
        script_args=["build_ext", "--inplace"],
        cmdclass={'build_ext': build_ext},
        ext_modules=[Extension(
            "pp_c",
            sources=["c_ext/pp_c.pyx", "c_ext/functions.c"],
            include_dirs=[np.get_include()],
            extra_compile_args=["-fopenmp"],
            extra_link_args=["-O2", "-fopenmp"]
        )]
    )
    print("Done building pp_c!")

if __name__ == "__main__":
    build_pp_c()

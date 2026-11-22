"""Native extension configuration; normal installation uses pyproject.toml."""

import os
from pathlib import Path

from setuptools import Extension, setup

DEFAULT_EIGEN_INCLUDE = "/usr/include/eigen3"
EIGEN_INCLUDE = os.environ.get("EIGEN_INCLUDE", DEFAULT_EIGEN_INCLUDE)
HEADER_FILES = sorted(str(path) for path in Path("include/langevin").glob("*.hpp"))
COMPILE_FLAGS = ["-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic"]

setup(
    license="MIT",
    license_files=["LICENSE"],
    ext_modules=[
        Extension(
            "langevin_monte_carlo._core",
            sources=["python/bindings.cpp"],
            include_dirs=["include", EIGEN_INCLUDE],
            depends=HEADER_FILES,
            language="c++",
            extra_compile_args=COMPILE_FLAGS,
        )
    ]
)

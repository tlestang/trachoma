from setuptools import Extension, setup

setup(
    ext_modules=[
        Extension(
            name="ntd-trachoma.libtrachoma",
            sources=["src/libtrachoma/trachoma.c"],
        ),
    ]
)
